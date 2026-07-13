"""Reusable equilibrium contexts backed by the GLISS C ABI."""

import ctypes
import hashlib
import operator
import os
import tempfile
from pathlib import Path
from typing import Any, Optional, Tuple, Union

import numpy as np

from . import _load_library, _require_symbols

_C_INT_MAX = 2 ** (ctypes.sizeof(ctypes.c_int) * 8 - 1) - 1
_ERROR_CAPACITY = 512
PathLike = Union[str, os.PathLike]
_FileIdentity = Tuple[int, int, int, int, int]


class GlissError(RuntimeError):
    """Base class for failures reported by the GLISS native library."""


class GlissIOError(GlissError):
    """An equilibrium or result file could not be read or written."""


class GlissComputationError(GlissError):
    """A native numerical operation failed."""


class GlissCapacityError(GlissError):
    """A caller-owned output buffer is too small."""


class GlissArgumentError(GlissError):
    """The native ABI rejected an argument."""


class GlissAllocationError(GlissError):
    """The native library could not allocate an internal object."""


class GlissInternalError(GlissError):
    """The native library encountered an internal lifecycle failure."""


_STATUS_EXCEPTIONS = {
    1: GlissIOError,
    2: GlissComputationError,
    3: GlissCapacityError,
    4: GlissArgumentError,
    5: GlissAllocationError,
    6: GlissInternalError,
}


def _file_identity(path: Path) -> _FileIdentity:
    metadata = path.stat()
    return (
        metadata.st_dev,
        metadata.st_ino,
        metadata.st_size,
        metadata.st_mtime_ns,
        metadata.st_ctime_ns,
    )


def _stable_file_digest(
    path: Path, expected_identity: Optional[_FileIdentity] = None
) -> Tuple[int, str]:
    try:
        before = _file_identity(path)
        if expected_identity is not None and before != expected_identity:
            raise GlissIOError("equilibrium export changed after loading")
        digest = hashlib.sha256()
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
        if _file_identity(path) != before:
            raise GlissIOError("equilibrium export changed while checksumming")
    except OSError as error:
        message = "equilibrium export changed while checksumming"
        if expected_identity is not None:
            message = "equilibrium export changed after loading"
        raise GlissIOError(message) from error
    return before[2], digest.hexdigest()


def _resolution(value: Any, name: str) -> int:
    if isinstance(value, (bool, np.bool_)):
        raise TypeError(f"{name} must be an integer")
    try:
        result = operator.index(value)
    except TypeError as error:
        raise TypeError(f"{name} must be an integer") from error
    if not 1 <= result <= _C_INT_MAX:
        raise ValueError(f"{name} must be between 1 and {_C_INT_MAX}")
    return result


def _export_path(path: PathLike) -> Tuple[Path, bytes]:
    try:
        value = os.fspath(path)
    except TypeError as error:
        raise TypeError("path must be a string or path-like object") from error
    if not isinstance(value, str):
        raise TypeError("path must resolve to a string")
    if "\0" in value:
        raise ValueError("path contains a null byte")
    export = Path(value)
    if not export.exists():
        raise FileNotFoundError(f"equilibrium export does not exist: {export}")
    if not export.is_file():
        raise ValueError(f"equilibrium export is not a file: {export}")
    encoded = os.fsencode(export)
    if len(encoded) > ctypes.c_size_t(-1).value:
        raise ValueError("encoded path is too long for the GLISS C API")
    return export, encoded


def _output_path(path: PathLike) -> Path:
    try:
        value = os.fspath(path)
    except TypeError as error:
        raise TypeError("output path must be a string or path-like object") from error
    if not isinstance(value, str):
        raise TypeError("output path must resolve to a string")
    if "\0" in value:
        raise ValueError("output path contains a null byte")
    output = Path(value)
    if not output.parent.is_dir():
        raise FileNotFoundError(f"output directory does not exist: {output.parent}")
    if output.exists() and not output.is_file():
        raise ValueError(f"output path is a directory or special file: {output}")
    return output


def _bind(library: Any) -> None:
    _require_symbols(
        library,
        (
            "gliss_equilibrium_create",
            "gliss_equilibrium_destroy",
            "gliss_equilibrium_surface_count",
            "gliss_equilibrium_schema_version",
            "gliss_equilibrium_write",
            "gliss_mercier_profile_context",
        ),
        "equilibrium context",
    )
    library.gliss_equilibrium_create.argtypes = (
        ctypes.c_char_p,
        ctypes.c_size_t,
        ctypes.POINTER(ctypes.c_void_p),
        ctypes.c_void_p,
        ctypes.c_size_t,
    )
    library.gliss_equilibrium_create.restype = ctypes.c_int
    library.gliss_equilibrium_destroy.argtypes = (
        ctypes.POINTER(ctypes.c_void_p),
        ctypes.c_void_p,
        ctypes.c_size_t,
    )
    library.gliss_equilibrium_destroy.restype = ctypes.c_int
    library.gliss_equilibrium_surface_count.argtypes = (
        ctypes.c_void_p,
        ctypes.POINTER(ctypes.c_size_t),
        ctypes.c_void_p,
        ctypes.c_size_t,
    )
    library.gliss_equilibrium_surface_count.restype = ctypes.c_int
    library.gliss_equilibrium_schema_version.argtypes = (
        ctypes.c_void_p,
        ctypes.POINTER(ctypes.c_int32),
        ctypes.c_void_p,
        ctypes.c_size_t,
    )
    library.gliss_equilibrium_schema_version.restype = ctypes.c_int
    library.gliss_equilibrium_write.argtypes = (
        ctypes.c_void_p,
        ctypes.c_char_p,
        ctypes.c_size_t,
        ctypes.c_void_p,
        ctypes.c_size_t,
    )
    library.gliss_equilibrium_write.restype = ctypes.c_int
    library.gliss_mercier_profile_context.argtypes = (
        ctypes.c_void_p,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_size_t,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_size_t),
        ctypes.c_void_p,
        ctypes.c_size_t,
    )
    library.gliss_mercier_profile_context.restype = ctypes.c_int


def _error_buffer() -> ctypes.Array:
    return ctypes.create_string_buffer(_ERROR_CAPACITY)


def _raise_for_status(status: int, error: ctypes.Array, operation: str) -> None:
    if status == 0:
        return
    message = error.value.decode("utf-8", errors="replace") or "no error detail"
    exception = _STATUS_EXCEPTIONS.get(status, GlissError)
    raise exception(f"{operation} failed with status {status}: {message}")


class Equilibrium:
    """Loaded GVEC/CAS3D equilibrium with explicit native lifetime."""

    def __init__(self, path: PathLike):
        self.path, encoded = _export_path(path)
        source_identity = _file_identity(self.path)
        self._library = _load_library()
        _bind(self._library)
        self._handle = ctypes.c_void_p()
        error = _error_buffer()
        status = self._library.gliss_equilibrium_create(
            encoded,
            len(encoded),
            ctypes.byref(self._handle),
            error,
            len(error),
        )
        if status != 0 and self._handle.value is not None:
            self._library.gliss_equilibrium_destroy(ctypes.byref(self._handle), None, 0)
        _raise_for_status(status, error, "gliss_equilibrium_create")
        if self._handle.value is None:
            raise GlissInternalError("GLISS returned a null equilibrium handle")
        try:
            source_changed = _file_identity(self.path) != source_identity
        except OSError as error:
            self._library.gliss_equilibrium_destroy(ctypes.byref(self._handle), None, 0)
            raise GlissIOError("equilibrium export changed while loading") from error
        if source_changed:
            self._library.gliss_equilibrium_destroy(ctypes.byref(self._handle), None, 0)
            raise GlissIOError("equilibrium export changed while loading")
        self._source_identity = source_identity

    @property
    def closed(self) -> bool:
        """Whether the native equilibrium has been released."""
        return self._handle.value is None

    @property
    def schema_version(self) -> int:
        """NetCDF export schema: 0 for legacy inputs and 1 for GLISS exports."""
        self._require_open()
        version = ctypes.c_int32()
        error = _error_buffer()
        status = self._library.gliss_equilibrium_schema_version(
            self._handle, ctypes.byref(version), error, len(error)
        )
        _raise_for_status(status, error, "gliss_equilibrium_schema_version")
        if version.value not in (0, 1):
            raise GlissInternalError(
                f"GLISS returned unsupported equilibrium schema {version.value}"
            )
        return version.value

    def close(self) -> None:
        """Release the native equilibrium; repeated calls are safe."""
        if self.closed:
            return
        error = _error_buffer()
        status = self._library.gliss_equilibrium_destroy(
            ctypes.byref(self._handle), error, len(error)
        )
        _raise_for_status(status, error, "gliss_equilibrium_destroy")
        if not self.closed:
            raise GlissInternalError("GLISS did not clear the equilibrium handle")

    def __enter__(self) -> "Equilibrium":
        if self.closed:
            raise RuntimeError("Equilibrium is closed")
        return self

    def __exit__(self, exc_type: Any, exc_value: Any, traceback: Any) -> None:
        self.close()

    def __repr__(self) -> str:
        state = "closed" if self.closed else "open"
        return f"<gliss.Equilibrium(path={self.path.name!r}, state={state!r})>"

    def _surface_count(self) -> int:
        self._require_open()
        surfaces = ctypes.c_size_t()
        error = _error_buffer()
        status = self._library.gliss_equilibrium_surface_count(
            self._handle, ctypes.byref(surfaces), error, len(error)
        )
        _raise_for_status(status, error, "gliss_equilibrium_surface_count")
        if surfaces.value > np.iinfo(np.intp).max:
            raise GlissCapacityError(
                f"surface count {surfaces.value} exceeds NumPy's index limit"
            )
        return surfaces.value

    def mercier_profile(
        self, n_theta: int = 64, n_zeta: int = 64
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Return radial coordinates and Mercier discriminants."""
        self._require_open()
        n_theta = _resolution(n_theta, "n_theta")
        n_zeta = _resolution(n_zeta, "n_zeta")
        count = self._surface_count()
        s_values = np.empty(count, dtype=np.float64)
        d_mercier = np.empty(count, dtype=np.float64)
        written = ctypes.c_size_t()
        error = _error_buffer()
        status = self._library.gliss_mercier_profile_context(
            self._handle,
            n_theta,
            n_zeta,
            count,
            s_values.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            d_mercier.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            ctypes.byref(written),
            error,
            len(error),
        )
        _raise_for_status(status, error, "gliss_mercier_profile_context")
        if written.value != count:
            raise GlissCapacityError(
                f"GLISS wrote {written.value} surfaces; expected {count}"
            )
        return s_values, d_mercier

    def write(self, path: PathLike) -> Path:
        """Atomically write a version-1 equilibrium export."""
        self._require_open()
        destination = _output_path(path)
        descriptor, temporary_name = tempfile.mkstemp(
            dir=destination.parent,
            prefix=f".{destination.name}.",
            suffix=".tmp",
        )
        os.close(descriptor)
        temporary = Path(temporary_name)
        temporary.unlink()
        try:
            encoded = os.fsencode(temporary)
            if len(encoded) > ctypes.c_size_t(-1).value:
                raise ValueError("encoded output path is too long for the GLISS C API")
            error = _error_buffer()
            status = self._library.gliss_equilibrium_write(
                self._handle, encoded, len(encoded), error, len(error)
            )
            _raise_for_status(status, error, "gliss_equilibrium_write")
            os.replace(temporary, destination)
            return destination
        finally:
            if temporary.exists():
                temporary.unlink()

    def _require_open(self) -> None:
        if self.closed:
            raise RuntimeError("Equilibrium is closed")
