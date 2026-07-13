"""Reusable equilibrium contexts backed by the GLISS C ABI."""

import ctypes
import operator
import os
from pathlib import Path
from typing import Any, Tuple, Union

import numpy as np

from . import _load_library, _require_symbols

_C_INT_MAX = 2 ** (ctypes.sizeof(ctypes.c_int) * 8 - 1) - 1
_ERROR_CAPACITY = 512
PathLike = Union[str, os.PathLike]


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


def _bind(library: Any) -> None:
    _require_symbols(
        library,
        (
            "gliss_equilibrium_create",
            "gliss_equilibrium_destroy",
            "gliss_equilibrium_surface_count",
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

    @property
    def closed(self) -> bool:
        """Whether the native equilibrium has been released."""
        return self._handle.value is None

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
        return f"<gliss.Equilibrium(path={self.path!r}, state={state!r})>"

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

    def _require_open(self) -> None:
        if self.closed:
            raise RuntimeError("Equilibrium is closed")
