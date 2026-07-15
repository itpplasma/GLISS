"""Deterministic binary interchange for complete stability spectra."""

from __future__ import annotations

import json
import math
import os
import tempfile
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Tuple

import numpy as np

from ._result_schema import stability_result_from_dict, stability_result_to_dict
from ._schema_support import _unique_object, document_path, fields, schema
from .equilibrium import GlissAllocationError, GlissInternalError, PathLike
from .full_spectrum import (
    FullSpectrumResult,
    FullStabilityResult,
    _validate_full_spectrum,
)
from .schema import RunManifest, StabilityConfiguration, _create_run_manifest
from .stability import SpectrumResult, StabilityResult

FULL_RESULT_SCHEMA = "gliss.stability.full-result"
FULL_RUN_SCHEMA = "gliss.stability.full-run"
_METADATA_NAME = "metadata.json"
_ARRAY_ATTRIBUTES = tuple(
    "eigenvalues rayleigh_quotients residuals resolutions eigenvectors".split()
)
_STORAGE = dict(
    format="zip-npy",
    array_dtype="<f8",
    array_order="eigenpair-component",
    compression="stored",
)
_MAX_METADATA_BYTES = 16 * 1024 * 1024


def _certified_result(result: FullStabilityResult) -> StabilityResult:
    return StabilityResult(
        (result.classes[0].certified_lowest, result.classes[1].certified_lowest)
    )


def _unknown_count(certified: SpectrumResult, context: str) -> int:
    count = certified.normal_unknowns + certified.eta_unknowns + certified.mu_unknowns
    if count > math.isqrt(np.iinfo(np.intp).max):
        raise ValueError(
            f"{context} unknown count {count} exceeds NumPy's full-spectrum limit"
        )
    return count


def _expected_shapes(item: FullSpectrumResult) -> Dict[str, Tuple[int, ...]]:
    count = _unknown_count(
        item.certified_lowest, f"class {item.certified_lowest.parity_class}"
    )
    shapes: Dict[str, Tuple[int, ...]] = {name: (count,) for name in _ARRAY_ATTRIBUTES}
    shapes["eigenvectors"] = (count, count)
    return shapes


def _validate_full_result(result: FullStabilityResult) -> StabilityResult:
    if not isinstance(result, FullStabilityResult):
        raise TypeError("result must be a gliss.FullStabilityResult")
    if not isinstance(result.classes, tuple) or len(result.classes) != 2:
        raise ValueError("full result classes must contain parity classes 1 then 2")
    if any(not isinstance(item, FullSpectrumResult) for item in result.classes):
        raise TypeError("full result classes must be gliss.FullSpectrumResult objects")
    certified = _certified_result(result)
    stability_result_to_dict(certified)
    if tuple(item.certified_lowest.parity_class for item in result.classes) != (1, 2):
        raise ValueError("full result parity classes must be 1 then 2")
    for index, item in enumerate(result.classes, start=1):
        shapes = _expected_shapes(item)
        arrays = []
        for name in _ARRAY_ATTRIBUTES:
            array = getattr(item, name)
            if not isinstance(array, np.ndarray):
                raise TypeError(f"class {index} {name} must be a NumPy array")
            if array.dtype != np.dtype(np.float64):
                raise ValueError(f"class {index} {name} dtype must be float64")
            if array.shape != shapes[name]:
                raise ValueError(
                    f"class {index} {name} shape is {array.shape}; "
                    f"expected {shapes[name]}"
                )
            arrays.append(array)
        try:
            _validate_full_spectrum(item.certified_lowest, *arrays)
        except GlissInternalError as error:
            raise ValueError(
                f"class {index} full spectrum is inconsistent: {error}"
            ) from error
    return certified


def _full_result_metadata(result: FullStabilityResult) -> Dict[str, Any]:
    certified = _validate_full_result(result)
    certified_document = stability_result_to_dict(certified)
    return {
        "schema": FULL_RESULT_SCHEMA,
        "schema_version": certified_document["schema_version"],
        "storage": dict(_STORAGE),
        "certified_result": certified_document,
    }


def _metadata_bytes(document: Mapping[str, Any]) -> bytes:
    try:
        text = json.dumps(
            document, allow_nan=False, indent=2, sort_keys=True, ensure_ascii=False
        )
    except (TypeError, ValueError) as error:
        raise ValueError(f"metadata is not finite JSON data: {error}") from error
    return (text + "\n").encode("utf-8")


def _entry_name(parity_class: int, attribute: str) -> str:
    return f"class-{parity_class}-{attribute.replace('_', '-')}.npy"


def _zip_info(name: str) -> zipfile.ZipInfo:
    info = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
    info.compress_type = zipfile.ZIP_STORED
    info.create_system = 3
    info.external_attr = 0o600 << 16
    return info


def _storage_array(array: np.ndarray, context: str) -> np.ndarray:
    try:
        return np.ascontiguousarray(array, dtype=np.dtype("<f8"))
    except MemoryError as error:
        raise GlissAllocationError(f"unable to store {context}") from error


def _write_archive(
    destination: Path, metadata: Mapping[str, Any], result: FullStabilityResult
) -> None:
    payload = _metadata_bytes(metadata)
    if len(payload) > _MAX_METADATA_BYTES:
        unit = "byte" if _MAX_METADATA_BYTES == 1 else "bytes"
        raise ValueError(f"metadata exceeds {_MAX_METADATA_BYTES} {unit}")
    descriptor = None
    temporary = None
    try:
        descriptor, name = tempfile.mkstemp(
            dir=destination.parent,
            prefix=f".{destination.name}.",
            suffix=".tmp",
        )
        temporary = Path(name)
        os.close(descriptor)
        descriptor = None
        with zipfile.ZipFile(temporary, "w", compression=zipfile.ZIP_STORED) as archive:
            archive.writestr(_zip_info(_METADATA_NAME), payload)
            for parity_class, item in enumerate(result.classes, start=1):
                for attribute in _ARRAY_ATTRIBUTES:
                    array = _storage_array(
                        getattr(item, attribute),
                        f"class {parity_class} {attribute}",
                    )
                    info = _zip_info(_entry_name(parity_class, attribute))
                    with archive.open(info, "w") as stream:
                        np.lib.format.write_array(stream, array, allow_pickle=False)
        with temporary.open("rb") as output:
            os.fsync(output.fileno())
        os.replace(temporary, destination)
    finally:
        if descriptor is not None:
            os.close(descriptor)
        if temporary is not None and temporary.exists():
            temporary.unlink()


def _write_container(
    path: PathLike, metadata: Mapping[str, Any], result: FullStabilityResult
) -> None:
    _validate_full_result(result)
    destination = document_path(path, "output")
    if not destination.parent.is_dir():
        raise FileNotFoundError(
            f"output directory does not exist: {destination.parent}"
        )
    _write_archive(destination, metadata, result)


def _open_archive(path: PathLike) -> zipfile.ZipFile:
    source = document_path(path, "input")
    if not source.exists():
        raise FileNotFoundError(f"input file does not exist: {source}")
    if not source.is_file():
        raise ValueError(f"input path is not a file: {source}")
    try:
        archive = zipfile.ZipFile(source, "r")
    except (OSError, zipfile.BadZipFile) as error:
        raise ValueError(
            f"{source}: invalid or truncated full-spectrum container"
        ) from error
    try:
        _validate_archive_entries(archive, source)
    except Exception:
        archive.close()
        raise
    return archive


def _validate_archive_entries(archive: zipfile.ZipFile, source: Path) -> None:
    information = archive.infolist()
    names = [item.filename for item in information]
    if len(names) != len(set(names)):
        raise ValueError(f"{source}: full-spectrum container has duplicate entries")
    expected = {_METADATA_NAME} | {
        _entry_name(parity_class, attribute)
        for parity_class in (1, 2)
        for attribute in _ARRAY_ATTRIBUTES
    }
    unknown = sorted(set(names) - expected)
    if unknown:
        raise ValueError(
            f"{source}: full-spectrum container has unknown entry {unknown[0]!r}"
        )
    missing = sorted(expected - set(names))
    if missing:
        raise ValueError(
            f"{source}: full-spectrum container is missing entry {missing[0]!r}"
        )
    if archive.comment:
        raise ValueError(
            f"{source}: full-spectrum container has an unsupported comment"
        )
    for item in information:
        if item.compress_type != zipfile.ZIP_STORED:
            raise ValueError(f"{source}: entry {item.filename!r} must be stored")
        if item.flag_bits & 1:
            raise ValueError(f"{source}: encrypted entries are not supported")


def _read_metadata(archive: zipfile.ZipFile, source: Path) -> Dict[str, Any]:
    info = archive.getinfo(_METADATA_NAME)
    if info.file_size > _MAX_METADATA_BYTES:
        raise ValueError(f"{source}: metadata exceeds {_MAX_METADATA_BYTES} bytes")
    try:
        text = archive.read(info).decode("utf-8")
    except UnicodeDecodeError as error:
        raise ValueError(f"{source}: metadata is not valid UTF-8") from error
    except (OSError, zipfile.BadZipFile) as error:
        raise ValueError(
            f"{source}: invalid or truncated full-spectrum container"
        ) from error
    try:
        document = json.loads(text, object_pairs_hook=_unique_object)
    except json.JSONDecodeError as error:
        raise ValueError(
            f"{source}: invalid or truncated JSON at line {error.lineno}, "
            f"column {error.colno}"
        ) from error
    except ValueError as error:
        raise ValueError(f"{source}: {error}") from error
    if not isinstance(document, dict):
        raise ValueError(f"{source}: metadata must be an object")
    return document


def _parse_result_metadata(document: Any) -> StabilityResult:
    value = fields(
        document,
        {"schema", "schema_version", "storage", "certified_result"},
        "full result",
    )
    version = schema(value, FULL_RESULT_SCHEMA, "full result", (1, 2, 3))
    storage = fields(value["storage"], set(_STORAGE), "full result.storage")
    for name, expected in _STORAGE.items():
        if storage[name] != expected:
            raise ValueError(f"full result.storage.{name} must be {expected!r}")
    certified_document = value["certified_result"]
    if not isinstance(certified_document, dict):
        raise ValueError("full result.certified_result must be an object")
    if certified_document.get("schema_version") != version:
        raise ValueError("full result schema_version does not match certified result")
    return stability_result_from_dict(certified_document)


def _read_array(
    archive: zipfile.ZipFile, name: str, shape: Tuple[int, ...], source: Path
) -> np.ndarray:
    info = archive.getinfo(name)
    try:
        with archive.open(info) as stream:
            version = np.lib.format.read_magic(stream)
            if version != (1, 0):
                raise ValueError("NumPy format version must be 1.0")
            actual_shape, fortran_order, dtype = np.lib.format.read_array_header_1_0(
                stream
            )
            data_offset = stream.tell()
        expected_size = math.prod(shape) * 8
        if dtype.str != "<f8":
            raise ValueError(f"dtype is {dtype.str!r}; expected '<f8'")
        if fortran_order:
            raise ValueError("Fortran-order arrays are not supported")
        if actual_shape != shape:
            raise ValueError(f"shape is {actual_shape}; expected {shape}")
        if info.file_size != data_offset + expected_size:
            raise ValueError("payload size does not match its header")
        with archive.open(info) as stream:
            array = np.lib.format.read_array(stream, allow_pickle=False)
    except MemoryError as error:
        raise GlissAllocationError(
            f"unable to allocate {name} with shape {shape}"
        ) from error
    except (EOFError, OSError, ValueError, zipfile.BadZipFile) as error:
        raise ValueError(
            f"{source}: {name}: invalid or truncated NumPy array: {error}"
        ) from error
    array = np.asarray(array, dtype=np.float64, order="C")
    array.setflags(write=False)
    return array


def _read_full_result(
    archive: zipfile.ZipFile, metadata: Any, source: Path
) -> FullStabilityResult:
    certified = _parse_result_metadata(metadata)
    classes = []
    for parity_class, item in enumerate(certified.classes, start=1):
        count = _unknown_count(item, f"class {parity_class}")
        shapes: Dict[str, Tuple[int, ...]] = {
            name: (count,) for name in _ARRAY_ATTRIBUTES
        }
        shapes["eigenvectors"] = (count, count)
        arrays = tuple(
            _read_array(
                archive,
                _entry_name(parity_class, name),
                shapes[name],
                source,
            )
            for name in _ARRAY_ATTRIBUTES
        )
        classes.append(FullSpectrumResult(item, *arrays))
    result = FullStabilityResult((classes[0], classes[1]))
    _validate_full_result(result)
    return result


def write_full_result(result: FullStabilityResult, path: PathLike) -> None:
    """Atomically write a deterministic versioned full-spectrum container."""
    metadata = _full_result_metadata(result)
    _write_container(path, metadata, result)


def read_full_result(path: PathLike) -> FullStabilityResult:
    """Read and strictly validate a versioned full-spectrum container."""
    source = document_path(path, "input")
    with _open_archive(source) as archive:
        return _read_full_result(archive, _read_metadata(archive, source), source)


@dataclass(frozen=True, eq=False)
class FullRunManifest:
    """Inputs, complete spectra, checksums, and software provenance for one run."""

    equilibrium_filename: str
    equilibrium_size_bytes: int
    equilibrium_sha256: str
    equilibrium_schema_version: int
    configuration: StabilityConfiguration
    result: FullStabilityResult
    gliss_python_version: str
    gliss_native_version: str
    gliss_abi_version: int
    numpy_version: str
    python_version: str

    def __post_init__(self) -> None:
        _validate_full_result(self.result)
        self._base_manifest()

    def _base_manifest(self) -> RunManifest:
        return RunManifest(
            equilibrium_filename=self.equilibrium_filename,
            equilibrium_size_bytes=self.equilibrium_size_bytes,
            equilibrium_sha256=self.equilibrium_sha256,
            equilibrium_schema_version=self.equilibrium_schema_version,
            configuration=self.configuration,
            result=_certified_result(self.result),
            gliss_python_version=self.gliss_python_version,
            gliss_native_version=self.gliss_native_version,
            gliss_abi_version=self.gliss_abi_version,
            numpy_version=self.numpy_version,
            python_version=self.python_version,
        )

    @classmethod
    def _from_base(
        cls, base: RunManifest, result: FullStabilityResult
    ) -> "FullRunManifest":
        return cls(
            equilibrium_filename=base.equilibrium_filename,
            equilibrium_size_bytes=base.equilibrium_size_bytes,
            equilibrium_sha256=base.equilibrium_sha256,
            equilibrium_schema_version=base.equilibrium_schema_version,
            configuration=base.configuration,
            result=result,
            gliss_python_version=base.gliss_python_version,
            gliss_native_version=base.gliss_native_version,
            gliss_abi_version=base.gliss_abi_version,
            numpy_version=base.numpy_version,
            python_version=base.python_version,
        )

    def _metadata(self) -> Dict[str, Any]:
        document = self._base_manifest().to_dict()
        document["schema"] = FULL_RUN_SCHEMA
        document["result"] = _full_result_metadata(self.result)
        return document

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, FullRunManifest):
            return NotImplemented
        if self._metadata() != other._metadata():
            return False
        return all(
            np.array_equal(getattr(left, name), getattr(right, name))
            for left, right in zip(self.result.classes, other.result.classes)
            for name in _ARRAY_ATTRIBUTES
        )

    def write(self, path: PathLike) -> None:
        """Atomically write this deterministic full-run container."""
        _write_container(path, self._metadata(), self.result)

    @classmethod
    def read(cls, path: PathLike) -> "FullRunManifest":
        """Read and strictly validate a versioned full-run container."""
        source = document_path(path, "input")
        with _open_archive(source) as archive:
            metadata = _read_metadata(archive, source)
            value = fields(
                metadata,
                {
                    "schema",
                    "schema_version",
                    "equilibrium",
                    "software",
                    "configuration",
                    "result",
                },
                "full run",
            )
            schema(value, FULL_RUN_SCHEMA, "full run", (1, 2, 3))
            result = _read_full_result(archive, value["result"], source)
            base_document = dict(value)
            base_document["schema"] = "gliss.stability.run"
            base_document["result"] = value["result"]["certified_result"]
            base = RunManifest.from_dict(base_document)
            return cls._from_base(base, result)

    def verify_equilibrium(self, path: PathLike) -> None:
        """Verify a candidate input against the recorded size and SHA-256."""
        self._base_manifest().verify_equilibrium(path)


def _write_full_run_manifest(
    path: PathLike,
    equilibrium_path: PathLike,
    configuration: StabilityConfiguration,
    result: FullStabilityResult,
    expected_equilibrium: Optional[Tuple[int, int, str]] = None,
) -> FullRunManifest:
    certified = _validate_full_result(result)
    base = _create_run_manifest(
        equilibrium_path,
        configuration,
        certified,
        expected_equilibrium,
    )
    manifest = FullRunManifest._from_base(base, result)
    manifest.write(path)
    return manifest


def write_full_run_manifest(
    path: PathLike,
    equilibrium_path: PathLike,
    configuration: StabilityConfiguration,
    result: FullStabilityResult,
) -> FullRunManifest:
    """Create and write a self-contained manifest with both complete spectra."""
    return _write_full_run_manifest(path, equilibrium_path, configuration, result)
