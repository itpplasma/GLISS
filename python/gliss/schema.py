"""Versioned configuration, result, and run-manifest schemas."""

import hashlib
import platform
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Mapping, Tuple

import numpy as np

from ._result_schema import (
    stability_result_from_dict,
    stability_result_to_dict,
)
from ._schema_support import (
    SCHEMA_VERSION,
    fields,
    integer,
    read_json,
    schema,
    string,
    write_json,
)
from .equilibrium import PathLike, _export_path
from .stability import (
    _QUADRATURE,
    _real_parameter,
    _validate_modes,
    StabilityProblem,
    StabilityResult,
)

_CONFIGURATION_SCHEMA = "gliss.stability.configuration"
_RUN_SCHEMA = "gliss.stability.run"
_EQUILIBRIUM_FORMAT = "gvec-cas3d-netcdf"


@dataclass(frozen=True)
class StabilityConfiguration:
    """Validated fixed-boundary problem inputs with JSON interchange."""

    modes: Tuple[Tuple[int, int], ...]
    adiabatic_index: float = 5.0 / 3.0
    density_kg_m3: float = 1.0
    zero_floor: float = 1.0
    radial_quadrature: str = "midpoint"

    def __post_init__(self) -> None:
        if not isinstance(self.radial_quadrature, str):
            raise TypeError("radial_quadrature must be a string")
        if self.radial_quadrature not in _QUADRATURE:
            raise ValueError("radial_quadrature must be 'midpoint' or 'gauss2'")
        object.__setattr__(self, "modes", _validate_modes(self.modes))
        object.__setattr__(
            self,
            "adiabatic_index",
            _real_parameter(self.adiabatic_index, "adiabatic_index", True),
        )
        object.__setattr__(
            self,
            "density_kg_m3",
            _real_parameter(self.density_kg_m3, "density_kg_m3", False),
        )
        floor = _real_parameter(self.zero_floor, "zero_floor", False)
        if floor > 0.125 * np.finfo(np.float64).max:
            raise ValueError("zero_floor is too large for spectrum certification")
        object.__setattr__(self, "zero_floor", floor)

    def create_problem(self, equilibrium: Any) -> StabilityProblem:
        """Create an assembled problem from this immutable configuration."""
        return StabilityProblem(
            equilibrium,
            self.modes,
            self.adiabatic_index,
            self.density_kg_m3,
            self.zero_floor,
            self.radial_quadrature,
        )

    def to_dict(self) -> Dict[str, Any]:
        """Return the canonical version-1 configuration document."""
        return {
            "schema": _CONFIGURATION_SCHEMA,
            "schema_version": SCHEMA_VERSION,
            "boundary_condition": "fixed",
            "modes": [list(mode) for mode in self.modes],
            "adiabatic_index": self.adiabatic_index,
            "density_kg_m3": self.density_kg_m3,
            "zero_floor": self.zero_floor,
            "radial_quadrature": self.radial_quadrature,
        }

    @classmethod
    def from_dict(cls, document: Mapping[str, Any]) -> "StabilityConfiguration":
        """Validate and construct a version-1 configuration document."""
        expected = {
            "schema",
            "schema_version",
            "boundary_condition",
            "modes",
            "adiabatic_index",
            "density_kg_m3",
            "zero_floor",
            "radial_quadrature",
        }
        value = fields(document, expected, "configuration")
        schema(value, _CONFIGURATION_SCHEMA, "configuration")
        if value["boundary_condition"] != "fixed":
            raise ValueError("configuration.boundary_condition must be 'fixed'")
        try:
            return cls(
                modes=value["modes"],
                adiabatic_index=value["adiabatic_index"],
                density_kg_m3=value["density_kg_m3"],
                zero_floor=value["zero_floor"],
                radial_quadrature=value["radial_quadrature"],
            )
        except (TypeError, ValueError) as error:
            raise ValueError(f"configuration: {error}") from error

    def write(self, path: PathLike) -> None:
        """Atomically write the canonical JSON configuration."""
        write_json(path, self.to_dict())

    @classmethod
    def read(cls, path: PathLike) -> "StabilityConfiguration":
        """Read and strictly validate a JSON configuration."""
        return cls.from_dict(read_json(path))


@dataclass(frozen=True, eq=False)
class RunManifest:
    """Portable inputs, outputs, checksums, and software provenance for one run."""

    equilibrium_filename: str
    equilibrium_size_bytes: int
    equilibrium_sha256: str
    configuration: StabilityConfiguration
    result: StabilityResult
    gliss_python_version: str
    gliss_native_version: str
    gliss_abi_version: int
    numpy_version: str
    python_version: str

    def __post_init__(self) -> None:
        filename = string(self.equilibrium_filename, "equilibrium_filename")
        if "\0" in filename:
            raise ValueError("equilibrium_filename contains a null byte")
        if Path(filename).name != filename or "/" in filename or "\\" in filename:
            raise ValueError("equilibrium_filename must be a base name without a path")
        integer(self.equilibrium_size_bytes, "equilibrium_size_bytes")
        digest = string(self.equilibrium_sha256, "equilibrium_sha256")
        if re.fullmatch(r"[0-9a-f]{64}", digest) is None:
            raise ValueError("equilibrium_sha256 must be 64 lowercase hex digits")
        if not isinstance(self.configuration, StabilityConfiguration):
            raise TypeError("configuration must be a gliss.StabilityConfiguration")
        if not isinstance(self.result, StabilityResult):
            raise TypeError("result must be a gliss.StabilityResult")
        string(self.gliss_python_version, "gliss_python_version")
        string(self.gliss_native_version, "gliss_native_version")
        abi = integer(self.gliss_abi_version, "gliss_abi_version", 1)
        if abi != 1:
            raise ValueError("gliss_abi_version is incompatible; expected 1")
        string(self.numpy_version, "numpy_version")
        string(self.python_version, "python_version")
        stability_result_to_dict(self.result)
        _validate_result_configuration(self.configuration, self.result)

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, RunManifest):
            return NotImplemented
        return self.to_dict() == other.to_dict()

    def to_dict(self) -> Dict[str, Any]:
        """Return the canonical version-1 run document."""
        return {
            "schema": _RUN_SCHEMA,
            "schema_version": SCHEMA_VERSION,
            "equilibrium": {
                "format": _EQUILIBRIUM_FORMAT,
                "schema_version": 1,
                "filename": self.equilibrium_filename,
                "size_bytes": self.equilibrium_size_bytes,
                "sha256": self.equilibrium_sha256,
            },
            "software": {
                "gliss_python": self.gliss_python_version,
                "gliss_native": self.gliss_native_version,
                "gliss_abi": self.gliss_abi_version,
                "numpy": self.numpy_version,
                "python": self.python_version,
            },
            "configuration": self.configuration.to_dict(),
            "result": stability_result_to_dict(self.result),
        }

    def write(self, path: PathLike) -> None:
        """Atomically write this canonical run manifest."""
        write_json(path, self.to_dict())

    @classmethod
    def read(cls, path: PathLike) -> "RunManifest":
        """Read and strictly validate a version-1 run manifest."""
        document = read_json(path)
        expected = {
            "schema",
            "schema_version",
            "equilibrium",
            "software",
            "configuration",
            "result",
        }
        value = fields(document, expected, "run")
        schema(value, _RUN_SCHEMA, "run")
        equilibrium = fields(
            value["equilibrium"],
            {"format", "schema_version", "filename", "size_bytes", "sha256"},
            "run.equilibrium",
        )
        if equilibrium["format"] != _EQUILIBRIUM_FORMAT:
            raise ValueError(f"run.equilibrium.format must be {_EQUILIBRIUM_FORMAT!r}")
        if equilibrium["schema_version"] != 1:
            raise ValueError("run.equilibrium.schema_version must be 1")
        digest = string(equilibrium["sha256"], "run.equilibrium.sha256")
        if re.fullmatch(r"[0-9a-f]{64}", digest) is None:
            raise ValueError("run.equilibrium.sha256 must be 64 lowercase hex digits")
        software = fields(
            value["software"],
            {"gliss_python", "gliss_native", "gliss_abi", "numpy", "python"},
            "run.software",
        )
        abi = integer(software["gliss_abi"], "run.software.gliss_abi", 1)
        if abi != 1:
            raise ValueError("run.software.gliss_abi is incompatible; expected 1")
        configuration = StabilityConfiguration.from_dict(value["configuration"])
        result = stability_result_from_dict(value["result"])
        _validate_result_configuration(configuration, result)
        return cls(
            equilibrium_filename=string(
                equilibrium["filename"], "run.equilibrium.filename"
            ),
            equilibrium_size_bytes=integer(
                equilibrium["size_bytes"], "run.equilibrium.size_bytes"
            ),
            equilibrium_sha256=digest,
            configuration=configuration,
            result=result,
            gliss_python_version=string(
                software["gliss_python"], "run.software.gliss_python"
            ),
            gliss_native_version=string(
                software["gliss_native"], "run.software.gliss_native"
            ),
            gliss_abi_version=abi,
            numpy_version=string(software["numpy"], "run.software.numpy"),
            python_version=string(software["python"], "run.software.python"),
        )

    def verify_equilibrium(self, path: PathLike) -> None:
        """Verify a candidate input against the recorded size and SHA-256."""
        export, _ = _export_path(path)
        if export.stat().st_size != self.equilibrium_size_bytes:
            raise ValueError("equilibrium size does not match the run manifest")
        if _sha256(export) != self.equilibrium_sha256:
            raise ValueError("equilibrium SHA-256 does not match the run manifest")


def _validate_result_configuration(
    configuration: StabilityConfiguration, result: StabilityResult
) -> None:
    reference = result.classes[0]
    if reference.modes != configuration.modes:
        raise ValueError("result modes do not match configuration modes")
    names = ("adiabatic_index", "density_kg_m3", "zero_floor", "radial_quadrature")
    for name in names:
        if getattr(reference, name) != getattr(configuration, name):
            raise ValueError(f"result {name} does not match configuration {name}")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _native_version() -> str:
    from . import version

    return version()


def write_run_manifest(
    path: PathLike,
    equilibrium_path: PathLike,
    configuration: StabilityConfiguration,
    result: StabilityResult,
) -> RunManifest:
    """Create and atomically write a portable manifest for one stability run."""
    if not isinstance(configuration, StabilityConfiguration):
        raise TypeError("configuration must be a gliss.StabilityConfiguration")
    if not isinstance(result, StabilityResult):
        raise TypeError("result must be a gliss.StabilityResult")
    stability_result_to_dict(result)
    _validate_result_configuration(configuration, result)
    export, _ = _export_path(equilibrium_path)
    from . import __version__

    manifest = RunManifest(
        equilibrium_filename=export.name,
        equilibrium_size_bytes=export.stat().st_size,
        equilibrium_sha256=_sha256(export),
        configuration=configuration,
        result=result,
        gliss_python_version=__version__,
        gliss_native_version=_native_version(),
        gliss_abi_version=1,
        numpy_version=np.__version__,
        python_version=platform.python_version(),
    )
    manifest.write(path)
    return manifest
