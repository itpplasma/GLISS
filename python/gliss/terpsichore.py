"""Fixed- and free-boundary TERPSICHORE compatibility solves."""

import ctypes
import math
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Union

from . import _load_library, _require_symbols
from .equilibrium import GlissInternalError, _error_buffer, _raise_for_status
from ._stability_input import mode_integer

PathLike = Union[str, os.PathLike]


@dataclass(frozen=True)
class TerpsichoreFixedBoundaryResult:
    """Certified lowest negative eigenpair from an IVAC=0 FORT.23 file."""

    unknowns: int
    negative_count: int
    eigenvalue: float
    certificate: float
    residual: float
    resolution: float
    reference_eigenvalue: float
    reference_potential: float
    computed_potential: float
    reference_kinetic: float
    computed_kinetic: float
    reference_residual: float
    mode_overlap: float


@dataclass(frozen=True)
class TerpsichorePseudoplasmaResult:
    """Certified IVAC>0 eigenpair and TERPSICHORE reference diagnostics."""

    unknowns: int
    negative_count: int
    eigenvalue: float
    certificate: float
    residual: float
    resolution: float
    growth_rate: float
    reference_eigenvalue: float
    reference_potential: float
    computed_potential: float
    reference_kinetic: float
    computed_kinetic: float
    reference_residual: float
    mode_overlap: float


class _TerpsichoreFixedBoundaryResult(ctypes.Structure):
    _fields_ = [
        ("struct_size", ctypes.c_size_t),
        ("unknowns", ctypes.c_size_t),
        ("negative_count", ctypes.c_size_t),
        ("eigenvalue", ctypes.c_double),
        ("certificate", ctypes.c_double),
        ("residual", ctypes.c_double),
        ("resolution", ctypes.c_double),
        ("reference_eigenvalue", ctypes.c_double),
        ("reference_potential", ctypes.c_double),
        ("computed_potential", ctypes.c_double),
        ("reference_kinetic", ctypes.c_double),
        ("computed_kinetic", ctypes.c_double),
        ("reference_residual", ctypes.c_double),
        ("mode_overlap", ctypes.c_double),
    ]


class _TerpsichorePseudoplasmaResult(ctypes.Structure):
    _fields_ = [
        ("struct_size", ctypes.c_size_t),
        ("unknowns", ctypes.c_size_t),
        ("negative_count", ctypes.c_size_t),
        ("eigenvalue", ctypes.c_double),
        ("certificate", ctypes.c_double),
        ("residual", ctypes.c_double),
        ("resolution", ctypes.c_double),
        ("growth_rate", ctypes.c_double),
        ("reference_eigenvalue", ctypes.c_double),
        ("reference_potential", ctypes.c_double),
        ("computed_potential", ctypes.c_double),
        ("reference_kinetic", ctypes.c_double),
        ("computed_kinetic", ctypes.c_double),
        ("reference_residual", ctypes.c_double),
        ("mode_overlap", ctypes.c_double),
    ]


def _bind(library: Any) -> None:
    _require_symbols(
        library,
        ("gliss_terpsichore_fixed_boundary",),
        "TERPSICHORE fixed-boundary solver",
    )
    function = library.gliss_terpsichore_fixed_boundary
    function.argtypes = (
        ctypes.c_char_p,
        ctypes.c_size_t,
        ctypes.POINTER(_TerpsichoreFixedBoundaryResult),
        ctypes.c_void_p,
        ctypes.c_size_t,
    )
    function.restype = ctypes.c_int


def _bind_pseudoplasma(library: Any) -> None:
    _require_symbols(
        library,
        ("gliss_terpsichore_pseudoplasma",),
        "TERPSICHORE pseudoplasma solver",
    )
    function = library.gliss_terpsichore_pseudoplasma
    function.argtypes = (
        ctypes.c_char_p,
        ctypes.c_size_t,
        ctypes.c_int32,
        ctypes.c_char_p,
        ctypes.c_size_t,
        ctypes.POINTER(_TerpsichorePseudoplasmaResult),
        ctypes.c_void_p,
        ctypes.c_size_t,
    )
    function.restype = ctypes.c_int


def _fixture_path(path: PathLike, label: str = "FORT.23") -> bytes:
    try:
        value = os.fspath(path)
    except TypeError as error:
        raise TypeError("path must be a string or path-like object") from error
    if not isinstance(value, str):
        raise TypeError("path must resolve to a string")
    if "\0" in value:
        raise ValueError("path contains a null byte")
    fixture = Path(value)
    if not fixture.exists():
        raise FileNotFoundError(f"TERPSICHORE {label} does not exist: {fixture}")
    if not fixture.is_file():
        raise ValueError(f"TERPSICHORE {label} is not a file: {fixture}")
    encoded = os.fsencode(fixture)
    if len(encoded) > ctypes.c_size_t(-1).value:
        raise ValueError("encoded path is too long for the GLISS C API")
    return encoded


def _result(native: _TerpsichoreFixedBoundaryResult) -> TerpsichoreFixedBoundaryResult:
    values = tuple(getattr(native, name) for name, _ in native._fields_[3:])
    if (
        native.unknowns < 1
        or native.negative_count < 1
        or native.eigenvalue >= 0.0
        or native.certificate < 0.0
        or native.residual < 0.0
        or native.resolution < 0.0
        or native.reference_kinetic <= 0.0
        or native.computed_kinetic <= 0.0
        or native.reference_residual < 0.0
        or native.mode_overlap < 0.0
        or native.mode_overlap > 1.0 + 1.0e-12
        or not all(math.isfinite(value) for value in values)
    ):
        raise GlissInternalError("GLISS returned an invalid TERPSICHORE result")
    return TerpsichoreFixedBoundaryResult(
        unknowns=native.unknowns,
        negative_count=native.negative_count,
        eigenvalue=native.eigenvalue,
        certificate=native.certificate,
        residual=native.residual,
        resolution=native.resolution,
        reference_eigenvalue=native.reference_eigenvalue,
        reference_potential=native.reference_potential,
        computed_potential=native.computed_potential,
        reference_kinetic=native.reference_kinetic,
        computed_kinetic=native.computed_kinetic,
        reference_residual=native.reference_residual,
        mode_overlap=native.mode_overlap,
    )


def _pseudoplasma_result(
    native: _TerpsichorePseudoplasmaResult,
) -> TerpsichorePseudoplasmaResult:
    values = tuple(getattr(native, name) for name, _ in native._fields_[3:])
    invalid = (
        native.unknowns < 1
        or native.negative_count < 1
        or native.eigenvalue >= 0.0
        or native.certificate < 0.0
        or native.residual < 0.0
        or native.resolution < 0.0
        or native.growth_rate <= 0.0
        or native.reference_kinetic <= 0.0
        or native.computed_kinetic <= 0.0
        or native.reference_residual < 0.0
        or native.mode_overlap < 0.0
        or native.mode_overlap > 1.0 + 1.0e-12
        or not all(math.isfinite(value) for value in values)
    )
    if invalid:
        raise GlissInternalError(
            "GLISS returned an invalid TERPSICHORE pseudoplasma result"
        )
    return TerpsichorePseudoplasmaResult(
        unknowns=native.unknowns,
        negative_count=native.negative_count,
        eigenvalue=native.eigenvalue,
        certificate=native.certificate,
        residual=native.residual,
        resolution=native.resolution,
        growth_rate=native.growth_rate,
        reference_eigenvalue=native.reference_eigenvalue,
        reference_potential=native.reference_potential,
        computed_potential=native.computed_potential,
        reference_kinetic=native.reference_kinetic,
        computed_kinetic=native.computed_kinetic,
        reference_residual=native.reference_residual,
        mode_overlap=native.mode_overlap,
    )


def solve_terpsichore_fixed_boundary(
    path: PathLike,
) -> TerpsichoreFixedBoundaryResult:
    """Solve the lowest negative IVAC=0, MODELK=0 TERPSICHORE eigenpair."""

    encoded = _fixture_path(path)
    library = _load_library()
    _bind(library)
    native = _TerpsichoreFixedBoundaryResult(
        struct_size=ctypes.sizeof(_TerpsichoreFixedBoundaryResult)
    )
    error = _error_buffer()
    status = library.gliss_terpsichore_fixed_boundary(
        encoded, len(encoded), ctypes.byref(native), error, len(error)
    )
    _raise_for_status(status, error, "gliss_terpsichore_fixed_boundary")
    return _result(native)


def solve_terpsichore_pseudoplasma(
    matrix_path: PathLike,
    vacuum_intervals: int,
    vacuum_path: PathLike,
) -> TerpsichorePseudoplasmaResult:
    """Solve a MODELK=0 TERPSICHORE pressureless-pseudoplasma problem."""

    intervals = mode_integer(vacuum_intervals, "vacuum_intervals")
    if intervals <= 0:
        raise ValueError("vacuum_intervals must be positive")
    matrix = _fixture_path(matrix_path)
    vacuum = _fixture_path(vacuum_path, "FORT.24")
    library = _load_library()
    _bind_pseudoplasma(library)
    native = _TerpsichorePseudoplasmaResult(
        struct_size=ctypes.sizeof(_TerpsichorePseudoplasmaResult)
    )
    error = _error_buffer()
    status = library.gliss_terpsichore_pseudoplasma(
        matrix,
        len(matrix),
        intervals,
        vacuum,
        len(vacuum),
        ctypes.byref(native),
        error,
        len(error),
    )
    _raise_for_status(status, error, "gliss_terpsichore_pseudoplasma")
    return _pseudoplasma_result(native)
