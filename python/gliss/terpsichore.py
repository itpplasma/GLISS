"""Fixed-boundary TERPSICHORE compatibility solves."""

import ctypes
import math
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Union

from . import _load_library, _require_symbols
from .equilibrium import GlissInternalError, _error_buffer, _raise_for_status

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


class _TerpsichoreFixedBoundaryResult(ctypes.Structure):
    _fields_ = [
        ("struct_size", ctypes.c_size_t),
        ("unknowns", ctypes.c_size_t),
        ("negative_count", ctypes.c_size_t),
        ("eigenvalue", ctypes.c_double),
        ("certificate", ctypes.c_double),
        ("residual", ctypes.c_double),
        ("resolution", ctypes.c_double),
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


def _fixture_path(path: PathLike) -> bytes:
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
        raise FileNotFoundError(f"TERPSICHORE FORT.23 does not exist: {fixture}")
    if not fixture.is_file():
        raise ValueError(f"TERPSICHORE FORT.23 is not a file: {fixture}")
    encoded = os.fsencode(fixture)
    if len(encoded) > ctypes.c_size_t(-1).value:
        raise ValueError("encoded path is too long for the GLISS C API")
    return encoded


def _result(native: _TerpsichoreFixedBoundaryResult) -> TerpsichoreFixedBoundaryResult:
    values = (native.eigenvalue, native.certificate, native.residual, native.resolution)
    if (
        native.unknowns < 1
        or native.negative_count < 1
        or native.eigenvalue >= 0.0
        or native.certificate < 0.0
        or native.residual < 0.0
        or native.resolution < 0.0
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
