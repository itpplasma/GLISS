"""Axisymmetric fixed-boundary stability calculations."""

import ctypes
import math
from dataclasses import dataclass
from typing import Any, Optional

from . import _require_symbols
from ._stability_input import mode_integer
from .equilibrium import (
    Equilibrium,
    GlissInternalError,
    _error_buffer,
    _raise_for_status,
)

@dataclass(frozen=True)
class AxisymmetricResult:
    """Inertia and optional certified eigenpair for one axisymmetric family."""

    has_eigenpair: bool
    field_periods: int
    toroidal_mode: int
    poloidal_max: int
    mode_count: int
    radial_surfaces: int
    parity_class: int
    degree: int
    negative_count: int
    lowest_eigenvalue: Optional[float]
    certificate: Optional[float]
    eigenpair_residual: Optional[float]
    force_balance_residual: float


class _AxisymmetricResult(ctypes.Structure):
    _fields_ = [
        ("struct_size", ctypes.c_size_t),
        ("has_eigenpair", ctypes.c_int32),
        ("field_periods", ctypes.c_int32),
        ("toroidal_mode", ctypes.c_int32),
        ("poloidal_max", ctypes.c_int32),
        ("mode_count", ctypes.c_size_t),
        ("radial_surfaces", ctypes.c_size_t),
        ("parity_class", ctypes.c_int32),
        ("degree", ctypes.c_int32),
        ("negative_count", ctypes.c_size_t),
        ("lowest_eigenvalue", ctypes.c_double),
        ("certificate", ctypes.c_double),
        ("eigenpair_residual", ctypes.c_double),
        ("force_balance_residual", ctypes.c_double),
    ]


def _bind(library: Any) -> None:
    _require_symbols(
        library,
        ("gliss_axisymmetric_spectrum",),
        "axisymmetric spectrum solver",
    )
    function = library.gliss_axisymmetric_spectrum
    function.argtypes = (
        ctypes.c_void_p,
        ctypes.c_int32,
        ctypes.c_int32,
        ctypes.c_int32,
        ctypes.c_int32,
        ctypes.POINTER(_AxisymmetricResult),
        ctypes.c_void_p,
        ctypes.c_size_t,
    )
    function.restype = ctypes.c_int


def _positive_mode(value: Any, name: str) -> int:
    result = mode_integer(value, name)
    if result <= 0:
        raise ValueError(f"{name} must be positive")
    return result


def _degree(value: Any) -> int:
    result = mode_integer(value, "degree")
    if result < 1 or result > 4:
        raise ValueError("degree must be between 1 and 4")
    return result


def _result(
    native: _AxisymmetricResult,
    toroidal_mode: int,
    poloidal_max: int,
    degree: int,
    solve_eigenpair: bool,
) -> AxisymmetricResult:
    metadata_valid = (
        native.has_eigenpair == int(solve_eigenpair)
        and native.field_periods == 1
        and native.toroidal_mode == toroidal_mode
        and native.poloidal_max == poloidal_max
        and native.mode_count == 2 * poloidal_max + 1
        and native.radial_surfaces >= 1
        and native.parity_class == 1
        and native.degree == degree
        and math.isfinite(native.force_balance_residual)
        and native.force_balance_residual >= 0.0
    )
    eigenpair = (
        native.lowest_eigenvalue,
        native.certificate,
        native.eigenpair_residual,
    )
    if solve_eigenpair:
        pair_valid = (
            all(math.isfinite(value) for value in eigenpair)
            and native.certificate >= 0.0
            and native.eigenpair_residual >= 0.0
        )
    else:
        pair_valid = all(math.isnan(value) for value in eigenpair)
    if not metadata_valid or not pair_valid:
        raise GlissInternalError("GLISS returned an invalid axisymmetric result")
    return AxisymmetricResult(
        has_eigenpair=solve_eigenpair,
        field_periods=native.field_periods,
        toroidal_mode=native.toroidal_mode,
        poloidal_max=native.poloidal_max,
        mode_count=native.mode_count,
        radial_surfaces=native.radial_surfaces,
        parity_class=native.parity_class,
        degree=native.degree,
        negative_count=native.negative_count,
        lowest_eigenvalue=native.lowest_eigenvalue if solve_eigenpair else None,
        certificate=native.certificate if solve_eigenpair else None,
        eigenpair_residual=native.eigenpair_residual if solve_eigenpair else None,
        force_balance_residual=native.force_balance_residual,
    )


def _calculate(
    equilibrium: Equilibrium,
    toroidal_mode: int,
    poloidal_max: int,
    degree: int,
    solve_eigenpair: bool,
) -> AxisymmetricResult:
    if not isinstance(equilibrium, Equilibrium):
        raise TypeError("equilibrium must be a gliss.Equilibrium")
    equilibrium._require_open()
    toroidal_mode = _positive_mode(toroidal_mode, "toroidal_mode")
    poloidal_max = _positive_mode(poloidal_max, "poloidal_max")
    degree = _degree(degree)
    _bind(equilibrium._library)
    native = _AxisymmetricResult(struct_size=ctypes.sizeof(_AxisymmetricResult))
    error = _error_buffer()
    status = equilibrium._library.gliss_axisymmetric_spectrum(
        equilibrium._handle,
        toroidal_mode,
        poloidal_max,
        degree,
        int(solve_eigenpair),
        ctypes.byref(native),
        error,
        len(error),
    )
    _raise_for_status(status, error, "gliss_axisymmetric_spectrum")
    return _result(
        native, toroidal_mode, poloidal_max, degree, solve_eigenpair
    )


def axisymmetric_inertia(
    equilibrium: Equilibrium,
    toroidal_mode: int = 1,
    poloidal_max: int = 8,
    degree: int = 2,
) -> AxisymmetricResult:
    """Count negative eigenvalues in one axisymmetric Fourier family."""

    return _calculate(
        equilibrium,
        toroidal_mode,
        poloidal_max,
        degree,
        solve_eigenpair=False,
    )


def solve_axisymmetric(
    equilibrium: Equilibrium,
    toroidal_mode: int = 1,
    poloidal_max: int = 8,
    degree: int = 2,
) -> AxisymmetricResult:
    """Return inertia and the certified lowest axisymmetric eigenpair."""

    return _calculate(
        equilibrium,
        toroidal_mode,
        poloidal_max,
        degree,
        solve_eigenpair=True,
    )
