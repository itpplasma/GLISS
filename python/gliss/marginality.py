"""CAS3D-compatible incompressible marginal-stability calculations."""

import ctypes
import math
from dataclasses import dataclass
from typing import Any, Optional, Sequence, Tuple

from . import _require_symbols
from ._stability_input import mode_integer, validate_modes
from .equilibrium import (
    Equilibrium,
    GlissInternalError,
    _error_buffer,
    _raise_for_status,
)

_QUADRATURE = {"midpoint": 1}


@dataclass(frozen=True)
class Cas3dMarginalityResult:
    """Inertia and optional eigenpair in the artificial CAS3D norm."""

    has_eigenpair: bool
    field_periods: int
    modes: Tuple[Tuple[int, int], ...]
    radial_surfaces: int
    parity_class: int
    radial_quadrature: str
    angular_resolution: Tuple[int, int]
    negative_count: int
    lowest_eigenvalue: Optional[float]
    certificate: Optional[float]
    eigenpair_residual: Optional[float]
    force_balance_residual: float
    normalization: str = (
        "CAS3D artificial L2 norm of transformed normal and tangential components"
    )
    interpretation: str = "stability and marginality only; not a physical growth rate"
    boundary_condition: str = "fixed"
    coordinate_handedness: str = "left-handed"
    fourier_convention: str = "2*pi*(m*theta - n*zeta/N_T)"


class _Cas3dMarginalityResult(ctypes.Structure):
    _fields_ = [
        ("struct_size", ctypes.c_size_t),
        ("has_eigenpair", ctypes.c_int32),
        ("field_periods", ctypes.c_int32),
        ("mode_count", ctypes.c_size_t),
        ("radial_surfaces", ctypes.c_size_t),
        ("parity_class", ctypes.c_int32),
        ("radial_quadrature", ctypes.c_int32),
        ("angular_theta", ctypes.c_int32),
        ("angular_zeta", ctypes.c_int32),
        ("negative_count", ctypes.c_size_t),
        ("lowest_eigenvalue", ctypes.c_double),
        ("certificate", ctypes.c_double),
        ("eigenpair_residual", ctypes.c_double),
        ("force_balance_residual", ctypes.c_double),
    ]


def _bind(library: Any) -> None:
    _require_symbols(
        library,
        ("gliss_cas3d_marginality",),
        "CAS3D marginality solver",
    )
    function = library.gliss_cas3d_marginality
    function.argtypes = (
        ctypes.c_void_p,
        ctypes.c_size_t,
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.c_int32,
        ctypes.c_int32,
        ctypes.c_int32,
        ctypes.c_int32,
        ctypes.c_int32,
        ctypes.POINTER(_Cas3dMarginalityResult),
        ctypes.c_void_p,
        ctypes.c_size_t,
    )
    function.restype = ctypes.c_int


def _angular_resolution(value: Any, name: str) -> int:
    result = mode_integer(value, name)
    if result < 4:
        raise ValueError(f"{name} must be at least 4")
    return result


def _parity_class(value: Any) -> int:
    result = mode_integer(value, "parity_class")
    if result not in (1, 2):
        raise ValueError("parity_class must be 1 or 2")
    return result


def _quadrature(value: Any) -> int:
    if not isinstance(value, str):
        raise TypeError("radial_quadrature must be a string")
    try:
        return _QUADRATURE[value]
    except KeyError as error:
        raise ValueError("radial_quadrature must be midpoint") from error


def _result(
    native: _Cas3dMarginalityResult,
    modes: Tuple[Tuple[int, int], ...],
    parity_class: int,
    radial_quadrature: int,
    angular_resolution: Tuple[int, int],
    solve_eigenpair: bool,
) -> Cas3dMarginalityResult:
    metadata_valid = (
        native.has_eigenpair == int(solve_eigenpair)
        and native.field_periods >= 1
        and native.mode_count == len(modes)
        and native.radial_surfaces >= 2
        and native.parity_class == parity_class
        and native.radial_quadrature == radial_quadrature
        and (native.angular_theta, native.angular_zeta) == angular_resolution
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
        raise GlissInternalError("GLISS returned an invalid CAS3D marginality result")
    quadrature = {value: name for name, value in _QUADRATURE.items()}
    return Cas3dMarginalityResult(
        has_eigenpair=solve_eigenpair,
        field_periods=native.field_periods,
        modes=modes,
        radial_surfaces=native.radial_surfaces,
        parity_class=native.parity_class,
        radial_quadrature=quadrature[native.radial_quadrature],
        angular_resolution=angular_resolution,
        negative_count=native.negative_count,
        lowest_eigenvalue=native.lowest_eigenvalue if solve_eigenpair else None,
        certificate=native.certificate if solve_eigenpair else None,
        eigenpair_residual=(native.eigenpair_residual if solve_eigenpair else None),
        force_balance_residual=native.force_balance_residual,
    )


def _calculate(
    equilibrium: Equilibrium,
    modes: Sequence[Tuple[int, int]],
    parity_class: int,
    radial_quadrature: str,
    angular_theta: int,
    angular_zeta: int,
    solve_eigenpair: bool,
) -> Cas3dMarginalityResult:
    if not isinstance(equilibrium, Equilibrium):
        raise TypeError("equilibrium must be a gliss.Equilibrium")
    equilibrium._require_open()
    validated_modes = validate_modes(modes)
    validated_parity = _parity_class(parity_class)
    quadrature = _quadrature(radial_quadrature)
    resolution = (
        _angular_resolution(angular_theta, "angular_theta"),
        _angular_resolution(angular_zeta, "angular_zeta"),
    )
    _bind(equilibrium._library)
    count = len(validated_modes)
    integers = ctypes.c_int32 * count
    mode_m = integers(*(mode[0] for mode in validated_modes))
    mode_n = integers(*(mode[1] for mode in validated_modes))
    native = _Cas3dMarginalityResult(struct_size=ctypes.sizeof(_Cas3dMarginalityResult))
    error = _error_buffer()
    status = equilibrium._library.gliss_cas3d_marginality(
        equilibrium._handle,
        count,
        mode_m,
        mode_n,
        validated_parity,
        quadrature,
        resolution[0],
        resolution[1],
        int(solve_eigenpair),
        ctypes.byref(native),
        error,
        len(error),
    )
    _raise_for_status(status, error, "gliss_cas3d_marginality")
    return _result(
        native,
        validated_modes,
        validated_parity,
        quadrature,
        resolution,
        solve_eigenpair,
    )


def cas3d_marginality_inertia(
    equilibrium: Equilibrium,
    modes: Sequence[Tuple[int, int]],
    parity_class: int = 1,
    radial_quadrature: str = "midpoint",
    angular_theta: int = 64,
    angular_zeta: int = 64,
) -> Cas3dMarginalityResult:
    """Count negative directions in the artificial CAS3D normalization."""

    return _calculate(
        equilibrium,
        modes,
        parity_class,
        radial_quadrature,
        angular_theta,
        angular_zeta,
        solve_eigenpair=False,
    )


def solve_cas3d_marginality(
    equilibrium: Equilibrium,
    modes: Sequence[Tuple[int, int]],
    parity_class: int = 1,
    radial_quadrature: str = "midpoint",
    angular_theta: int = 64,
    angular_zeta: int = 64,
) -> Cas3dMarginalityResult:
    """Return inertia and the lowest artificial-norm CAS3D eigenpair."""

    return _calculate(
        equilibrium,
        modes,
        parity_class,
        radial_quadrature,
        angular_theta,
        angular_zeta,
        solve_eigenpair=True,
    )
