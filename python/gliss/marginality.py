"""CAS3D-compatible incompressible marginal-stability calculations."""

import ctypes
import math
import numbers
import sys
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

_PERPENDICULAR_L2 = "perpendicular_l2"
_CAS3D2MN_COEFFICIENT = "cas3d2mn_coefficient"
_RADIAL_QUADRATURES = {"gauss5": 1, "cas3d_midpoint": 2}
_NORMALIZATION_LABELS = {
    _PERPENDICULAR_L2: "compatible perpendicular L2 norm on physical modes",
    _CAS3D2MN_COEFFICIENT: (
        "Schwab CAS3D2MN labeled-envelope coefficient norm: half identity, "
        "physical-sideband stiffness pullback, f_l=1"
    ),
}


@dataclass(frozen=True)
class Cas3dMarginalityResult:
    """Inertia and optional eigenpair from the compatible FEEC problem."""

    has_eigenpair: bool
    field_periods: int
    modes: Tuple[Tuple[int, int], ...]
    radial_surfaces: int
    parity_class: int
    degree: int
    angular_resolution: Tuple[int, int]
    negative_count: int
    lowest_eigenvalue: Optional[float]
    certificate: Optional[float]
    eigenpair_residual: Optional[float]
    force_balance_residual: float
    normalization: str = (
        "compatible perpendicular L2 norm of normal and tangential components"
    )
    interpretation: str = "stability and marginality only; not a physical growth rate"
    boundary_condition: str = "fixed"
    coordinate_handedness: str = "left-handed"
    fourier_convention: str = "2*pi*(m*theta - n*zeta/N_T)"


@dataclass(frozen=True)
class Cas3dPhaseEnvelopeResult:
    """Inertia and optional eigenpair for a CAS3D2MN envelope."""

    has_eigenpair: bool
    field_periods: int
    base_mode: Tuple[int, int]
    envelope_modes: Tuple[Tuple[int, int], ...]
    labeled_sideband_count: int
    radial_surfaces: int
    parity_class: int
    degree: int
    angular_resolution: Tuple[int, int]
    negative_count: int
    lowest_eigenvalue: Optional[float]
    certificate: Optional[float]
    eigenpair_residual: Optional[float]
    force_balance_residual: float
    coefficient_angular_resolution: Optional[Tuple[int, int]] = None
    reference_length: Optional[float] = None
    radial_quadrature: str = "gauss5"
    inertia_zero_floor: float = 1.0e-12
    normalization: str = "compatible perpendicular L2 norm on physical modes"
    interpretation: str = "stability and marginality only; not a physical growth rate"
    boundary_condition: str = "fixed"
    coordinate_handedness: str = "left-handed"
    base_fourier_convention: str = "2*pi*(M*theta - N*zeta/N_T)"
    envelope_fourier_convention: str = "2*pi*(m*theta - n*zeta)"


class _Cas3dMarginalityResult(ctypes.Structure):
    _fields_ = [
        ("struct_size", ctypes.c_size_t),
        ("has_eigenpair", ctypes.c_int32),
        ("field_periods", ctypes.c_int32),
        ("mode_count", ctypes.c_size_t),
        ("radial_surfaces", ctypes.c_size_t),
        ("parity_class", ctypes.c_int32),
        ("degree", ctypes.c_int32),
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


def _bind_phase_envelope(library: Any, symbol: str) -> None:
    _require_symbols(
        library,
        (symbol,),
        "CAS3D2MN phase-envelope solver",
    )
    function = getattr(library, symbol)
    argument_types = (
        ctypes.c_void_p,
        ctypes.c_int32,
        ctypes.c_int32,
        ctypes.c_size_t,
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.c_int32,
        ctypes.c_int32,
        ctypes.c_int32,
        ctypes.c_int32,
    )
    if symbol == "gliss_cas3d2mn_phase_envelope":
        argument_types += (
            ctypes.c_int32,
            ctypes.c_int32,
            ctypes.c_double,
            ctypes.c_int32,
        )
    function.argtypes = argument_types + (
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


def _degree(value: Any) -> int:
    result = mode_integer(value, "degree")
    if result < 1 or result > 4:
        raise ValueError("degree must be between 1 and 4")
    return result


def _normalization(value: Any) -> str:
    if not isinstance(value, str):
        raise TypeError("normalization must be a string")
    if value not in _NORMALIZATION_LABELS:
        choices = ", ".join(repr(choice) for choice in _NORMALIZATION_LABELS)
        raise ValueError(f"normalization must be one of {choices}")
    return value


def _coefficient_resolution(value: Any) -> Tuple[int, int]:
    if not isinstance(value, (tuple, list)) or len(value) != 2:
        raise TypeError("coefficient_angular_resolution must be a pair")
    result = (
        mode_integer(value[0], "coefficient angular theta"),
        mode_integer(value[1], "coefficient angular zeta"),
    )
    if result[0] < 1 or result[1] < 1:
        raise ValueError("coefficient angular resolutions must be positive")
    return result


def _reference_length(value: Any) -> float:
    if isinstance(value, bool) or not isinstance(value, numbers.Real):
        raise TypeError("reference_length must be a real number")
    result = float(value)
    if not math.isfinite(result) or result <= 0.0:
        raise ValueError("reference_length must be finite and positive")
    cube_limit = sys.float_info.max ** (1.0 / 3.0)
    cube_floor = sys.float_info.min ** (1.0 / 3.0)
    if result < cube_floor or result > cube_limit:
        raise ValueError("reference_length cube must remain in the normal binary64 range")
    return result


def _radial_quadrature(value: Any) -> Tuple[str, int]:
    if not isinstance(value, str):
        raise TypeError("radial_quadrature must be a string")
    if value not in _RADIAL_QUADRATURES:
        choices = ", ".join(repr(choice) for choice in _RADIAL_QUADRATURES)
        raise ValueError(f"radial_quadrature must be one of {choices}")
    return value, _RADIAL_QUADRATURES[value]


def _result(
    native: _Cas3dMarginalityResult,
    modes: Tuple[Tuple[int, int], ...],
    parity_class: int,
    degree: int,
    angular_resolution: Tuple[int, int],
    solve_eigenpair: bool,
) -> Cas3dMarginalityResult:
    metadata_valid = (
        native.has_eigenpair == int(solve_eigenpair)
        and native.field_periods >= 1
        and native.mode_count == len(modes)
        and native.radial_surfaces >= 2
        and native.parity_class == parity_class
        and native.degree == degree
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
    return Cas3dMarginalityResult(
        has_eigenpair=solve_eigenpair,
        field_periods=native.field_periods,
        modes=modes,
        radial_surfaces=native.radial_surfaces,
        parity_class=native.parity_class,
        degree=native.degree,
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
    degree: int,
    angular_theta: int,
    angular_zeta: int,
    solve_eigenpair: bool,
) -> Cas3dMarginalityResult:
    if not isinstance(equilibrium, Equilibrium):
        raise TypeError("equilibrium must be a gliss.Equilibrium")
    equilibrium._require_open()
    validated_modes = validate_modes(modes)
    validated_parity = _parity_class(parity_class)
    degree = _degree(degree)
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
        degree,
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
        degree,
        resolution,
        solve_eigenpair,
    )


def _mode_pair(value: Any, name: str) -> Tuple[int, int]:
    if not isinstance(value, (tuple, list)) or len(value) != 2:
        raise TypeError(f"{name} must be an (m, n) pair")
    return (
        mode_integer(value[0], f"{name} poloidal mode"),
        mode_integer(value[1], f"{name} toroidal mode"),
    )


def _validate_phase_envelope(
    base_mode: Tuple[int, int],
    envelope_modes: Sequence[Tuple[int, int]],
) -> Tuple[Tuple[int, int], Tuple[Tuple[int, int], ...]]:
    base = _mode_pair(base_mode, "base_mode")
    if base[0] < 0:
        raise ValueError("base_mode poloidal mode must be nonnegative")
    if base[0] == 0 and base[1] < 0:
        raise ValueError("base_mode axis mode requires nonnegative n")
    try:
        entries = tuple(envelope_modes)
    except TypeError as error:
        raise TypeError("envelope_modes must be a sequence of (m, n) pairs") from error
    if not entries:
        raise ValueError("envelope_modes must not be empty")
    result = []
    seen = set()
    for index, entry in enumerate(entries):
        pair = _mode_pair(entry, f"envelope_modes[{index}]")
        if pair[0] < 0:
            raise ValueError(
                f"envelope_modes[{index}] poloidal mode must be nonnegative"
            )
        if pair in seen:
            raise ValueError(f"duplicate envelope mode {pair!r}")
        seen.add(pair)
        result.append(pair)
    if result[0] != (0, 0):
        raise ValueError("envelope_modes must begin with (0, 0)")
    return base, tuple(result)


def _phase_envelope_result(
    native: _Cas3dMarginalityResult,
    base_mode: Tuple[int, int],
    envelope_modes: Tuple[Tuple[int, int], ...],
    parity_class: int,
    degree: int,
    angular_resolution: Tuple[int, int],
    solve_eigenpair: bool,
    normalization: str,
    coefficient_resolution: Optional[Tuple[int, int]],
    reference_length: Optional[float],
    radial_quadrature: str,
) -> Cas3dPhaseEnvelopeResult:
    sideband_count = 2 * len(envelope_modes) - 1
    metadata_valid = (
        native.has_eigenpair == int(solve_eigenpair)
        and native.field_periods >= 1
        and native.mode_count == sideband_count
        and native.radial_surfaces >= 2
        and native.parity_class == parity_class
        and native.degree == degree
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
        raise GlissInternalError("GLISS returned an invalid CAS3D2MN result")
    return Cas3dPhaseEnvelopeResult(
        has_eigenpair=solve_eigenpair,
        field_periods=native.field_periods,
        base_mode=base_mode,
        envelope_modes=envelope_modes,
        labeled_sideband_count=sideband_count,
        radial_surfaces=native.radial_surfaces,
        parity_class=native.parity_class,
        degree=native.degree,
        angular_resolution=angular_resolution,
        negative_count=native.negative_count,
        lowest_eigenvalue=native.lowest_eigenvalue if solve_eigenpair else None,
        certificate=native.certificate if solve_eigenpair else None,
        eigenpair_residual=(native.eigenpair_residual if solve_eigenpair else None),
        force_balance_residual=native.force_balance_residual,
        normalization=_NORMALIZATION_LABELS[normalization],
        coefficient_angular_resolution=coefficient_resolution,
        reference_length=reference_length,
        radial_quadrature=radial_quadrature,
    )


def _calculate_phase_envelope(
    equilibrium: Equilibrium,
    base_mode: Tuple[int, int],
    envelope_modes: Sequence[Tuple[int, int]],
    parity_class: int,
    degree: int,
    angular_theta: int,
    angular_zeta: int,
    solve_eigenpair: bool,
    normalization: str,
    coefficient_angular_resolution: Optional[Tuple[int, int]],
    reference_length: Optional[float],
    radial_quadrature: str,
) -> Cas3dPhaseEnvelopeResult:
    if not isinstance(equilibrium, Equilibrium):
        raise TypeError("equilibrium must be a gliss.Equilibrium")
    equilibrium._require_open()
    base, envelopes = _validate_phase_envelope(base_mode, envelope_modes)
    validated_parity = _parity_class(parity_class)
    degree = _degree(degree)
    normalization = _normalization(normalization)
    quadrature_name, quadrature_policy = _radial_quadrature(radial_quadrature)
    if normalization == _CAS3D2MN_COEFFICIENT and degree != 1:
        raise ValueError("cas3d2mn_coefficient normalization requires degree=1")
    if normalization == _CAS3D2MN_COEFFICIENT:
        coefficient_resolution = _coefficient_resolution(coefficient_angular_resolution)
        length_scale = _reference_length(reference_length)
    else:
        if coefficient_angular_resolution is not None:
            raise ValueError(
                "coefficient_angular_resolution requires cas3d2mn_coefficient "
                "normalization"
            )
        if reference_length is not None:
            raise ValueError(
                "reference_length requires cas3d2mn_coefficient normalization"
            )
        coefficient_resolution = None
        length_scale = None
        if quadrature_name != "gauss5":
            raise ValueError(
                "cas3d_midpoint radial quadrature requires "
                "cas3d2mn_coefficient normalization"
            )
    resolution = (
        _angular_resolution(angular_theta, "angular_theta"),
        _angular_resolution(angular_zeta, "angular_zeta"),
    )
    if normalization == _PERPENDICULAR_L2:
        symbol = "gliss_cas3d_phase_envelope"
    else:
        symbol = "gliss_cas3d2mn_phase_envelope"
    _bind_phase_envelope(equilibrium._library, symbol)
    count = len(envelopes)
    integers = ctypes.c_int32 * count
    envelope_m = integers(*(mode[0] for mode in envelopes))
    envelope_n = integers(*(mode[1] for mode in envelopes))
    native = _Cas3dMarginalityResult(struct_size=ctypes.sizeof(_Cas3dMarginalityResult))
    error = _error_buffer()
    arguments = [
        equilibrium._handle,
        base[0],
        base[1],
        count,
        envelope_m,
        envelope_n,
        validated_parity,
        degree,
        resolution[0],
        resolution[1],
    ]
    if coefficient_resolution is not None:
        arguments.extend(
            [
                coefficient_resolution[0],
                coefficient_resolution[1],
                length_scale,
                quadrature_policy,
            ]
        )
    arguments.extend([int(solve_eigenpair), ctypes.byref(native), error, len(error)])
    status = getattr(equilibrium._library, symbol)(*arguments)
    _raise_for_status(status, error, symbol)
    return _phase_envelope_result(
        native,
        base,
        envelopes,
        validated_parity,
        degree,
        resolution,
        solve_eigenpair,
        normalization,
        coefficient_resolution,
        length_scale,
        quadrature_name,
    )


def cas3d_phase_envelope_inertia(
    equilibrium: Equilibrium,
    base_mode: Tuple[int, int],
    envelope_modes: Sequence[Tuple[int, int]],
    parity_class: int = 1,
    degree: int = 2,
    angular_theta: int = 64,
    angular_zeta: int = 64,
    normalization: str = _PERPENDICULAR_L2,
    coefficient_angular_resolution: Optional[Tuple[int, int]] = None,
    reference_length: Optional[float] = None,
    radial_quadrature: str = "gauss5",
) -> Cas3dPhaseEnvelopeResult:
    """Count negative directions in a CAS3D2MN phase envelope."""

    return _calculate_phase_envelope(
        equilibrium,
        base_mode,
        envelope_modes,
        parity_class,
        degree,
        angular_theta,
        angular_zeta,
        solve_eigenpair=False,
        normalization=normalization,
        coefficient_angular_resolution=coefficient_angular_resolution,
        reference_length=reference_length,
        radial_quadrature=radial_quadrature,
    )


def solve_cas3d_phase_envelope(
    equilibrium: Equilibrium,
    base_mode: Tuple[int, int],
    envelope_modes: Sequence[Tuple[int, int]],
    parity_class: int = 1,
    degree: int = 2,
    angular_theta: int = 64,
    angular_zeta: int = 64,
    normalization: str = _PERPENDICULAR_L2,
    coefficient_angular_resolution: Optional[Tuple[int, int]] = None,
    reference_length: Optional[float] = None,
    radial_quadrature: str = "gauss5",
) -> Cas3dPhaseEnvelopeResult:
    """Return inertia and the lowest pair for a CAS3D2MN phase envelope."""

    return _calculate_phase_envelope(
        equilibrium,
        base_mode,
        envelope_modes,
        parity_class,
        degree,
        angular_theta,
        angular_zeta,
        solve_eigenpair=True,
        normalization=normalization,
        coefficient_angular_resolution=coefficient_angular_resolution,
        reference_length=reference_length,
        radial_quadrature=radial_quadrature,
    )


def cas3d_marginality_inertia(
    equilibrium: Equilibrium,
    modes: Sequence[Tuple[int, int]],
    parity_class: int = 1,
    degree: int = 2,
    angular_theta: int = 64,
    angular_zeta: int = 64,
) -> Cas3dMarginalityResult:
    """Count negative directions in the compatible FEEC discretization."""

    return _calculate(
        equilibrium,
        modes,
        parity_class,
        degree,
        angular_theta,
        angular_zeta,
        solve_eigenpair=False,
    )


def solve_cas3d_marginality(
    equilibrium: Equilibrium,
    modes: Sequence[Tuple[int, int]],
    parity_class: int = 1,
    degree: int = 2,
    angular_theta: int = 64,
    angular_zeta: int = 64,
) -> Cas3dMarginalityResult:
    """Return inertia and the lowest compatible FEEC eigenpair."""

    return _calculate(
        equilibrium,
        modes,
        parity_class,
        degree,
        angular_theta,
        angular_zeta,
        solve_eigenpair=True,
    )
