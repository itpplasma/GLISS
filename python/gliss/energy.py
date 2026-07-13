"""Physical energy decomposition for fixed-boundary displacement vectors."""

import ctypes
import math
from dataclasses import dataclass
from typing import Any, Tuple, TYPE_CHECKING

import numpy as np

from . import _require_symbols
from .equilibrium import (
    GlissAllocationError,
    GlissInternalError,
    _error_buffer,
    _raise_for_status,
)
from ._stability_input import mode_integer

if TYPE_CHECKING:
    from .stability import StabilityProblem


class _EnergyTerms(ctypes.Structure):
    _fields_ = [
        ("struct_size", ctypes.c_size_t),
        ("field_line_bending", ctypes.c_double),
        ("magnetic_shear", ctypes.c_double),
        ("magnetic_compression", ctypes.c_double),
        ("pressure_drive", ctypes.c_double),
        ("plasma_compressibility", ctypes.c_double),
        ("potential_energy", ctypes.c_double),
        ("kinetic_energy", ctypes.c_double),
        ("rayleigh_quotient", ctypes.c_double),
        ("closure_error", ctypes.c_double),
        ("closure_tolerance", ctypes.c_double),
    ]


@dataclass(frozen=True)
class EnergyTerms:
    """Quadratic-form contributions for one displacement vector."""

    field_line_bending: float
    magnetic_shear: float
    magnetic_compression: float
    pressure_drive: float
    plasma_compressibility: float
    potential_energy: float
    kinetic_energy: float
    rayleigh_quotient: float
    closure_error: float
    closure_tolerance: float
    potential_form: str = "x.T @ K @ x"
    kinetic_form: str = "x.T @ M @ x"

    @property
    def components(self) -> Tuple[float, float, float, float, float]:
        """Terms in the order used by the native stiffness assembly."""
        return (
            self.field_line_bending,
            self.magnetic_shear,
            self.magnetic_compression,
            self.pressure_drive,
            self.plasma_compressibility,
        )


def _bind(library: Any) -> None:
    _require_symbols(
        library,
        ("gliss_stability_problem_energy",),
        "stability energy decomposition",
    )
    library.gliss_stability_problem_energy.argtypes = (
        ctypes.c_void_p,
        ctypes.c_int32,
        ctypes.c_size_t,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(_EnergyTerms),
        ctypes.c_void_p,
        ctypes.c_size_t,
    )
    library.gliss_stability_problem_energy.restype = ctypes.c_int


def diagnose_energy(
    problem: "StabilityProblem", parity_class: int, vector: Any
) -> EnergyTerms:
    """Evaluate energy terms through the native assembled problem."""
    problem._require_open()
    parity = mode_integer(parity_class, "parity_class")
    if parity not in (1, 2):
        raise ValueError("parity_class must be 1 or 2")
    count = problem._unknown_count(parity)
    values = _energy_vector(vector, count)
    _bind(problem._library)
    native = _EnergyTerms(struct_size=ctypes.sizeof(_EnergyTerms))
    error = _error_buffer()
    status = problem._library.gliss_stability_problem_energy(
        problem._handle,
        parity,
        count,
        values.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ctypes.byref(native),
        error,
        len(error),
    )
    _raise_for_status(status, error, "gliss_stability_problem_energy")
    return _validated_terms(native)


def _energy_vector(vector: Any, count: int) -> np.ndarray:
    try:
        try:
            raw = np.asarray(vector)
        except (TypeError, ValueError) as error:
            raise TypeError("energy vector must be an array of real numbers") from error
        if np.issubdtype(raw.dtype, np.bool_) or np.iscomplexobj(raw):
            raise TypeError("energy vector entries must be real numbers")
        try:
            values = np.asarray(raw, dtype=np.float64)
        except (TypeError, ValueError) as error:
            raise TypeError("energy vector entries must be real numbers") from error
        if values.ndim != 1:
            raise ValueError("energy vector must be one-dimensional")
        if values.size != count:
            raise ValueError(f"energy vector must contain {count} entries")
        if not np.all(np.isfinite(values)):
            raise ValueError("energy vector must contain only finite values")
        if not np.any(values):
            raise ValueError("energy vector must have positive kinetic norm")
        return np.ascontiguousarray(values)
    except MemoryError as error:
        raise GlissAllocationError(
            f"failed to allocate energy vector with {count} entries"
        ) from error


def _validated_terms(native: _EnergyTerms) -> EnergyTerms:
    result = EnergyTerms(
        native.field_line_bending,
        native.magnetic_shear,
        native.magnetic_compression,
        native.pressure_drive,
        native.plasma_compressibility,
        native.potential_energy,
        native.kinetic_energy,
        native.rayleigh_quotient,
        native.closure_error,
        native.closure_tolerance,
    )
    scalars = (*result.components, result.potential_energy, result.kinetic_energy,
               result.rayleigh_quotient, result.closure_error,
               result.closure_tolerance)
    if not all(math.isfinite(value) for value in scalars):
        raise GlissInternalError("GLISS returned nonfinite energy terms")
    if result.kinetic_energy <= 0.0:
        raise GlissInternalError("GLISS returned nonpositive kinetic energy")
    if result.closure_error < 0.0 or result.closure_tolerance < 0.0:
        raise GlissInternalError("GLISS returned an invalid energy closure bound")
    observed = abs(result.potential_energy - math.fsum(result.components))
    if observed > result.closure_tolerance:
        raise GlissInternalError("GLISS energy components do not close")
    if abs(observed - result.closure_error) > result.closure_tolerance:
        raise GlissInternalError("GLISS returned an inconsistent closure error")
    quotient = result.potential_energy / result.kinetic_energy
    quotient_tolerance = 64.0 * np.finfo(np.float64).eps * max(
        1.0, abs(quotient), abs(result.rayleigh_quotient)
    )
    if abs(quotient - result.rayleigh_quotient) > quotient_tolerance:
        raise GlissInternalError("GLISS returned an inconsistent Rayleigh quotient")
    for name, value in zip(
        ("field-line bending", "magnetic shear", "magnetic compression",
         "plasma compressibility"),
        (result.field_line_bending, result.magnetic_shear,
         result.magnetic_compression, result.plasma_compressibility),
    ):
        if value < -result.closure_tolerance:
            raise GlissInternalError(f"GLISS returned negative {name} energy")
    return result
