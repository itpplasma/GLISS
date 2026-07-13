"""Exact derivative actions for assembled fixed-boundary objectives."""

import ctypes
import math
import numbers
from typing import Any, TYPE_CHECKING

import numpy as np

from . import _require_symbols
from ._stability_input import mode_integer
from .energy import _coefficient_vector
from .equilibrium import (
    GlissInternalError,
    _empty_float64,
    _error_buffer,
    _raise_for_status,
)

if TYPE_CHECKING:
    from .stability import StabilityProblem


def _bind(library: Any) -> None:
    _require_symbols(
        library,
        ("gliss_stability_problem_rayleigh_vjp",),
        "Rayleigh derivative actions",
    )
    library.gliss_stability_problem_rayleigh_vjp.argtypes = (
        ctypes.c_void_p,
        ctypes.c_int32,
        ctypes.c_size_t,
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_double,
        ctypes.c_size_t,
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_void_p,
        ctypes.c_size_t,
    )
    library.gliss_stability_problem_rayleigh_vjp.restype = ctypes.c_int


def rayleigh_vjp(
    problem: "StabilityProblem",
    parity_class: int,
    vector: Any,
    cotangent: Any = 1.0,
) -> np.ndarray:
    """Reverse derivative of the Rayleigh quotient with respect to ``vector``."""
    parity, count, primal = _validated_primal(problem, parity_class, vector)
    weight = _cotangent(cotangent)
    return _native_vjp(problem, parity, count, primal, weight)


def rayleigh_jvp(
    problem: "StabilityProblem", parity_class: int, vector: Any, tangent: Any
) -> float:
    """Forward derivative of the Rayleigh quotient along ``tangent``."""
    parity, count, primal = _validated_primal(problem, parity_class, vector)
    direction = _coefficient_vector(tangent, count, "Rayleigh tangent", False)
    gradient = _native_vjp(problem, parity, count, primal, 1.0)
    return float(np.dot(gradient, direction))


def _validated_primal(
    problem: "StabilityProblem", parity_class: int, vector: Any
) -> tuple[int, int, np.ndarray]:
    problem._require_open()
    parity = mode_integer(parity_class, "parity_class")
    if parity not in (1, 2):
        raise ValueError("parity_class must be 1 or 2")
    count = problem._unknown_count(parity)
    primal = _coefficient_vector(vector, count, "Rayleigh vector", True)
    return parity, count, primal


def _cotangent(value: Any) -> float:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, numbers.Real):
        raise TypeError("cotangent must be a real number")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError("cotangent must be finite")
    return result


def _native_vjp(
    problem: "StabilityProblem",
    parity: int,
    count: int,
    primal: np.ndarray,
    cotangent: float,
) -> np.ndarray:
    _bind(problem._library)
    gradient = _empty_float64(count, f"Rayleigh gradient with {count} entries")
    error = _error_buffer()
    status = problem._library.gliss_stability_problem_rayleigh_vjp(
        problem._handle,
        parity,
        count,
        primal.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        cotangent,
        count,
        gradient.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        error,
        len(error),
    )
    _raise_for_status(status, error, "gliss_stability_problem_rayleigh_vjp")
    if not np.all(np.isfinite(gradient)):
        raise GlissInternalError("GLISS returned a nonfinite Rayleigh gradient")
    gradient.setflags(write=False)
    return gradient
