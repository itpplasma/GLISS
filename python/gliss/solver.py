"""Validated stopping controls for certified fixed-boundary solves."""

import ctypes
import math
import numbers
from dataclasses import dataclass
from typing import Any, Mapping


@dataclass(frozen=True)
class SolverTolerances:
    """Stopping tolerances and iteration limits for the block solver."""

    eigenvalue_relative: float = 1.0e-13
    residual_relative: float = 1.0e-12
    negative_bracket_relative: float = 1.0e-9
    negative_bracket_floor: float = 1.0e-3
    inverse_iteration_limit: int = 500
    bracket_iteration_limit: int = 200

    def __post_init__(self) -> None:
        for name in (
            "eigenvalue_relative",
            "residual_relative",
            "negative_bracket_relative",
            "negative_bracket_floor",
        ):
            object.__setattr__(self, name, _positive_real(getattr(self, name), name))
        for name in ("inverse_iteration_limit", "bracket_iteration_limit"):
            object.__setattr__(self, name, _iteration_limit(getattr(self, name), name))

    @classmethod
    def historical_defaults(cls) -> "SolverTolerances":
        """Return the stopping controls used before controls were public."""
        return cls()

    def to_dict(self) -> dict:
        """Return canonical finite JSON data for schema version 2."""
        return {
            "eigenvalue_relative": self.eigenvalue_relative,
            "residual_relative": self.residual_relative,
            "negative_bracket_relative": self.negative_bracket_relative,
            "negative_bracket_floor": self.negative_bracket_floor,
            "inverse_iteration_limit": self.inverse_iteration_limit,
            "bracket_iteration_limit": self.bracket_iteration_limit,
        }

    @classmethod
    def from_dict(cls, document: Mapping[str, Any]) -> "SolverTolerances":
        """Strictly validate solver controls from schema version 2."""
        if not isinstance(document, dict):
            raise ValueError("solver_tolerances must be an object")
        expected = set(cls.__dataclass_fields__)
        unknown = sorted(set(document) - expected)
        missing = sorted(expected - set(document))
        if unknown:
            raise ValueError(f"solver_tolerances has unknown field {unknown[0]!r}")
        if missing:
            raise ValueError(f"solver_tolerances is missing field {missing[0]!r}")
        try:
            return cls(**document)
        except (TypeError, ValueError) as error:
            raise ValueError(f"solver_tolerances: {error}") from error


class SolverTolerancesC(ctypes.Structure):
    """Private layout matching gliss_solver_tolerances."""

    _fields_ = [
        ("struct_size", ctypes.c_size_t),
        ("eigenvalue_relative", ctypes.c_double),
        ("residual_relative", ctypes.c_double),
        ("negative_bracket_relative", ctypes.c_double),
        ("negative_bracket_floor", ctypes.c_double),
        ("inverse_iteration_limit", ctypes.c_int32),
        ("bracket_iteration_limit", ctypes.c_int32),
    ]


def solver_tolerances_c(value: SolverTolerances) -> SolverTolerancesC:
    return SolverTolerancesC(
        ctypes.sizeof(SolverTolerancesC),
        value.eigenvalue_relative,
        value.residual_relative,
        value.negative_bracket_relative,
        value.negative_bracket_floor,
        value.inverse_iteration_limit,
        value.bracket_iteration_limit,
    )


def bind_solver_tolerances(library: Any) -> None:
    if not hasattr(library, "gliss_stability_problem_set_solver_tolerances"):
        return
    function = library.gliss_stability_problem_set_solver_tolerances
    function.argtypes = (
        ctypes.c_void_p,
        ctypes.POINTER(SolverTolerancesC),
        ctypes.c_void_p,
        ctypes.c_size_t,
    )
    function.restype = ctypes.c_int


def configure_solver_tolerances(
    library: Any, handle: ctypes.c_void_p, value: SolverTolerances
) -> None:
    if value == SolverTolerances.historical_defaults():
        return
    if not hasattr(library, "gliss_stability_problem_set_solver_tolerances"):
        raise OSError("custom solver tolerances require a matching native library")
    from .equilibrium import _error_buffer, _raise_for_status

    controls = solver_tolerances_c(value)
    error = _error_buffer()
    status = library.gliss_stability_problem_set_solver_tolerances(
        handle, ctypes.byref(controls), error, len(error)
    )
    _raise_for_status(
        status, error, "gliss_stability_problem_set_solver_tolerances"
    )


def _positive_real(value: Any, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, numbers.Real):
        raise TypeError(f"{name} must be a real number")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    if result <= 0.0:
        raise ValueError(f"{name} must be positive")
    return result


def _iteration_limit(value: Any, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, numbers.Integral):
        raise TypeError(f"{name} must be an integer")
    result = int(value)
    if result < 1:
        raise ValueError(f"{name} must be at least 1")
    if result > 2**31 - 1:
        raise ValueError(f"{name} exceeds the signed 32-bit range")
    return result
