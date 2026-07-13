"""Fixed-boundary ideal-MHD stability problems backed by the GLISS C ABI."""

import ctypes
import math
import numbers
import operator
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence, Tuple, TYPE_CHECKING

import numpy as np

from . import _load_library, _require_symbols
from .equilibrium import (
    Equilibrium,
    GlissCapacityError,
    GlissInternalError,
    _error_buffer,
    _raise_for_status,
)

if TYPE_CHECKING:
    from .schema import RunManifest, StabilityConfiguration

_INT32_MIN = -(2**31)
_INT32_MAX = 2**31 - 1
_QUADRATURE = {"midpoint": 1, "gauss2": 2}


class _SpectrumSummary(ctypes.Structure):
    _fields_ = [
        ("struct_size", ctypes.c_size_t),
        ("has_chart_metric", ctypes.c_int32),
        ("has_eigenvector", ctypes.c_int32),
        ("field_periods", ctypes.c_int32),
        ("parity_class", ctypes.c_int32),
        ("radial_quadrature", ctypes.c_int32),
        ("angular_theta", ctypes.c_int32),
        ("angular_zeta", ctypes.c_int32),
        ("mode_count", ctypes.c_size_t),
        ("unknowns", ctypes.c_size_t),
        ("normal_unknowns", ctypes.c_size_t),
        ("eta_unknowns", ctypes.c_size_t),
        ("mu_unknowns", ctypes.c_size_t),
        ("negative_count", ctypes.c_size_t),
        ("floor_count", ctypes.c_size_t),
        ("adiabatic_index", ctypes.c_double),
        ("density_kg_m3", ctypes.c_double),
        ("zero_floor", ctypes.c_double),
        ("lowest_eigenvalue", ctypes.c_double),
        ("certificate", ctypes.c_double),
        ("eigenpair_residual", ctypes.c_double),
        ("eigenpair_resolution", ctypes.c_double),
        ("inertia_interval", ctypes.c_double),
    ]


def _bind(library: Any) -> None:
    _require_symbols(
        library,
        (
            "gliss_stability_problem_create",
            "gliss_stability_problem_destroy",
            "gliss_stability_problem_unknown_count",
            "gliss_stability_problem_solve_class",
        ),
        "stability problem",
    )
    library.gliss_stability_problem_create.argtypes = (
        ctypes.c_void_p,
        ctypes.c_double,
        ctypes.c_double,
        ctypes.c_double,
        ctypes.c_size_t,
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.c_int32,
        ctypes.POINTER(ctypes.c_void_p),
        ctypes.c_void_p,
        ctypes.c_size_t,
    )
    library.gliss_stability_problem_create.restype = ctypes.c_int
    library.gliss_stability_problem_destroy.argtypes = (
        ctypes.POINTER(ctypes.c_void_p),
        ctypes.c_void_p,
        ctypes.c_size_t,
    )
    library.gliss_stability_problem_destroy.restype = ctypes.c_int
    library.gliss_stability_problem_unknown_count.argtypes = (
        ctypes.c_void_p,
        ctypes.c_int32,
        ctypes.POINTER(ctypes.c_size_t),
        ctypes.c_void_p,
        ctypes.c_size_t,
    )
    library.gliss_stability_problem_unknown_count.restype = ctypes.c_int
    library.gliss_stability_problem_solve_class.argtypes = (
        ctypes.c_void_p,
        ctypes.c_int32,
        ctypes.c_size_t,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_size_t),
        ctypes.POINTER(_SpectrumSummary),
        ctypes.c_void_p,
        ctypes.c_size_t,
    )
    library.gliss_stability_problem_solve_class.restype = ctypes.c_int


def _real_parameter(value: Any, name: str, allow_zero: bool) -> float:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, numbers.Real):
        raise TypeError(f"{name} must be a real number")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    if allow_zero:
        if result < 0.0:
            raise ValueError(f"{name} must be nonnegative")
    elif result <= 0.0:
        raise ValueError(f"{name} must be positive")
    return result


def _mode_integer(value: Any, name: str) -> int:
    if isinstance(value, (bool, np.bool_)):
        raise TypeError(f"{name} must be an integer")
    try:
        result = operator.index(value)
    except TypeError as error:
        raise TypeError(f"{name} must be an integer") from error
    if not _INT32_MIN <= result <= _INT32_MAX:
        raise ValueError(f"{name} must fit a signed 32-bit integer")
    return result


def _validate_modes(modes: Sequence[Tuple[int, int]]) -> Tuple[Tuple[int, int], ...]:
    try:
        entries = tuple(modes)
    except TypeError as error:
        raise TypeError("modes must be a sequence of (m, n) pairs") from error
    if not entries:
        raise ValueError("modes must not be empty")
    result = []
    seen = set()
    for index, entry in enumerate(entries):
        if not isinstance(entry, (tuple, list)) or len(entry) != 2:
            raise TypeError(f"modes[{index}] must be an (m, n) pair")
        m_value = _mode_integer(entry[0], f"modes[{index}] poloidal mode")
        n_value = _mode_integer(entry[1], f"modes[{index}] toroidal mode")
        if m_value < 0:
            raise ValueError(f"modes[{index}] poloidal mode must be nonnegative")
        if m_value == 0 and n_value < 0:
            raise ValueError(f"modes[{index}] axis mode requires nonnegative n")
        pair = (m_value, n_value)
        if pair in seen:
            raise ValueError(f"duplicate mode {pair!r}")
        seen.add(pair)
        result.append(pair)
    return tuple(result)


@dataclass(frozen=True)
class SpectrumResult:
    """Certified lowest eigenpair for one stellarator-symmetry parity class."""

    parity_class: int
    field_periods: int
    modes: Tuple[Tuple[int, int], ...]
    radial_quadrature: str
    angular_resolution: Tuple[int, int]
    adiabatic_index: float
    density_kg_m3: float
    zero_floor: float
    negative_count: int
    floor_count: int
    lowest_eigenvalue: float
    certificate: float
    eigenpair_residual: float
    eigenpair_resolution: float
    inertia_interval: float
    eigenvector: np.ndarray
    normal_unknowns: int
    eta_unknowns: int
    mu_unknowns: int
    has_chart_metric: bool
    has_eigenvector: bool
    eigenvalue_unit: str = "s^-2"
    boundary_condition: str = "fixed"
    normalization: str = "x.T @ M @ x = 1"
    coordinate_handedness: str = "left-handed"
    fourier_convention: str = "2*pi*(m*theta - n*zeta/N_T)"

    @property
    def normal(self) -> np.ndarray:
        """Fixed-edge normal coefficients in dynamic-layout order."""
        return self.eigenvector[: self.normal_unknowns]

    @property
    def eta(self) -> np.ndarray:
        """Tangential eta coefficients in dynamic-layout order."""
        start = self.normal_unknowns
        return self.eigenvector[start : start + self.eta_unknowns]

    @property
    def mu(self) -> np.ndarray:
        """Compressional mu coefficients in dynamic-layout order."""
        start = self.normal_unknowns + self.eta_unknowns
        return self.eigenvector[start : start + self.mu_unknowns]


@dataclass(frozen=True)
class StabilityResult:
    """Certified results for both decoupled parity classes."""

    classes: Tuple[SpectrumResult, SpectrumResult]

    @property
    def lowest(self) -> SpectrumResult:
        """Class with the algebraically lowest computed eigenvalue."""
        return min(self.classes, key=lambda result: result.lowest_eigenvalue)

    def to_dict(self) -> dict:
        """Return the canonical version-1 result document."""
        from ._result_schema import stability_result_to_dict

        return stability_result_to_dict(self)

    def write(self, path: Any) -> None:
        """Atomically write the canonical JSON result document."""
        from ._result_schema import write_stability_result

        write_stability_result(self, path)

    @classmethod
    def read(cls, path: Any) -> "StabilityResult":
        """Read and strictly validate a JSON result document."""
        from ._result_schema import read_stability_result

        return read_stability_result(path)


class StabilityProblem:
    """Reusable assembled fixed-boundary ideal-MHD eigenproblem."""

    def __init__(
        self,
        equilibrium: Equilibrium,
        modes: Sequence[Tuple[int, int]],
        adiabatic_index: float = 5.0 / 3.0,
        density_kg_m3: float = 1.0,
        zero_floor: float = 1.0,
        radial_quadrature: str = "midpoint",
    ):
        if not isinstance(equilibrium, Equilibrium):
            raise TypeError("equilibrium must be a gliss.Equilibrium")
        if equilibrium.closed:
            raise RuntimeError("Equilibrium is closed")
        if not isinstance(radial_quadrature, str):
            raise TypeError("radial_quadrature must be a string")
        if radial_quadrature not in _QUADRATURE:
            raise ValueError("radial_quadrature must be 'midpoint' or 'gauss2'")
        self.modes = _validate_modes(modes)
        self.adiabatic_index = _real_parameter(
            adiabatic_index, "adiabatic_index", allow_zero=True
        )
        self.density_kg_m3 = _real_parameter(
            density_kg_m3, "density_kg_m3", allow_zero=False
        )
        self.zero_floor = _real_parameter(zero_floor, "zero_floor", allow_zero=False)
        if self.zero_floor > 0.125 * np.finfo(np.float64).max:
            raise ValueError("zero_floor is too large for spectrum certification")
        self.radial_quadrature = radial_quadrature
        self.equilibrium_path = Path(equilibrium.path)
        self._library = _load_library()
        _bind(self._library)
        self._handle = ctypes.c_void_p()
        self._create(equilibrium)

    def _create(self, equilibrium: Equilibrium) -> None:
        count = len(self.modes)
        integers = ctypes.c_int32 * count
        mode_m = integers(*(mode[0] for mode in self.modes))
        mode_n = integers(*(mode[1] for mode in self.modes))
        error = _error_buffer()
        status = self._library.gliss_stability_problem_create(
            equilibrium._handle,
            self.adiabatic_index,
            self.density_kg_m3,
            self.zero_floor,
            count,
            mode_m,
            mode_n,
            _QUADRATURE[self.radial_quadrature],
            ctypes.byref(self._handle),
            error,
            len(error),
        )
        if status != 0 and self._handle.value is not None:
            self._library.gliss_stability_problem_destroy(
                ctypes.byref(self._handle), None, 0
            )
        _raise_for_status(status, error, "gliss_stability_problem_create")
        if self._handle.value is None:
            raise GlissInternalError("GLISS returned a null stability problem handle")

    @property
    def closed(self) -> bool:
        """Whether the assembled native problem has been released."""
        return self._handle.value is None

    def close(self) -> None:
        """Release the assembled native problem; repeated calls are safe."""
        if self.closed:
            return
        error = _error_buffer()
        status = self._library.gliss_stability_problem_destroy(
            ctypes.byref(self._handle), error, len(error)
        )
        _raise_for_status(status, error, "gliss_stability_problem_destroy")
        if not self.closed:
            raise GlissInternalError("GLISS did not clear the stability problem handle")

    def __enter__(self) -> "StabilityProblem":
        self._require_open()
        return self

    def __exit__(self, exc_type: Any, exc_value: Any, traceback: Any) -> None:
        self.close()

    def __repr__(self) -> str:
        state = "closed" if self.closed else "open"
        return (
            f"<gliss.StabilityProblem(path={self.equilibrium_path!r}, "
            f"modes={len(self.modes)}, state={state!r})>"
        )

    @property
    def configuration(self) -> "StabilityConfiguration":
        """Immutable serializable inputs for this assembled problem."""
        from .schema import StabilityConfiguration

        return StabilityConfiguration(
            self.modes,
            self.adiabatic_index,
            self.density_kg_m3,
            self.zero_floor,
            self.radial_quadrature,
        )

    def write_manifest(self, path: Any, result: StabilityResult) -> "RunManifest":
        """Write a portable manifest for a result produced by this problem."""
        from .schema import write_run_manifest

        return write_run_manifest(
            path, self.equilibrium_path, self.configuration, result
        )

    def solve(self) -> StabilityResult:
        """Solve and certify the lowest eigenpair in both parity classes."""
        self._require_open()
        return StabilityResult((self.solve_class(1), self.solve_class(2)))

    def solve_class(self, parity_class: int) -> SpectrumResult:
        """Solve and certify one parity class, numbered 1 or 2."""
        self._require_open()
        parity_class = _mode_integer(parity_class, "parity_class")
        if parity_class not in (1, 2):
            raise ValueError("parity_class must be 1 or 2")
        count = self._unknown_count(parity_class)
        vector = np.empty(count, dtype=np.float64)
        written = ctypes.c_size_t()
        summary = _SpectrumSummary(struct_size=ctypes.sizeof(_SpectrumSummary))
        error = _error_buffer()
        status = self._library.gliss_stability_problem_solve_class(
            self._handle,
            parity_class,
            count,
            vector.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            ctypes.byref(written),
            ctypes.byref(summary),
            error,
            len(error),
        )
        _raise_for_status(status, error, "gliss_stability_problem_solve_class")
        vector = self._validate_native_result(
            parity_class, count, written.value, summary, vector
        )
        vector.setflags(write=False)
        return self._result_from_summary(summary, vector)

    def _unknown_count(self, parity_class: int) -> int:
        count = ctypes.c_size_t()
        error = _error_buffer()
        status = self._library.gliss_stability_problem_unknown_count(
            self._handle, parity_class, ctypes.byref(count), error, len(error)
        )
        _raise_for_status(status, error, "gliss_stability_problem_unknown_count")
        if count.value > np.iinfo(np.intp).max:
            raise GlissCapacityError(
                f"stability problem size {count.value} exceeds NumPy's index limit"
            )
        return count.value

    def _validate_native_result(
        self,
        parity_class: int,
        expected: int,
        written: int,
        summary: _SpectrumSummary,
        vector: np.ndarray,
    ) -> np.ndarray:
        if summary.unknowns != expected:
            raise GlissCapacityError(
                f"GLISS returned {summary.unknowns} unknowns; expected {expected}"
            )
        if summary.parity_class != parity_class:
            raise GlissInternalError(
                f"GLISS returned parity class {summary.parity_class}; "
                f"expected {parity_class}"
            )
        if summary.has_eigenvector not in (0, 1):
            raise GlissInternalError(
                f"GLISS returned invalid has_eigenvector={summary.has_eigenvector}"
            )
        expected_written = expected if summary.has_eigenvector else 0
        if written != expected_written:
            raise GlissCapacityError(
                f"GLISS wrote {written} vector entries; expected {expected_written}"
            )
        components = (
            summary.normal_unknowns + summary.eta_unknowns + summary.mu_unknowns
        )
        if components != expected:
            raise GlissInternalError(
                f"GLISS component sizes total {components}; expected {expected}"
            )
        if not summary.has_eigenvector:
            return np.empty(0, dtype=np.float64)
        return vector

    def _result_from_summary(
        self, summary: _SpectrumSummary, vector: np.ndarray
    ) -> SpectrumResult:
        quadrature = {value: key for key, value in _QUADRATURE.items()}
        if summary.radial_quadrature not in quadrature:
            raise GlissInternalError(
                f"GLISS returned unknown radial quadrature {summary.radial_quadrature}"
            )
        return SpectrumResult(
            parity_class=summary.parity_class,
            field_periods=summary.field_periods,
            modes=self.modes,
            radial_quadrature=quadrature[summary.radial_quadrature],
            angular_resolution=(summary.angular_theta, summary.angular_zeta),
            adiabatic_index=summary.adiabatic_index,
            density_kg_m3=summary.density_kg_m3,
            zero_floor=summary.zero_floor,
            negative_count=summary.negative_count,
            floor_count=summary.floor_count,
            lowest_eigenvalue=summary.lowest_eigenvalue,
            certificate=summary.certificate,
            eigenpair_residual=summary.eigenpair_residual,
            eigenpair_resolution=summary.eigenpair_resolution,
            inertia_interval=summary.inertia_interval,
            eigenvector=vector,
            normal_unknowns=summary.normal_unknowns,
            eta_unknowns=summary.eta_unknowns,
            mu_unknowns=summary.mu_unknowns,
            has_chart_metric=bool(summary.has_chart_metric),
            has_eigenvector=bool(summary.has_eigenvector),
        )

    def _require_open(self) -> None:
        if self.closed:
            raise RuntimeError("StabilityProblem is closed")
