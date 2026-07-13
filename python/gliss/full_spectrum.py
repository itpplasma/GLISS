"""Opt-in dense fixed-boundary spectra."""

from __future__ import annotations

import ctypes
import math
from dataclasses import dataclass
from typing import Any, Optional, Tuple, TYPE_CHECKING

import numpy as np

from . import _require_symbols
from .equilibrium import (
    GlissCapacityError,
    GlissInternalError,
    _empty_float64,
    _error_buffer,
    _raise_for_status,
)
from ._stability_input import mode_integer as _mode_integer
from .stability import SpectrumResult

if TYPE_CHECKING:
    from .stability import StabilityProblem


def _bind_full_spectrum(library: Any) -> None:
    _require_symbols(
        library,
        ("gliss_stability_problem_full_spectrum",),
        "full stability spectrum",
    )
    library.gliss_stability_problem_full_spectrum.argtypes = (
        ctypes.c_void_p,
        ctypes.c_int32,
        ctypes.c_size_t,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_size_t,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_size_t),
        ctypes.POINTER(ctypes.c_size_t),
        ctypes.c_void_p,
        ctypes.c_size_t,
    )
    library.gliss_stability_problem_full_spectrum.restype = ctypes.c_int


@dataclass(frozen=True)
class FullSpectrumResult:
    """Every mass-normalized eigenpair for one parity class."""

    certified_lowest: SpectrumResult
    eigenvalues: np.ndarray
    rayleigh_quotients: np.ndarray
    residuals: np.ndarray
    resolutions: np.ndarray
    eigenvectors: np.ndarray

    @property
    def certified_index(self) -> Optional[int]:
        """Index of the independently certified active eigenpair."""
        if not self.certified_lowest.has_eigenvector:
            return None
        if self.certified_lowest.negative_count:
            return 0
        return self.certified_lowest.floor_count

    @property
    def normal(self) -> np.ndarray:
        """Fixed-edge normal coefficients for every eigenvector."""
        return self.eigenvectors[:, : self.certified_lowest.normal_unknowns]

    @property
    def eta(self) -> np.ndarray:
        """Tangential eta coefficients for every eigenvector."""
        start = self.certified_lowest.normal_unknowns
        return self.eigenvectors[:, start : start + self.certified_lowest.eta_unknowns]

    @property
    def mu(self) -> np.ndarray:
        """Compressional mu coefficients for every eigenvector."""
        start = (
            self.certified_lowest.normal_unknowns + self.certified_lowest.eta_unknowns
        )
        return self.eigenvectors[:, start:]


@dataclass(frozen=True)
class FullStabilityResult:
    """Complete spectra for both decoupled parity classes."""

    classes: Tuple[FullSpectrumResult, FullSpectrumResult]

    @property
    def lowest(self) -> FullSpectrumResult:
        """Class with the algebraically lowest certified active eigenvalue."""
        return min(
            self.classes,
            key=lambda result: result.certified_lowest.lowest_eigenvalue,
        )

    def write(self, path: Any) -> None:
        """Atomically write a deterministic full-spectrum container."""
        from .full_schema import write_full_result

        write_full_result(self, path)

    @classmethod
    def read(cls, path: Any) -> "FullStabilityResult":
        """Read and strictly validate a full-spectrum container."""
        from .full_schema import read_full_result

        return read_full_result(path)


def solve_full_spectrum(problem: StabilityProblem) -> FullStabilityResult:
    problem._require_open()
    return FullStabilityResult(
        (
            solve_full_spectrum_class(problem, 1),
            solve_full_spectrum_class(problem, 2),
        )
    )


def solve_full_spectrum_class(
    problem: StabilityProblem, parity_class: int
) -> FullSpectrumResult:
    problem._require_open()
    parity_class = _mode_integer(parity_class, "parity_class")
    if parity_class not in (1, 2):
        raise ValueError("parity_class must be 1 or 2")
    _bind_full_spectrum(problem._library)
    certified = problem.solve_class(parity_class)
    count = problem._unknown_count(parity_class)
    if count > math.isqrt(np.iinfo(np.intp).max):
        raise GlissCapacityError(
            f"full spectrum for {count} unknowns exceeds NumPy's index limit"
        )
    arrays = _allocate_full_spectrum(count)
    _call_full_spectrum(problem, parity_class, arrays)
    _validate_full_spectrum(certified, *arrays)
    for array in arrays:
        array.setflags(write=False)
    return FullSpectrumResult(certified, *arrays)


def _allocate_full_spectrum(count: int) -> Tuple[np.ndarray, ...]:
    size = count * count
    return (
        _empty_float64(count, f"{count} full-spectrum eigenvalues"),
        _empty_float64(count, f"{count} Rayleigh quotients"),
        _empty_float64(count, f"{count} eigenpair residuals"),
        _empty_float64(count, f"{count} eigenpair resolutions"),
        _empty_float64(
            (count, count),
            f"full-spectrum eigenvectors with {size} entries",
        ),
    )


def _call_full_spectrum(
    problem: StabilityProblem, parity_class: int, arrays: Tuple[np.ndarray, ...]
) -> None:
    eigenvalues, rayleigh_quotients, residuals, resolutions, eigenvectors = arrays
    count = eigenvalues.size
    eigenvalues_written = ctypes.c_size_t()
    eigenvectors_written = ctypes.c_size_t()
    error = _error_buffer()
    pointer = ctypes.POINTER(ctypes.c_double)
    status = problem._library.gliss_stability_problem_full_spectrum(
        problem._handle,
        parity_class,
        count,
        eigenvalues.ctypes.data_as(pointer),
        residuals.ctypes.data_as(pointer),
        resolutions.ctypes.data_as(pointer),
        rayleigh_quotients.ctypes.data_as(pointer),
        count * count,
        eigenvectors.ctypes.data_as(pointer),
        ctypes.byref(eigenvalues_written),
        ctypes.byref(eigenvectors_written),
        error,
        len(error),
    )
    _raise_for_status(status, error, "gliss_stability_problem_full_spectrum")
    if eigenvalues_written.value != count:
        raise GlissCapacityError(
            f"GLISS wrote {eigenvalues_written.value} eigenvalues; expected {count}"
        )
    if eigenvectors_written.value != count * count:
        raise GlissCapacityError(
            f"GLISS wrote {eigenvectors_written.value} eigenvector entries; "
            f"expected {count * count}"
        )


def _validate_full_spectrum(
    certified: SpectrumResult,
    eigenvalues: np.ndarray,
    rayleigh_quotients: np.ndarray,
    residuals: np.ndarray,
    resolutions: np.ndarray,
    eigenvectors: np.ndarray,
) -> None:
    if not np.all(np.isfinite(eigenvalues)):
        raise GlissInternalError("GLISS returned nonfinite full-spectrum values")
    if not np.all(np.isfinite(rayleigh_quotients)):
        raise GlissInternalError("GLISS returned nonfinite Rayleigh quotients")
    if not np.all(np.isfinite(residuals)) or np.any(residuals < 0.0):
        raise GlissInternalError("GLISS returned invalid full-spectrum residuals")
    if not np.all(np.isfinite(resolutions)) or np.any(resolutions < 0.0):
        raise GlissInternalError("GLISS returned invalid spectrum resolutions")
    if not np.all(np.isfinite(eigenvectors)):
        raise GlissInternalError("GLISS returned nonfinite full-spectrum vectors")
    if np.any(eigenvalues[1:] < eigenvalues[:-1]):
        raise GlissInternalError("GLISS returned an unsorted full spectrum")
    negative_count = np.count_nonzero(eigenvalues < -certified.zero_floor)
    floor_count = np.count_nonzero(np.abs(eigenvalues) <= certified.zero_floor)
    if negative_count != certified.negative_count:
        raise GlissInternalError("full-spectrum and certified negative counts differ")
    if floor_count != certified.floor_count:
        raise GlissInternalError("full-spectrum and certified floor counts differ")
    _validate_certified_pair(certified, eigenvalues)


def _validate_certified_pair(
    certified: SpectrumResult, eigenvalues: np.ndarray
) -> None:
    if not certified.has_eigenvector:
        return
    index = 0 if certified.negative_count else certified.floor_count
    if index >= eigenvalues.size:
        raise GlissInternalError(
            "certified eigenpair index lies outside the full spectrum"
        )
    tolerance = certified.certificate + 16.0 * np.finfo(np.float64).eps * max(
        1.0, abs(certified.lowest_eigenvalue)
    )
    if abs(eigenvalues[index] - certified.lowest_eigenvalue) > tolerance:
        raise GlissInternalError(
            "full spectrum disagrees with the certified active eigenvalue"
        )
