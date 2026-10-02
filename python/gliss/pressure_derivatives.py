"""Pressure-sample spectral derivatives at fixed imported geometry.

Samples are immutable pressure values in Pa on the imported equilibrium's s
grid. Geometry, magnetic flux/current profiles, density, gamma, modes, radial
and angular quadrature, and resonance topology are fixed. This partial map
does not impose force balance or differentiate an external equilibrium solve.
"""

import ctypes
from dataclasses import dataclass

import numpy as np

from . import _require_symbols
from ._stability_input import mode_integer
from .derivatives import _cotangent
from .energy import _coefficient_vector
from .equilibrium import (
    GlissCapacityError, GlissInternalError, _empty_float64,
    _error_buffer, _raise_for_status,
)
from .parameter_derivatives import _admit_cluster, _finite

PRESSURE_CONTRACT_VERSION = 1


def _bind(library):
    names = (
        "gliss_stability_problem_pressure_samples",
        "gliss_stability_problem_pressure_trace_jvp",
        "gliss_stability_problem_pressure_trace_vjp",
    )
    _require_symbols(library, names, "fixed-geometry pressure derivatives")
    pointer = ctypes.POINTER(ctypes.c_double)
    library.gliss_stability_problem_pressure_samples.argtypes = (
        ctypes.c_void_p, ctypes.c_size_t, pointer, pointer,
        ctypes.POINTER(ctypes.c_size_t), ctypes.c_void_p, ctypes.c_size_t,
    )
    library.gliss_stability_problem_pressure_trace_jvp.argtypes = (
        ctypes.c_void_p, ctypes.c_int32, ctypes.c_size_t, ctypes.c_size_t,
        pointer, ctypes.c_size_t, pointer, pointer, ctypes.c_void_p, ctypes.c_size_t,
    )
    library.gliss_stability_problem_pressure_trace_vjp.argtypes = (
        ctypes.c_void_p, ctypes.c_int32, ctypes.c_size_t, ctypes.c_size_t,
        pointer, ctypes.c_double, ctypes.c_size_t, pointer,
        ctypes.c_void_p, ctypes.c_size_t,
    )
    for name in names:
        getattr(library, name).restype = ctypes.c_int


def pressure_samples(problem):
    """Return owned read-only arrays (s coordinates, pressure samples in Pa).

    The problem retains the imported data until close(), so its originating
    Equilibrium may already be closed. Returned arrays survive problem.close().
    """
    problem._require_open()
    _bind(problem._library)
    function = problem._library.gliss_stability_problem_pressure_samples
    written = ctypes.c_size_t()
    error = _error_buffer()
    status = function(problem._handle, 0, None, None, ctypes.byref(written), error, len(error))
    if status != 3:
        _raise_for_status(status, error, "gliss_stability_problem_pressure_samples")
        raise GlissInternalError("GLISS returned an invalid pressure sample query")
    count = written.value
    if count < 2 or count > np.iinfo(np.intp).max:
        raise GlissCapacityError("invalid pressure sample count")
    nodes = _empty_float64(count, "pressure sample coordinates")
    pressure = _empty_float64(count, "pressure samples")
    status = function(
        problem._handle, count, _pointer(nodes), _pointer(pressure),
        ctypes.byref(written), error, len(error),
    )
    _raise_for_status(status, error, "gliss_stability_problem_pressure_samples")
    if written.value != count:
        raise GlissInternalError("GLISS changed its pressure sample count")
    if not np.all(np.isfinite(nodes)) or not np.all(np.isfinite(pressure)):
        raise GlissInternalError("GLISS returned nonfinite pressure samples")
    if np.any(np.diff(nodes) <= 0):
        raise GlissInternalError("GLISS returned unordered pressure coordinates")
    nodes.setflags(write=False)
    pressure.setflags(write=False)
    return nodes, pressure


@dataclass(frozen=True)
class PressureSensitivity:
    """Derivative of a fixed-index spectral sum with respect to pressure samples.

    gradient has units s^-2 Pa^-1, and sample directions use Pa. Whole
    unresolved eigenvalue clusters must be selected; their basis rotations
    leave the trace unchanged. The exterior-gap gate is the numerical
    admission diagnostic used for material sensitivities, not a mesh or
    equilibrium convergence certificate. Constructing the gradient currently
    costs one exact tangent assembly per pressure sample.
    """

    value: float
    gradient: np.ndarray
    pressure_pa: np.ndarray
    surfaces_s: np.ndarray
    start: int
    stop: int
    gap: float
    contract_version: int = PRESSURE_CONTRACT_VERSION

    def jvp(self, tangent):
        direction = _coefficient_vector(
            tangent, self.gradient.size, "pressure tangent", False
        )
        return _finite(float(np.dot(self.gradient, direction)))

    def vjp(self, cotangent=1.0):
        gradient = self.gradient * _cotangent(cotangent)
        if not np.all(np.isfinite(gradient)):
            raise GlissInternalError("nonfinite pressure VJP")
        gradient.setflags(write=False)
        return gradient


def _selection(problem, parity_class, start, stop, gap):
    problem._require_open()
    if problem.boundary_condition != "fixed":
        raise ValueError("pressure derivatives require a fixed-boundary problem")
    parity = problem._parity_class(parity_class)
    start = mode_integer(start, "start")
    stop = mode_integer(stop, "stop")
    if start < 0 or stop <= start:
        raise ValueError("require 0 <= start < stop <= spectrum size")
    gap = _cotangent(gap)
    if gap <= 0:
        raise ValueError("gap must be positive")
    spectrum = problem.solve_full_spectrum_class(parity)
    values = spectrum.eigenvalues
    if stop > values.size:
        raise ValueError("require 0 <= start < stop <= spectrum size")
    exterior_gap = _admit_cluster(
        values, spectrum.residuals + spectrum.resolutions, start, stop, gap
    )
    nodes, pressure = pressure_samples(problem)
    basis = np.ascontiguousarray(spectrum.eigenvectors[start:stop])
    return parity, start, stop, exterior_gap, values, basis, nodes, pressure


def spectral_pressure_sensitivity(problem, parity_class, start, stop, *, gap):
    """Return a pressure gradient for an isolated eigenvalue or complete cluster.

    The pressure domain is strictly positive at samples and assembly points;
    zero-width resonances require a valid full sample-direction domain.
    All other imported data and the assembly topology remain fixed. This
    exact reverse action currently uses ns tangent assemblies. For a single
    direction, spectral_pressure_jvp uses one assembly instead.
    """
    parity, start, stop, gap, values, basis, nodes, pressure = _selection(
        problem, parity_class, start, stop, gap
    )
    gradient = _empty_float64(pressure.size, "pressure gradient")
    error = _error_buffer()
    status = problem._library.gliss_stability_problem_pressure_trace_vjp(
        problem._handle, parity, basis.shape[1], basis.shape[0], _pointer(basis),
        1.0, gradient.size, _pointer(gradient), error, len(error),
    )
    _raise_for_status(status, error, "gliss_stability_problem_pressure_trace_vjp")
    if not np.all(np.isfinite(gradient)):
        raise GlissInternalError("GLISS returned a nonfinite pressure gradient")
    gradient.setflags(write=False)
    return PressureSensitivity(
        _finite(float(np.sum(values[start:stop]))), gradient, pressure, nodes,
        start, stop, gap,
    )


def spectral_pressure_jvp(problem, parity_class, start, stop, tangent, *, gap):
    """Differentiate a spectral sum along a pressure direction using one assembly."""
    parity, _, _, _, _, basis, _, pressure = _selection(
        problem, parity_class, start, stop, gap
    )
    direction = _coefficient_vector(tangent, pressure.size, "pressure tangent", False)
    derivative = ctypes.c_double()
    error = _error_buffer()
    status = problem._library.gliss_stability_problem_pressure_trace_jvp(
        problem._handle, parity, basis.shape[1], basis.shape[0], _pointer(basis),
        pressure.size, _pointer(direction), ctypes.byref(derivative), error, len(error),
    )
    _raise_for_status(status, error, "gliss_stability_problem_pressure_trace_jvp")
    return _finite(derivative.value)


def _pointer(array):
    return array.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
