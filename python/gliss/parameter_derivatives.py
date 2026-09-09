"""Version 1 derivatives of scalar material parameters at fixed equilibrium.

Parameter order is (density_kg_m3, adiabatic_index), in kg/m^3 and dimensionless
units. Both are strictly positive. Geometry, pressure, equilibrium, topology,
quadrature, and parity are fixed. These are not equilibrium-design derivatives.
Production assembly is affine in gamma and linear in rho:
    dK/dgamma = K_compressibility/gamma, dM/drho = M/rho.
The native energy API returns quadratic forms without factors of one half.
"""

from dataclasses import dataclass

import numpy as np

from ._stability_input import mode_integer
from .derivatives import _cotangent
from .energy import _coefficient_vector
from .equilibrium import GlissInternalError

PARAMETER_CONTRACT_VERSION = 1
PARAMETER_NAMES = ("density_kg_m3", "adiabatic_index")
PARAMETER_UNITS = ("kg/m^3", "1")


@dataclass(frozen=True)
class ParameterSensitivity:
    """Local derivative of a fixed-index spectral sum, in parameter order.

    ``start:stop`` is frozen and includes whole unresolved eigenvalue clusters.
    Internal crossings are allowed; crossings at either boundary are rejected.
    ``gap`` is the smaller exterior gap after numerical diagnostic allowances.
    Backward residuals do not rigorously enclose eigenvalue errors for arbitrary
    mass matrices. Admission is neither a rigorous gap nor convergence certificate.
    """

    value: float
    gradient: np.ndarray
    start: int
    stop: int
    gap: float
    parameters: tuple[float, float]
    contract_version: int = PARAMETER_CONTRACT_VERSION

    def jvp(self, tangent):
        direction = _coefficient_vector(tangent, 2, "parameter tangent", False)
        return _finite(float(np.dot(self.gradient, direction)))

    def vjp(self, cotangent=1.0):
        result = self.gradient * _cotangent(cotangent)
        if not np.all(np.isfinite(result)):
            raise GlissInternalError("nonfinite parameter VJP")
        result.setflags(write=False)
        return result


def spectral_parameter_sensitivity(problem, parity_class, start, stop, *, gap):
    """Differentiate the sum of eigenvalues in the half-open range [start, stop).

    Set stop=start+1 for an isolated eigenvalue. ``gap`` is a required positive
    absolute separation in s^-2, fixed by the caller before evaluation. Boundary
    gaps must exceed it plus both neighboring residuals and resolutions. The
    full spectrum and its existing independent lowest-pair certificate are
    recomputed from the native problem; caller-supplied spectra are not accepted.
    The caller must preserve indices and dimension across optimization steps.
    The residual-based boundary check is a numerical admission diagnostic;
    rigorous gap enclosures require additional conditioning or inertia bounds.
    """
    start = mode_integer(start, "start")
    stop = mode_integer(stop, "stop")
    gap = _cotangent(gap)
    if gap <= 0:
        raise ValueError("gap must be positive")
    spectrum = problem.solve_full_spectrum_class(parity_class)
    values = spectrum.eigenvalues
    if start < 0 or stop <= start or stop > len(values):
        raise ValueError("require 0 <= start < stop <= spectrum size")
    uncertainty = spectrum.resolutions + spectrum.residuals
    exterior_gap = _admit_cluster(values, uncertainty, start, stop, gap)
    summary = spectrum.certified_lowest
    parameters = (summary.density_kg_m3, summary.adiabatic_index)
    if not all(np.isfinite(p) and p > 0 for p in parameters):
        raise GlissInternalError("invalid native material parameters")
    gradient = _subspace_gradient(
        problem, parity_class, spectrum.eigenvectors[start:stop], parameters
    )
    return ParameterSensitivity(
        _finite(float(np.sum(values[start:stop]))),
        gradient,
        start,
        stop,
        exterior_gap,
        parameters,
    )


def _admit_cluster(values, resolutions, start, stop, gap):
    if not np.all(np.isfinite(values)) or not np.all(np.isfinite(resolutions)):
        raise GlissInternalError("nonfinite spectral selection data")
    if np.any(resolutions < 0) or np.any(np.diff(values) < 0):
        raise GlissInternalError("invalid spectral selection data")
    exterior_gap = float("inf")
    for left, right in ((start - 1, start), (stop - 1, stop)):
        if left < 0 or right >= len(values):
            continue
        separation = (
            values[right] - values[left] - resolutions[right] - resolutions[left]
        )
        if separation <= gap:
            raise ValueError(
                "unresolved cluster boundary: enlarge the cluster or resolve its gap"
            )
        exterior_gap = min(exterior_gap, float(separation))
    return exterior_gap


def _subspace_gradient(problem, parity, vectors, parameters):
    # The native eigensolver supplies an M-orthonormal basis. The trace of each
    # projected derivative is invariant under any orthogonal basis rotation.
    rho, gamma = parameters
    gradient = np.zeros(2)
    for vector in vectors:
        energy = problem.energy(parity, vector)
        gradient += (
            -energy.rayleigh_quotient / rho,
            energy.plasma_compressibility / (gamma * energy.kinetic_energy),
        )
    if not np.all(np.isfinite(gradient)):
        raise GlissInternalError("nonfinite parameter gradient")
    gradient.setflags(write=False)
    return gradient


def _finite(value):
    if not np.isfinite(value):
        raise GlissInternalError("nonfinite parameter derivative")
    return value
