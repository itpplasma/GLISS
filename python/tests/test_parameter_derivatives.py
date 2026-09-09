"""Independent small-pencil oracle for parameter spectral selection/actions."""

from types import SimpleNamespace

import numpy as np
import pytest
from scipy.linalg import eigh

from gliss.parameter_derivatives import (
    _admit_cluster,
    _subspace_gradient,
    spectral_parameter_sensitivity,
)


class Pencil:
    # Off-diagonal gamma perturbation mixes eigenvectors and splits the doublet.
    compression = np.array([[2.0, 1.0, 0.0], [1.0, 3.0, 0.0], [0.0, 0.0, 4.0]])
    base = np.diag([2.0, 2.0, 12.0]) - 1.5 * compression

    def __init__(self, rho=2.0, gamma=1.5):
        self.rho, self.gamma = rho, gamma
        self.k = self.base + gamma * self.compression

    def solve_full_spectrum_class(self, parity):
        values, vectors = eigh(self.k, self.rho * np.eye(3))
        return SimpleNamespace(
            eigenvalues=values,
            eigenvectors=vectors.T,
            resolutions=np.zeros(3),
            residuals=np.zeros(3),
            certified_lowest=SimpleNamespace(
                density_kg_m3=self.rho,
                adiabatic_index=self.gamma,
            ),
        )

    def energy(self, parity, x):
        mass = self.rho * x @ x
        return SimpleNamespace(
            kinetic_energy=mass,
            rayleigh_quotient=x @ self.k @ x / mass,
            plasma_compressibility=self.gamma * x @ self.compression @ x,
        )


def test_cluster_derivative_plateau_duality_scaling_and_rotation():
    p = Pencil()
    result = spectral_parameter_sensitivity(p, 1, 0, 2, gap=1.0)
    np.testing.assert_allclose(result.gradient, [-1.0, 2.5], atol=1e-14)
    tangent = np.array([0.7, -0.3])
    assert result.jvp(tangent) == pytest.approx(tangent @ result.vjp(2.0) / 2.0)
    for h in [1e-2, 3e-3, 1e-3]:
        vals = [
            sum(
                Pencil(2 + s * h * tangent[0], 1.5 + s * h * tangent[1])
                .solve_full_spectrum_class(1)
                .eigenvalues[:2]
            )
            for s in [-1, 1]
        ]
        assert (vals[1] - vals[0]) / (2 * h) == pytest.approx(
            result.jvp(tangent), rel=2e-5
        )
    vectors = p.solve_full_spectrum_class(1).eigenvectors[:2]
    rotation = np.array([[0.6, -0.8], [0.8, 0.6]])
    np.testing.assert_allclose(
        _subspace_gradient(p, 1, rotation @ vectors, (2.0, 1.5)),
        result.gradient,
        atol=1e-14,
    )
    scaled = spectral_parameter_sensitivity(Pencil(8.0), 1, 0, 2, gap=0.1)
    assert scaled.value == pytest.approx(result.value / 4.0)
    with pytest.raises(ValueError, match="boundary"):
        spectral_parameter_sensitivity(p, 1, 0, 1, gap=0.1)


def test_isolated_and_guards():
    result = spectral_parameter_sensitivity(Pencil(), 1, 2, 3, gap=1.0)
    np.testing.assert_allclose(result.gradient, [-3.0, 2.0])
    for h in [1e-2, 3e-3, 1e-3]:
        values = [
            Pencil(gamma=1.5 + s * h).solve_full_spectrum_class(1).eigenvalues[2]
            for s in [-1, 1]
        ]
        assert (values[1] - values[0]) / (2 * h) == pytest.approx(result.gradient[1])
    for tangent in ([1.0], [1.0, np.nan], [[1.0, 2.0]]):
        with pytest.raises(ValueError):
            result.jvp(tangent)
    with pytest.raises(ValueError):
        result.vjp(float("nan"))
    with pytest.raises(ValueError):
        spectral_parameter_sensitivity(Pencil(), 1, 0, 1, gap=0.0)
    with pytest.raises(ValueError, match="boundary"):
        _admit_cluster(np.array([0.0, 1.0]), np.array([0.6, 0.6]), 0, 1, 0.1)


def test_nonfinite_native_selection_and_energy_fail_closed():
    from gliss.equilibrium import GlissInternalError

    with pytest.raises(GlissInternalError, match="nonfinite"):
        _admit_cluster(np.array([0.0, np.nan]), np.zeros(2), 0, 1, 0.1)
    p = Pencil()
    p.compression = np.full((3, 3), np.nan)
    with pytest.raises(GlissInternalError, match="nonfinite"):
        spectral_parameter_sensitivity(p, 1, 0, 2, gap=1.0)
