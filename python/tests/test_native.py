"""Native checks of libgliss_c against independent oracles.

Every test here loads the real library; none substitutes a fake. The
oracles are the GPEC/DCON Newcomb results for the public Solov'ev fixtures
(benchmarks/solovev) and exact identities of the Rayleigh quotient.
"""

import numpy as np
import pytest

import gliss

pytestmark = pytest.mark.native


def test_native_solovev_mercier_sign_matches_dcon(native_library, test_data):
    # Independent GPEC/DCON Newcomb runs of this analytic Solov'ev family near
    # q0=1.04 find the core Mercier-unstable (D_I > 0) for s below about 0.16
    # and stable outside.  Positive D_Mercier is stable in GLISS.
    path = test_data / "solovev_q1.035.nc"
    s, d_mercier = gliss.mercier_profile(path)
    assert d_mercier[0] < 0.0
    assert np.all(d_mercier[s > 0.25] > 0.0)
    crossing = s[np.argmax(d_mercier > 0.0)]
    assert 0.1 < crossing < 0.3
    assert gliss.mercier_objective(path) == pytest.approx(-d_mercier.min())
    assert gliss.mercier_objective(path) > 0.0


def test_native_gvec_export_reports_left_handed_chart(native_library, test_data):
    with gliss.Equilibrium(test_data / "solovev_q1.035.nc") as equilibrium:
        assert equilibrium.coordinate_handedness == "left-handed"


@pytest.mark.parametrize(("name", "unstable"), [("1.035", True), ("1.045", False)])
def test_native_solovev_n1_stability_matches_dcon(
    native_library, test_data, name, unstable
):
    # GPEC/DCON (tolerance 1e-8) puts the n=1 fixed-boundary marginal point of
    # this family between q0=1.0391 (unstable) and q0=1.0398 (stable). The
    # conforming axis space bounds eigenvalues from above; degree 3 resolves
    # the near-marginal kink on these 16-surface exports.
    with gliss.Equilibrium(test_data / f"solovev_q{name}.nc") as equilibrium:
        result = gliss.solve_axisymmetric(equilibrium, poloidal_max=6, degree=3)
    assert result.negative_count == int(unstable)
    assert (result.lowest_eigenvalue < 0.0) == unstable
    assert result.certificate < abs(result.lowest_eigenvalue)
    assert 0.0 < result.inertia_zero_floor < abs(result.lowest_eigenvalue)


def test_native_rayleigh_gradient_matches_finite_differences(
    native_library, test_data
):
    # The Rayleigh quotient R(x) = x.K.x / x.M.x has gradient
    # 2 (K x - R M x) / x.M.x; it matches central differences of the native
    # energy for any x and vanishes at the lowest eigenvector.
    generator = np.random.default_rng(7)
    with gliss.Equilibrium(test_data / "solovev_q1.045.nc") as equilibrium:
        with gliss.StabilityProblem(
            equilibrium,
            modes=[(0, 1), (1, -1), (1, 1)],
            degree=1,
            angular_theta=32,
            angular_zeta=8,
        ) as problem:
            size = problem._unknown_count(1)
            vector = generator.standard_normal(size)
            tangent = generator.standard_normal(size)
            gradient = problem.rayleigh_vjp(1, vector)
            step = 1.0e-6 * np.linalg.norm(vector) / np.linalg.norm(tangent)
            plus = problem.energy(1, vector + step * tangent).rayleigh_quotient
            minus = problem.energy(1, vector - step * tangent).rayleigh_quotient
            difference = (plus - minus) / (2.0 * step)
            assert problem.rayleigh_jvp(1, vector, tangent) == pytest.approx(
                difference, rel=1.0e-6
            )
            assert np.dot(gradient, tangent) == pytest.approx(difference, rel=1.0e-6)

            lowest = problem.solve().classes[0]
            stationary = problem.rayleigh_vjp(1, lowest.eigenvector)
            scale = np.linalg.norm(gradient) * np.linalg.norm(vector)
            assert np.linalg.norm(stationary) * np.linalg.norm(
                lowest.eigenvector
            ) < 1.0e-6 * scale
