"""Independent physical-dilation oracle for the complete native spectral path.

For dimensionless a > 0, X_a = a X, Phi_a = a**2 Phi and chi_a = a**2 chi,
with pressure (Pa), magnetic-field strength (T), density (kg/m**3) and gamma
fixed as functions of the unchanged flux labels. Then B_a(a X) = B(X), and
both J cross B and grad(p) scale as a**-1. An exactly force-balanced baseline
therefore remains force balanced. These coarse analytic Solov'ev exports retain
their original equilibrium/discretization errors; dilation cannot certify them.

The imported Cartesian positions have units m, fluxes T m**2, metric m**2,
Jacobian m**3, second fundamental form m, covariant B components T m, and
contravariant angular B components T/m. Angles, s, rho, iota, beta, field
periods, winding and the Boozer/Cartesian frame are fixed. No external solver
or coordinate conversion is differentiated by this test.

At fixed radial/angular topology, fixed boundary and fixed density/gamma,
the complete spectrum obeys omega**2(a) = omega**2(1)/a**2. Traces over whole
eigenvalue clusters and each physical-energy component in an M-orthonormal
cluster basis obey the same law. This checks the actual import, primitive
reconstruction, mass/stiffness assembly and eigensolver against physical
similarity, without assuming coordinate-basis matrices scale uniformly.
"""

import shutil

import numpy as np
import pytest
from scipy.io import netcdf_file

import gliss

pytestmark = pytest.mark.native

# This is test-side physical dimensional analysis, not production derivative
# code. Listing all fields explicitly prevents silently freezing new inputs.
_POWERS = {
    "N_FP": 0, "winding": 0, "s": 0, "rho": 0, "m": 0, "n": 0,
    "p": 0, "iota": 0, "beta_avg": 0, "mod_B": 0,
    "Phi": 2, "chi": 2, "B_theta_avg": 1, "B_zeta_avg": 1,
    "xhat": 1, "yhat": 1, "zhat": 1, "Jac": 3,
    "g_tt": 2, "g_tz": 2, "g_zz": 2, "g_st": 2, "g_sz": 2,
    "II_tt": 1, "II_tz": 1, "II_zz": 1,
    "B_contra_t": -1, "B_contra_z": -1,
}


def _dilate(source, destination, scale, *, scale_flux=True):
    shutil.copyfile(source, destination)
    with netcdf_file(destination, "a", mmap=False) as file:
        for name, variable in file.variables.items():
            field = name[:-4] if name.endswith(("_mnc", "_mns")) else name
            power = _POWERS[field]
            if not scale_flux and field in ("Phi", "chi"):
                power = 0
            if power:
                variable[...] = variable.data.copy() * scale**power


def _evaluate(path, radial_cells):
    with gliss.Equilibrium(path) as equilibrium:
        with gliss.StabilityProblem(
            equilibrium, [(1, 1), (2, 1)], degree=2,
            angular_theta=64, angular_zeta=16, radial_cells=radial_cells,
        ) as problem:
            result = {}
            for parity in problem.parity_classes:
                spectrum = problem.solve_full_spectrum_class(parity)
                # The translated fixture uses a coupled representation with
                # a doublet. Summing the whole cluster removes basis choices.
                count = 2 if parity == 0 else 1
                energies = [
                    problem.energy(parity, vector)
                    for vector in spectrum.eigenvectors[-count:]
                ]
                components = np.sum([energy.components for energy in energies], axis=0)
                result[parity] = (
                    spectrum.eigenvalues,
                    float(np.sum(spectrum.eigenvalues[-count:])),
                    components,
                )
            return result


@pytest.mark.parametrize(
    ("fixture", "radial_cells"),
    [("solovev_q1.045.nc", 2), ("solovev_q1.035_shifted.nc", 4)],
)
def test_native_physical_dilation_spectrum_and_energy(
    native_library, test_data, tmp_path, fixture, radial_cells,
):
    source = test_data / fixture
    baseline = _evaluate(source, radial_cells)
    for scale in (0.7, 1.6):
        path = tmp_path / f"dilation-{scale}.nc"
        _dilate(source, path, scale)
        scaled = _evaluate(path, radial_cells)
        assert scaled.keys() == baseline.keys()
        for parity, (values, trace, components) in baseline.items():
            actual_values, actual_trace, actual_components = scaled[parity]
            np.testing.assert_allclose(
                actual_values * scale**2, values, rtol=2e-8, atol=1e-5,
            )
            assert actual_trace * scale**2 == pytest.approx(trace, rel=2e-9)
            np.testing.assert_allclose(
                actual_components * scale**2, components, rtol=2e-8, atol=1e-5,
            )
    # Freezing the fluxes under geometric dilation changes the physical B.
    # This wrong transformation must fail the oracle even though it imports.
    path = tmp_path / "frozen-flux.nc"
    _dilate(source, path, 1.6, scale_flux=False)
    wrong = _evaluate(path, radial_cells)
    for parity, (_, trace, _) in baseline.items():
        assert abs(wrong[parity][1] * 1.6**2 / trace - 1.0) > 0.1


def test_native_dilation_derivative_step_study(native_library, test_data, tmp_path):
    source = test_data / "solovev_q1.035_shifted.nc"
    baseline = _evaluate(source, 3)[0][1]
    # This exact derivative is a physical oracle, not a finite-difference
    # production derivative or an arbitrary GVEC design-variable gradient.
    derivative = -2.0 * baseline
    errors = []
    for step in (1e-2, 3e-3, 1e-3):
        values = []
        for sign in (-1, 1):
            path = tmp_path / f"dilation-{step}-{sign}.nc"
            _dilate(source, path, 1.0 + sign * step)
            values.append(_evaluate(path, 3)[0][1])
        difference = (values[1] - values[0]) / (2.0 * step)
        error = abs(difference / derivative - 1.0)
        errors.append(error)
        assert error < 3.0 * step**2 + 2e-8
    assert errors[-1] < 0.05 * errors[0]
