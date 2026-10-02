"""Native pressure-chain checks against independently perturbed NetCDF inputs."""

import ctypes
import shutil

import numpy as np
import pytest
from scipy.io import netcdf_file

import gliss
from gliss.pressure_derivatives import _bind, _pointer

pytestmark = pytest.mark.native

OPTIONS = dict(degree=1, angular_theta=64, angular_zeta=16, radial_cells=2)
MODES = [(1, 1)]


def _perturb(source, destination, direction, step):
    shutil.copyfile(source, destination)
    with netcdf_file(destination, "a", mmap=False) as file:
        file.variables["p"][:] = file.variables["p"].data.copy() + step * direction


def _trace(path, parity, start, stop):
    with gliss.Equilibrium(path) as equilibrium:
        with gliss.StabilityProblem(equilibrium, MODES, **OPTIONS) as problem:
            spectrum = problem.solve_full_spectrum_class(parity)
            return float(np.sum(spectrum.eigenvalues[start:stop]))


@pytest.mark.parametrize("coupled", [False, True])
def test_pressure_spectral_chain_matches_step_study(
    native_library, test_data, tmp_path, coupled
):
    source = test_data / (
        "solovev_q1.035_shifted.nc" if coupled else "solovev_q1.045.nc"
    )
    equilibrium = gliss.Equilibrium(source)
    with gliss.StabilityProblem(equilibrium, MODES, **OPTIONS) as problem:
        equilibrium.close()
        # Closing the original owner cannot invalidate the retained samples.
        nodes, pressure = problem.pressure_samples()
        assert pressure.size == 16 and not pressure.flags.writeable
        parity = 0 if coupled else 1
        spectrum = problem.solve_full_spectrum_class(parity)
        stop = spectrum.eigenvalues.size
        start = stop - (2 if coupled else 1)
        if coupled:
            with pytest.raises(ValueError, match="unresolved cluster boundary"):
                problem.spectral_pressure_sensitivity(parity, stop - 1, stop, gap=1e-4)
        sensitivity = problem.spectral_pressure_sensitivity(parity, start, stop, gap=1e-4)
        direction = pressure * (0.25 + np.sin(4.0 * np.pi * nodes))
        jvp = problem.spectral_pressure_jvp(parity, start, stop, direction, gap=1e-4)
        assert sensitivity.jvp(direction) == pytest.approx(jvp, rel=2e-9, abs=1e-7)
        assert np.dot(sensitivity.vjp(0.37), direction) == pytest.approx(0.37 * jvp)
        assert not sensitivity.gradient.flags.writeable
    # Neither oracle reads the analytic tangent assembly: both reread the
    # independently changed p samples and reassemble/solve a fresh problem.
    differences = []
    for index, step in enumerate((1e-3, 3e-4, 1e-4)):
        plus = tmp_path / f"plus-{index}.nc"
        minus = tmp_path / f"minus-{index}.nc"
        _perturb(source, plus, direction, step)
        _perturb(source, minus, direction, -step)
        difference = (_trace(plus, parity, start, stop) - _trace(minus, parity, start, stop)) / (2 * step)
        differences.append(difference)
    np.testing.assert_allclose(differences, jvp, rtol=2e-5, atol=1e-2)
    assert sensitivity.jvp(direction) == pytest.approx(jvp)
    with pytest.raises(RuntimeError, match="closed"):
        problem.pressure_samples()


def test_pressure_trace_is_invariant_under_basis_rotation(native_library, test_data):
    with gliss.Equilibrium(test_data / "solovev_q1.035_shifted.nc") as equilibrium:
        with gliss.StabilityProblem(equilibrium, MODES, **OPTIONS) as problem:
            spectrum = problem.solve_full_spectrum_class(0)
            basis = np.ascontiguousarray(spectrum.eigenvectors[-2:])
            angle = 0.731
            rotation = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
            rotated = np.ascontiguousarray(rotation @ basis)
            nodes, pressure = problem.pressure_samples()
            direction = np.ascontiguousarray(pressure * np.cos(3.0 * nodes))
            _bind(problem._library)
            values = []
            for vectors in (basis, rotated):
                output = ctypes.c_double()
                error = ctypes.create_string_buffer(512)
                status = problem._library.gliss_stability_problem_pressure_trace_jvp(
                    problem._handle, 0, vectors.shape[1], 2, _pointer(vectors),
                    pressure.size, _pointer(direction), ctypes.byref(output), error, len(error),
                )
                assert status == 0, error.value
                values.append(output.value)
            assert values[0] == pytest.approx(values[1], rel=2e-10, abs=1e-6)


def test_pressure_native_failures_preserve_output(native_library, test_data):
    with gliss.Equilibrium(test_data / "solovev_q1.045.nc") as equilibrium:
        with gliss.StabilityProblem(equilibrium, MODES, **OPTIONS) as problem:
            nodes, pressure = problem.pressure_samples()
            basis = np.ascontiguousarray(problem.solve_full_spectrum_class(1).eigenvectors[-1:])
            output = ctypes.c_double(123.0)
            error = ctypes.create_string_buffer(512)
            direction = np.ones(pressure.size)
            direction[0] = np.nan
            status = problem._library.gliss_stability_problem_pressure_trace_jvp(
                problem._handle, 1, basis.shape[1], 1, _pointer(basis),
                pressure.size, _pointer(direction), ctypes.byref(output), error, len(error),
            )
            assert status != 0 and output.value == 123.0
            gradient = np.full(pressure.size, 321.0)
            status = problem._library.gliss_stability_problem_pressure_trace_vjp(
                problem._handle, 1, basis.shape[1], 1, _pointer(basis),
                1.0, pressure.size - 1, _pointer(gradient), error, len(error),
            )
            assert status == 3
            np.testing.assert_array_equal(gradient, np.full(pressure.size, 321.0))
            for tangent in (np.ones(pressure.size - 1), np.full(pressure.size, np.inf)):
                with pytest.raises(ValueError):
                    problem.spectral_pressure_jvp(1, 0, basis.shape[1], tangent, gap=1e-4)


def test_pressure_owner_is_retained_by_each_problem(native_library, test_data):
    equilibrium = gliss.Equilibrium(test_data / "solovev_q1.045.nc")
    first = gliss.StabilityProblem(equilibrium, MODES, **OPTIONS)
    try:
        with gliss.StabilityProblem(equilibrium, MODES, **OPTIONS) as second:
            equilibrium.close()
            first.close()
            nodes, pressure = second.pressure_samples()
            spectrum = second.solve_full_spectrum_class(1)
            stop = spectrum.eigenvalues.size
            derivative = second.spectral_pressure_jvp(
                1, stop - 1, stop, np.ones(pressure.size), gap=1e-4,
            )
            assert np.isfinite(derivative)
        # The snapshot has independent storage after the last native owner exits.
        assert np.all(np.isfinite(nodes)) and np.all(pressure > 0)
    finally:
        first.close()
        equilibrium.close()


def test_pressure_derivatives_reject_free_boundary(native_library, test_data):
    with gliss.Equilibrium(test_data / "solovev_q1.045.nc") as equilibrium:
        with gliss.StabilityProblem(
            equilibrium, MODES, **OPTIONS, vacuum=gliss.VacuumModel((8, 4), 0.1)
        ) as problem:
            nodes, pressure = problem.pressure_samples()
            with pytest.raises(ValueError, match="fixed-boundary"):
                problem.spectral_pressure_jvp(1, 0, 1, pressure, gap=1e-4)
            with pytest.raises(ValueError, match="fixed-boundary"):
                problem.spectral_pressure_sensitivity(1, 0, 1, gap=1e-4)
