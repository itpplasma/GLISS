import numpy as np
import pytest

import gliss
from gliss.stability import SpectrumResult, StabilityProblem, _bind


class FakeFunction:
    def __init__(self, callback):
        self.callback = callback
        self.argtypes = None
        self.restype = None

    def __call__(self, *args):
        return self.callback(*args)


class FakeLibrary:
    def __init__(self):
        self.problem_creates = 0
        self.problem_destroys = 0
        self.has_eigenvector = True
        self.gliss_equilibrium_create = FakeFunction(self.equilibrium_create)
        self.gliss_equilibrium_destroy = FakeFunction(self.equilibrium_destroy)
        self.gliss_equilibrium_surface_count = FakeFunction(self.surface_count)
        self.gliss_mercier_profile_context = FakeFunction(self.mercier_profile)
        self.gliss_stability_problem_create = FakeFunction(self.problem_create)
        self.gliss_stability_problem_destroy = FakeFunction(self.problem_destroy)
        self.gliss_stability_problem_unknown_count = FakeFunction(
            self.problem_unknown_count
        )
        self.gliss_stability_problem_solve_class = FakeFunction(
            self.problem_solve_class
        )

    def equilibrium_create(self, path, length, handle, error, error_capacity):
        handle._obj.value = 10
        return 0

    def equilibrium_destroy(self, handle, error, error_capacity):
        handle._obj.value = None
        return 0

    def surface_count(self, handle, surfaces, error, error_capacity):
        surfaces._obj.value = 3
        return 0

    def mercier_profile(self, *args):
        raise AssertionError("Mercier was not requested")

    def problem_create(
        self,
        equilibrium,
        gamma,
        density,
        floor,
        mode_count,
        mode_m,
        mode_n,
        radial_quadrature,
        handle,
        error,
        error_capacity,
    ):
        self.problem_creates += 1
        handle._obj.value = 100 + self.problem_creates
        error.value = b""
        return 0

    def problem_destroy(self, handle, error, error_capacity):
        self.problem_destroys += 1
        handle._obj.value = None
        if error is not None:
            error.value = b""
        return 0

    def problem_unknown_count(
        self, handle, parity_class, unknowns, error, error_capacity
    ):
        unknowns._obj.value = 6
        error.value = b""
        return 0

    def problem_solve_class(
        self,
        handle,
        parity_class,
        capacity,
        vector,
        written,
        summary,
        error,
        error_capacity,
    ):
        if self.has_eigenvector:
            np.ctypeslib.as_array(vector, shape=(capacity,))[:] = np.arange(capacity)
        written._obj.value = capacity if self.has_eigenvector else 0
        values = summary._obj
        values.has_chart_metric = 1
        values.has_eigenvector = int(self.has_eigenvector)
        values.field_periods = 3
        values.parity_class = parity_class
        values.radial_quadrature = 1
        values.angular_theta = 64
        values.angular_zeta = 64
        values.mode_count = 2
        values.unknowns = capacity
        values.normal_unknowns = 2
        values.eta_unknowns = 2
        values.mu_unknowns = 2
        values.negative_count = 1
        values.floor_count = 3
        values.adiabatic_index = 5.0 / 3.0
        values.density_kg_m3 = 2.0
        values.zero_floor = 1.0
        values.lowest_eigenvalue = -4.0 - parity_class
        values.eigenpair_residual = 1.0e-10
        values.eigenpair_resolution = 2.0e-10
        values.inertia_interval = 3.0e-10
        values.certificate = 6.0e-10
        error.value = b""
        return 0


@pytest.fixture
def contexts(monkeypatch, tmp_path):
    library = FakeLibrary()
    monkeypatch.setattr("gliss.equilibrium._load_library", lambda: library)
    monkeypatch.setattr("gliss.stability._load_library", lambda: library)
    export = tmp_path / "equilibrium.nc"
    export.touch()
    equilibrium = gliss.Equilibrium(export)
    yield library, equilibrium
    equilibrium.close()


def test_stability_problem_lifecycle_and_results(contexts):
    library, equilibrium = contexts
    problem = StabilityProblem(
        equilibrium,
        modes=[(1, 1), (2, 1)],
        adiabatic_index=5.0 / 3.0,
        density_kg_m3=2.0,
        zero_floor=1.0,
    )
    assert repr(problem).startswith("<gliss.StabilityProblem(")
    equilibrium.close()

    result = problem.solve()
    assert len(result.classes) == 2
    assert all(isinstance(item, SpectrumResult) for item in result.classes)
    first = result.classes[0]
    assert first.parity_class == 1
    assert first.eigenvalue_unit == "s^-2"
    assert first.boundary_condition == "fixed"
    assert first.normalization == "x.T @ M @ x = 1"
    assert first.coordinate_handedness == "left-handed"
    assert first.fourier_convention == "2*pi*(m*theta - n*zeta/N_T)"
    assert first.has_eigenvector
    assert not first.eigenvector.flags.writeable
    np.testing.assert_array_equal(first.normal, [0.0, 1.0])
    np.testing.assert_array_equal(first.eta, [2.0, 3.0])
    np.testing.assert_array_equal(first.mu, [4.0, 5.0])

    problem.close()
    problem.close()
    assert problem.closed
    assert library.problem_creates == 1
    assert library.problem_destroys == 1
    with pytest.raises(RuntimeError, match="closed"):
        problem.solve()


@pytest.mark.parametrize(
    ("keyword", "value", "match"),
    [
        ("modes", [], "modes"),
        ("modes", [(1, 1), (1, 1)], "duplicate"),
        ("modes", [(-1, 1)], "poloidal"),
        ("modes", [(0, -1)], "axis"),
        ("modes", [(True, 1)], "integer"),
        ("modes", [(2**31, 1)], "32-bit"),
        ("adiabatic_index", float("nan"), "finite"),
        ("adiabatic_index", "5/3", "real number"),
        ("density_kg_m3", 0.0, "positive"),
        ("zero_floor", -1.0, "positive"),
        ("zero_floor", np.finfo(np.float64).max, "too large"),
        ("radial_quadrature", "unknown", "midpoint"),
        ("radial_quadrature", 1, "string"),
    ],
)
def test_stability_problem_rejects_invalid_inputs(contexts, keyword, value, match):
    _, equilibrium = contexts
    arguments = {"modes": [(1, 1)]}
    arguments[keyword] = value
    with pytest.raises((TypeError, ValueError), match=match):
        StabilityProblem(equilibrium, **arguments)


def test_stability_problem_requires_open_equilibrium(contexts):
    _, equilibrium = contexts
    equilibrium.close()
    with pytest.raises(RuntimeError, match="closed"):
        StabilityProblem(equilibrium, modes=[(1, 1)])


def test_stability_problem_requires_equilibrium_object(tmp_path):
    export = tmp_path / "equilibrium.nc"
    export.touch()
    with pytest.raises(TypeError, match="gliss.Equilibrium"):
        StabilityProblem(export, modes=[(1, 1)])


def test_stability_problem_failed_partial_create_is_cleaned_up(contexts):
    library, equilibrium = contexts

    def fail_create(*args):
        handle = args[8]
        error = args[9]
        handle._obj.value = 1234
        error.value = b"assembly failed"
        return 2

    library.gliss_stability_problem_create = FakeFunction(fail_create)
    with pytest.raises(gliss.GlissComputationError, match="assembly failed"):
        StabilityProblem(equilibrium, modes=[(1, 1)])
    assert library.problem_destroys == 1


def test_stability_problem_rejects_invalid_parity(contexts):
    _, equilibrium = contexts
    with StabilityProblem(equilibrium, modes=[(1, 1)]) as problem:
        with pytest.raises(TypeError, match="integer"):
            problem.solve_class(True)
        with pytest.raises(ValueError, match="1 or 2"):
            problem.solve_class(3)


def test_stability_problem_represents_fully_floored_result(contexts):
    library, equilibrium = contexts
    library.has_eigenvector = False
    with StabilityProblem(equilibrium, modes=[(1, 1), (2, 1)]) as problem:
        result = problem.solve_class(1)
    assert not result.has_eigenvector
    assert result.eigenvector.shape == (0,)
    assert result.normal.shape == (0,)
    assert not result.eigenvector.flags.writeable


def test_stability_problem_is_public():
    assert gliss.StabilityProblem is StabilityProblem
    assert gliss.SpectrumResult is SpectrumResult


def test_stability_problem_exposes_configuration_and_manifest(
    contexts, tmp_path, monkeypatch
):
    _, equilibrium = contexts
    monkeypatch.setattr("gliss.schema._native_version", lambda: "0.0.1")
    path = tmp_path / "run.json"
    with StabilityProblem(
        equilibrium,
        modes=[(1, 1), (2, 1)],
        density_kg_m3=2.0,
    ) as problem:
        result = problem.solve()
        assert problem.configuration == gliss.StabilityConfiguration(
            modes=((1, 1), (2, 1)), density_kg_m3=2.0
        )
        manifest = problem.write_manifest(path, result)

    assert path.is_file()
    assert gliss.RunManifest.read(path) == manifest


def test_stability_bind_reports_missing_native_symbols():
    with pytest.raises(OSError, match="stability problem.*matching"):
        _bind(object())
