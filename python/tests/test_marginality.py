import ctypes
from pathlib import Path

import numpy as np
import pytest

import gliss
from gliss.marginality import _bind


class FakeFunction:
    def __init__(self, callback):
        self.callback = callback
        self.argtypes = None
        self.restype = None

    def __call__(self, *args):
        return self.callback(*args)


class FakeLibrary:
    def __init__(self):
        self.arguments = None
        self.coefficient_arguments = None
        self.gliss_cas3d_marginality = FakeFunction(self.solve)
        self.gliss_cas3d_phase_envelope = FakeFunction(self.solve_envelope)
        self.gliss_cas3d2mn_phase_envelope = FakeFunction(
            self.solve_envelope_coefficient
        )

    def solve(
        self,
        equilibrium,
        mode_count,
        mode_m,
        mode_n,
        parity_class,
        degree,
        angular_theta,
        angular_zeta,
        solve_eigenpair,
        result,
        error,
        error_capacity,
    ):
        self.arguments = (
            equilibrium.value,
            tuple(mode_m),
            tuple(mode_n),
            parity_class,
            degree,
            angular_theta,
            angular_zeta,
            solve_eigenpair,
        )
        native = result._obj
        native.has_eigenpair = solve_eigenpair
        native.field_periods = 5
        native.mode_count = mode_count
        native.radial_surfaces = 48
        native.parity_class = parity_class
        native.degree = degree
        native.angular_theta = angular_theta
        native.angular_zeta = angular_zeta
        native.negative_count = 2
        native.lowest_eigenvalue = -5.3e-4 if solve_eigenpair else np.nan
        native.certificate = 2.0e-10 if solve_eigenpair else np.nan
        native.eigenpair_residual = 4.0e-13 if solve_eigenpair else np.nan
        native.force_balance_residual = 7.0e-5
        error.value = b""
        return 0

    def solve_envelope(
        self,
        equilibrium,
        base_m,
        base_n,
        envelope_count,
        envelope_m,
        envelope_n,
        parity_class,
        degree,
        angular_theta,
        angular_zeta,
        solve_eigenpair,
        result,
        error,
        error_capacity,
    ):
        self.arguments = (
            equilibrium.value,
            (base_m, base_n),
            tuple(envelope_m),
            tuple(envelope_n),
            parity_class,
            degree,
            angular_theta,
            angular_zeta,
            solve_eigenpair,
        )
        native = result._obj
        native.has_eigenpair = solve_eigenpair
        native.field_periods = 5
        native.mode_count = 2 * envelope_count - 1
        native.radial_surfaces = 48
        native.parity_class = parity_class
        native.degree = degree
        native.angular_theta = angular_theta
        native.angular_zeta = angular_zeta
        native.negative_count = 3
        native.lowest_eigenvalue = -3.8e-4 if solve_eigenpair else np.nan
        native.certificate = 3.0e-10 if solve_eigenpair else np.nan
        native.eigenpair_residual = 5.0e-13 if solve_eigenpair else np.nan
        native.force_balance_residual = 7.0e-5
        error.value = b""
        return 0

    def solve_envelope_coefficient(self, *arguments):
        self.coefficient_arguments = arguments[10:14]
        base_arguments = arguments[:10] + arguments[14:]
        status = self.solve_envelope(*base_arguments)
        result = base_arguments[-3]._obj
        if result.has_eigenpair:
            result.lowest_eigenvalue = -0.371
        return status


def equilibrium(library):
    result = gliss.Equilibrium.__new__(gliss.Equilibrium)
    result.path = Path("w7x.nc")
    result._library = library
    result._handle = ctypes.c_void_p(23)
    return result


def test_cas3d_marginality_inertia_uses_explicit_general_modes():
    library = FakeLibrary()

    result = gliss.cas3d_marginality_inertia(
        equilibrium(library),
        modes=[(3, -2), (8, -7)],
        parity_class=2,
        angular_theta=72,
        angular_zeta=48,
    )

    assert library.arguments == (23, (3, 8), (-2, -7), 2, 2, 72, 48, 0)
    assert result.has_eigenpair is False
    assert result.field_periods == 5
    assert result.modes == ((3, -2), (8, -7))
    assert result.radial_surfaces == 48
    assert result.parity_class == 2
    assert result.degree == 2
    assert result.angular_resolution == (72, 48)
    assert result.negative_count == 2
    assert result.lowest_eigenvalue is None
    assert result.certificate is None
    assert result.eigenpair_residual is None
    assert result.force_balance_residual == pytest.approx(7.0e-5)


def test_solve_cas3d_marginality_labels_artificial_normalization():
    result = gliss.solve_cas3d_marginality(equilibrium(FakeLibrary()), modes=[(3, -2)])

    assert result.has_eigenpair is True
    assert result.lowest_eigenvalue == pytest.approx(-5.3e-4)
    assert result.certificate == pytest.approx(2.0e-10)
    assert result.eigenpair_residual == pytest.approx(4.0e-13)
    assert "compatible" in result.normalization
    assert "not a physical growth rate" in result.interpretation
    assert result.boundary_condition == "fixed"
    assert result.coordinate_handedness == "left-handed"
    assert result.fourier_convention == "2*pi*(m*theta - n*zeta/N_T)"


def test_solve_cas3d_phase_envelope_reports_labeled_input_count():
    library = FakeLibrary()

    result = gliss.solve_cas3d_phase_envelope(
        equilibrium(library),
        base_mode=(3, 2),
        envelope_modes=[(0, 0), (1, 0), (0, 1), (0, -1)],
        parity_class=2,
        angular_theta=72,
        angular_zeta=48,
    )

    assert library.arguments == (
        23,
        (3, 2),
        (0, 1, 0, 0),
        (0, 0, 1, -1),
        2,
        2,
        72,
        48,
        1,
    )
    assert result.base_mode == (3, 2)
    assert result.envelope_modes == ((0, 0), (1, 0), (0, 1), (0, -1))
    assert result.labeled_sideband_count == 7
    assert result.inertia_zero_floor == 1.0e-12
    assert result.negative_count == 3
    assert result.lowest_eigenvalue == pytest.approx(-3.8e-4)
    assert "physical modes" in result.normalization
    assert result.base_fourier_convention == "2*pi*(M*theta - N*zeta/N_T)"
    assert result.envelope_fourier_convention == "2*pi*(m*theta - n*zeta)"


def test_solve_cas3d_phase_envelope_selects_exact_coefficient_norm():
    library = FakeLibrary()
    result = gliss.solve_cas3d_phase_envelope(
        equilibrium(library),
        base_mode=(3, 2),
        envelope_modes=[(0, 0), (1, 0)],
        degree=1,
        angular_theta=36,
        angular_zeta=24,
        normalization="cas3d2mn_coefficient",
        coefficient_angular_resolution=(36, 24),
        reference_length=10.0,
        radial_quadrature="cas3d_midpoint",
    )

    assert result.lowest_eigenvalue == pytest.approx(-0.371)
    assert library.coefficient_arguments == (36, 24, 10.0, 2)
    assert "half identity" in result.normalization
    assert result.coefficient_angular_resolution == (36, 24)
    assert result.reference_length == 10.0
    assert result.radial_quadrature == "cas3d_midpoint"


@pytest.mark.parametrize(
    ("arguments", "exception", "match"),
    [
        ({"normalization": "unknown"}, ValueError, "one of"),
        ({"normalization": None}, TypeError, "string"),
        (
            {"normalization": "cas3d2mn_coefficient", "degree": 2},
            ValueError,
            "degree=1",
        ),
        (
            {"normalization": "cas3d2mn_coefficient", "degree": 1},
            TypeError,
            "must be a pair",
        ),
        (
            {
                "normalization": "cas3d2mn_coefficient",
                "degree": 1,
                "coefficient_angular_resolution": (36, 24),
            },
            TypeError,
            "real number",
        ),
        (
            {
                "normalization": "cas3d2mn_coefficient",
                "degree": 1,
                "coefficient_angular_resolution": (36, 24),
                "reference_length": float("nan"),
            },
            ValueError,
            "finite and positive",
        ),
        (
            {
                "normalization": "cas3d2mn_coefficient",
                "degree": 1,
                "coefficient_angular_resolution": (36, 24),
                "reference_length": 1.0e200,
            },
            ValueError,
            "normal binary64 range",
        ),
        (
            {
                "normalization": "cas3d2mn_coefficient",
                "degree": 1,
                "coefficient_angular_resolution": (36, 24),
                "reference_length": 1.0e-200,
            },
            ValueError,
            "normal binary64 range",
        ),
        (
            {"coefficient_angular_resolution": (36, 24)},
            ValueError,
            "requires cas3d2mn",
        ),
        (
            {"radial_quadrature": "cas3d_midpoint"},
            ValueError,
            "requires cas3d2mn",
        ),
        ({"radial_quadrature": "unknown"}, ValueError, "one of"),
    ],
)
def test_phase_envelope_rejects_invalid_normalization(arguments, exception, match):
    with pytest.raises(exception, match=match):
        gliss.solve_cas3d_phase_envelope(
            equilibrium(FakeLibrary()),
            base_mode=(3, 2),
            envelope_modes=[(0, 0)],
            **arguments,
        )


@pytest.mark.parametrize(
    ("keyword", "value", "exception", "match"),
    [
        ("base_mode", (-1, 2), ValueError, "nonnegative"),
        ("base_mode", (0, -1), ValueError, "axis mode"),
        ("base_mode", (True, 1), TypeError, "integer"),
        ("envelope_modes", [], ValueError, "must not be empty"),
        ("envelope_modes", [(1, 0)], ValueError, "begin with"),
        (
            "envelope_modes",
            [(0, 0), (0, 0)],
            ValueError,
            "duplicate envelope",
        ),
        ("envelope_modes", [(0, 0), (-1, 0)], ValueError, "nonnegative"),
    ],
)
def test_cas3d_phase_envelope_rejects_invalid_input(keyword, value, exception, match):
    arguments = {"base_mode": (3, 2), "envelope_modes": [(0, 0), (1, 0)]}
    arguments[keyword] = value

    with pytest.raises(exception, match=match):
        gliss.cas3d_phase_envelope_inertia(equilibrium(FakeLibrary()), **arguments)


@pytest.mark.parametrize(
    ("keyword", "value", "exception", "match"),
    [
        ("modes", [], ValueError, "must not be empty"),
        ("modes", [(1, 1), (1, 1)], ValueError, "duplicate"),
        ("modes", [(-1, 1)], ValueError, "nonnegative"),
        ("modes", [(0, -1)], ValueError, "axis mode"),
        ("modes", [(True, 1)], TypeError, "integer"),
        ("parity_class", True, TypeError, "integer"),
        ("parity_class", 0, ValueError, "1 or 2"),
        ("degree", True, TypeError, "integer"),
        ("degree", 0, ValueError, "between 1 and 4"),
        ("degree", 5, ValueError, "between 1 and 4"),
        ("angular_theta", 3, ValueError, "at least 4"),
        ("angular_theta", 2**31, ValueError, "32-bit"),
        ("angular_zeta", 4.5, TypeError, "integer"),
    ],
)
def test_cas3d_marginality_rejects_invalid_input(keyword, value, exception, match):
    arguments = {"modes": [(1, 1)]}
    arguments[keyword] = value

    with pytest.raises(exception, match=match):
        gliss.cas3d_marginality_inertia(equilibrium(FakeLibrary()), **arguments)


@pytest.mark.parametrize(
    "function",
    [gliss.cas3d_marginality_inertia, gliss.solve_cas3d_marginality],
)
def test_cas3d_marginality_rejects_closed_equilibrium(function):
    value = equilibrium(FakeLibrary())
    value._handle = ctypes.c_void_p()

    with pytest.raises(RuntimeError, match="Equilibrium is closed"):
        function(value, modes=[(1, 1)])


def test_cas3d_marginality_rejects_non_equilibrium():
    with pytest.raises(TypeError, match="gliss.Equilibrium"):
        gliss.cas3d_marginality_inertia(object(), modes=[(1, 1)])


def test_cas3d_marginality_rejects_invalid_native_result():
    library = FakeLibrary()

    def invalid(*arguments):
        result = arguments[9]._obj
        result.has_eigenpair = 0
        result.field_periods = 0
        return 0

    library.gliss_cas3d_marginality = FakeFunction(invalid)

    with pytest.raises(gliss.GlissInternalError, match="invalid CAS3D"):
        gliss.cas3d_marginality_inertia(equilibrium(library), modes=[(1, 1)])


def test_cas3d_marginality_propagates_native_error():
    library = FakeLibrary()

    def fail(*arguments):
        arguments[10].value = b"mode table aliases the angular quadrature"
        return 4

    library.gliss_cas3d_marginality = FakeFunction(fail)

    with pytest.raises(gliss.GlissArgumentError, match="aliases"):
        gliss.cas3d_marginality_inertia(equilibrium(library), modes=[(1, 1)])


def test_cas3d_marginality_binding_requires_symbol():
    with pytest.raises(OSError, match="gliss_cas3d_marginality"):
        _bind(object())


def test_cas3d_phase_envelope_binding_requires_symbol():
    with pytest.raises(OSError, match="gliss_cas3d_phase_envelope"):
        gliss.solve_cas3d_phase_envelope(
            equilibrium(object()),
            base_mode=(3, 2),
            envelope_modes=[(0, 0)],
        )
