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
        self.gliss_cas3d_marginality = FakeFunction(self.solve)

    def solve(
        self,
        equilibrium,
        mode_count,
        mode_m,
        mode_n,
        parity_class,
        radial_quadrature,
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
            radial_quadrature,
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
        native.radial_quadrature = radial_quadrature
        native.angular_theta = angular_theta
        native.angular_zeta = angular_zeta
        native.negative_count = 2
        native.lowest_eigenvalue = -5.3e-4 if solve_eigenpair else np.nan
        native.certificate = 2.0e-10 if solve_eigenpair else np.nan
        native.eigenpair_residual = 4.0e-13 if solve_eigenpair else np.nan
        native.force_balance_residual = 7.0e-5
        error.value = b""
        return 0


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

    assert library.arguments == (23, (3, 8), (-2, -7), 2, 1, 72, 48, 0)
    assert result.has_eigenpair is False
    assert result.field_periods == 5
    assert result.modes == ((3, -2), (8, -7))
    assert result.radial_surfaces == 48
    assert result.parity_class == 2
    assert result.radial_quadrature == "midpoint"
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
    assert "artificial" in result.normalization
    assert "not a physical growth rate" in result.interpretation
    assert result.boundary_condition == "fixed"
    assert result.coordinate_handedness == "left-handed"
    assert result.fourier_convention == "2*pi*(m*theta - n*zeta/N_T)"


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
        ("radial_quadrature", 1, TypeError, "string"),
        ("radial_quadrature", "gauss2", ValueError, "midpoint"),
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
