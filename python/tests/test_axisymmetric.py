import ctypes
from pathlib import Path

import numpy as np
import pytest

import gliss
from gliss.axisymmetric import _bind


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
        self.gliss_axisymmetric_spectrum = FakeFunction(self.solve)

    def solve(
        self,
        equilibrium,
        toroidal_mode,
        poloidal_max,
        radial_quadrature,
        solve_eigenpair,
        result,
        error,
        error_capacity,
    ):
        self.arguments = (
            equilibrium.value,
            toroidal_mode,
            poloidal_max,
            radial_quadrature,
            solve_eigenpair,
        )
        native = result._obj
        native.has_eigenpair = solve_eigenpair
        native.field_periods = 1
        native.toroidal_mode = toroidal_mode
        native.poloidal_max = poloidal_max
        native.mode_count = 2 * poloidal_max + 1
        native.radial_surfaces = 768
        native.parity_class = 1
        native.radial_quadrature = radial_quadrature
        native.negative_count = 1
        native.lowest_eigenvalue = -2.5e-7 if solve_eigenpair else np.nan
        native.certificate = 1.0e-11 if solve_eigenpair else np.nan
        native.eigenpair_residual = 3.0e-14 if solve_eigenpair else np.nan
        native.force_balance_residual = 4.0e-10
        error.value = b""
        return 0


def equilibrium(library):
    result = gliss.Equilibrium.__new__(gliss.Equilibrium)
    result.path = Path("solovev.nc")
    result._library = library
    result._handle = ctypes.c_void_p(17)
    return result


def test_axisymmetric_inertia_uses_loaded_equilibrium():
    library = FakeLibrary()

    result = gliss.axisymmetric_inertia(
        equilibrium(library),
        toroidal_mode=2,
        poloidal_max=7,
        radial_quadrature="gauss2",
    )

    assert library.arguments == (17, 2, 7, 2, 0)
    assert result.has_eigenpair is False
    assert result.field_periods == 1
    assert result.toroidal_mode == 2
    assert result.poloidal_max == 7
    assert result.mode_count == 15
    assert result.radial_surfaces == 768
    assert result.parity_class == 1
    assert result.radial_quadrature == "gauss2"
    assert result.negative_count == 1
    assert result.lowest_eigenvalue is None
    assert result.certificate is None
    assert result.eigenpair_residual is None
    assert result.force_balance_residual == pytest.approx(4.0e-10)


def test_solve_axisymmetric_returns_certified_pair():
    library = FakeLibrary()

    result = gliss.solve_axisymmetric(equilibrium(library))

    assert library.arguments == (17, 1, 8, 1, 1)
    assert result.has_eigenpair is True
    assert result.lowest_eigenvalue == pytest.approx(-2.5e-7)
    assert result.certificate == pytest.approx(1.0e-11)
    assert result.eigenpair_residual == pytest.approx(3.0e-14)


@pytest.mark.parametrize("function", [gliss.axisymmetric_inertia, gliss.solve_axisymmetric])
def test_axisymmetric_rejects_closed_equilibrium(function):
    value = equilibrium(FakeLibrary())
    value._handle = ctypes.c_void_p()

    with pytest.raises(RuntimeError, match="Equilibrium is closed"):
        function(value)


def test_axisymmetric_rejects_non_equilibrium():
    with pytest.raises(TypeError, match="gliss.Equilibrium"):
        gliss.axisymmetric_inertia(object())


@pytest.mark.parametrize(
    ("keyword", "value", "exception", "match"),
    [
        ("toroidal_mode", True, TypeError, "must be an integer"),
        ("toroidal_mode", 0, ValueError, "must be positive"),
        ("poloidal_max", 0, ValueError, "must be positive"),
        ("radial_quadrature", 1, TypeError, "must be a string"),
        ("radial_quadrature", "bad", ValueError, "midpoint or gauss2"),
    ],
)
def test_axisymmetric_rejects_invalid_input(keyword, value, exception, match):
    with pytest.raises(exception, match=match):
        gliss.axisymmetric_inertia(equilibrium(FakeLibrary()), **{keyword: value})


def test_axisymmetric_rejects_invalid_native_result():
    library = FakeLibrary()

    def invalid(*arguments):
        result = arguments[5]._obj
        result.has_eigenpair = 0
        result.field_periods = 2
        return 0

    library.gliss_axisymmetric_spectrum = FakeFunction(invalid)

    with pytest.raises(gliss.GlissInternalError, match="invalid axisymmetric result"):
        gliss.axisymmetric_inertia(equilibrium(library))


def test_axisymmetric_propagates_native_error():
    library = FakeLibrary()

    def fail(*arguments):
        arguments[6].value = b"equilibrium contains nonaxisymmetric harmonics"
        return 4

    library.gliss_axisymmetric_spectrum = FakeFunction(fail)

    with pytest.raises(gliss.GlissArgumentError, match="nonaxisymmetric harmonics"):
        gliss.axisymmetric_inertia(equilibrium(library))


def test_axisymmetric_binding_requires_symbol():
    with pytest.raises(OSError, match="gliss_axisymmetric_spectrum"):
        _bind(object())
