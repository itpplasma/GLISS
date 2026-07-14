import ctypes
import os

import pytest

import gliss
from gliss.terpsichore import _bind


class FakeFunction:
    def __init__(self, callback):
        self.callback = callback
        self.argtypes = None
        self.restype = None

    def __call__(self, *args):
        return self.callback(*args)


class FakeLibrary:
    def __init__(self):
        self.path = None
        self.gliss_terpsichore_fixed_boundary = FakeFunction(self.solve)

    def solve(self, path, length, result, error, error_capacity):
        self.path = ctypes.string_at(path, length)
        result._obj.unknowns = 24257
        result._obj.negative_count = 5
        result._obj.eigenvalue = -7.3862621032169963e-7
        result._obj.certificate = 1.9531250000095766e-11
        result._obj.residual = 5.1278099516263922e-15
        result._obj.resolution = 7.5810235738805062e-10
        error.value = b""
        return 0


def test_fixed_boundary_result(monkeypatch, tmp_path):
    library = FakeLibrary()
    monkeypatch.setattr("gliss.terpsichore._load_library", lambda: library)
    fixture = tmp_path / "fort.23"
    fixture.touch()

    result = gliss.solve_terpsichore_fixed_boundary(fixture)

    assert library.path == os.fsencode(fixture)
    assert result.unknowns == 24257
    assert result.negative_count == 5
    assert result.eigenvalue == pytest.approx(-7.3862621032169963e-7)
    assert result.certificate == pytest.approx(1.9531250000095766e-11)
    assert result.residual == pytest.approx(5.1278099516263922e-15)
    assert result.resolution == pytest.approx(7.5810235738805062e-10)


def test_fixed_boundary_native_error(monkeypatch, tmp_path):
    library = FakeLibrary()

    def fail(path, length, result, error, error_capacity):
        error.value = b"invalid TERPSICHORE FORT.23"
        return 1

    library.gliss_terpsichore_fixed_boundary = FakeFunction(fail)
    monkeypatch.setattr("gliss.terpsichore._load_library", lambda: library)
    fixture = tmp_path / "fort.23"
    fixture.touch()

    with pytest.raises(gliss.GlissIOError, match="invalid TERPSICHORE FORT.23"):
        gliss.solve_terpsichore_fixed_boundary(fixture)


def test_fixed_boundary_rejects_invalid_native_result(monkeypatch, tmp_path):
    library = FakeLibrary()

    def return_invalid(path, length, result, error, error_capacity):
        result._obj.unknowns = 1
        result._obj.negative_count = 0
        result._obj.eigenvalue = 1.0
        return 0

    library.gliss_terpsichore_fixed_boundary = FakeFunction(return_invalid)
    monkeypatch.setattr("gliss.terpsichore._load_library", lambda: library)
    fixture = tmp_path / "fort.23"
    fixture.touch()

    with pytest.raises(gliss.GlissInternalError, match="invalid TERPSICHORE"):
        gliss.solve_terpsichore_fixed_boundary(fixture)


@pytest.mark.parametrize(
    ("path", "exception", "match"),
    [
        (b"fort.23", TypeError, "resolve to a string"),
        ("bad\0fort.23", ValueError, "null byte"),
        ("missing", FileNotFoundError, "does not exist"),
    ],
)
def test_fixed_boundary_rejects_invalid_paths(path, exception, match):
    with pytest.raises(exception, match=match):
        gliss.solve_terpsichore_fixed_boundary(path)


def test_fixed_boundary_binding_requires_symbol():
    with pytest.raises(OSError, match="gliss_terpsichore_fixed_boundary"):
        _bind(object())
