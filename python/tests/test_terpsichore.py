import ctypes
import os

import pytest

import gliss
from gliss.terpsichore import _bind, _bind_pseudoplasma


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
        self.gliss_terpsichore_pseudoplasma = FakeFunction(self.solve_pseudoplasma)

    def solve(self, path, length, result, error, error_capacity):
        self.path = ctypes.string_at(path, length)
        result._obj.unknowns = 24257
        result._obj.negative_count = 5
        result._obj.eigenvalue = -7.3862621032169963e-7
        result._obj.certificate = 1.9531250000095766e-11
        result._obj.residual = 5.1278099516263922e-15
        result._obj.resolution = 7.5810235738805062e-10
        result._obj.reference_eigenvalue = -7.38626214e-7
        result._obj.reference_potential = -3.01081e-5
        result._obj.computed_potential = -3.01081001e-5
        result._obj.reference_kinetic = 40.8167
        result._obj.computed_kinetic = 40.8167
        result._obj.reference_residual = 3.8e-5
        result._obj.mode_overlap = 0.99999999
        error.value = b""
        return 0

    def solve_pseudoplasma(
        self,
        matrix_path,
        matrix_length,
        vacuum_intervals,
        vacuum_path,
        vacuum_length,
        result,
        error,
        error_capacity,
    ):
        self.matrix_path = ctypes.string_at(matrix_path, matrix_length)
        self.vacuum_path = ctypes.string_at(vacuum_path, vacuum_length)
        self.vacuum_intervals = vacuum_intervals
        result._obj.unknowns = 24448
        result._obj.negative_count = 7
        result._obj.eigenvalue = -4.870618368605979e-5
        result._obj.certificate = 1.2499999999993528e-9
        result._obj.residual = 1.1746531588612227e-13
        result._obj.resolution = 4.235602369785144e-9
        result._obj.growth_rate = 8.315076745275407e-3
        result._obj.reference_eigenvalue = -4.87061837182045e-5
        result._obj.reference_potential = -9.949905417601407e-5
        result._obj.computed_potential = -9.949905415935845e-5
        result._obj.reference_kinetic = 2.0428423370567863
        result._obj.computed_kinetic = 2.0428423370567965
        result._obj.reference_residual = 3.7920576355991595e-5
        result._obj.mode_overlap = 0.9999999992664245
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
    assert result.reference_eigenvalue == pytest.approx(-7.38626214e-7)
    assert result.reference_potential == pytest.approx(-3.01081e-5)
    assert result.computed_potential == pytest.approx(-3.01081001e-5)
    assert result.reference_kinetic == pytest.approx(40.8167)
    assert result.computed_kinetic == pytest.approx(40.8167)
    assert result.reference_residual == pytest.approx(3.8e-5)
    assert result.mode_overlap == pytest.approx(0.99999999)


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


def test_pseudoplasma_result(monkeypatch, tmp_path):
    library = FakeLibrary()
    monkeypatch.setattr("gliss.terpsichore._load_library", lambda: library)
    matrix = tmp_path / "fort.23"
    vacuum = tmp_path / "fort.24"
    matrix.touch()
    vacuum.touch()

    result = gliss.solve_terpsichore_pseudoplasma(matrix, 16, vacuum)

    assert library.matrix_path == os.fsencode(matrix)
    assert library.vacuum_path == os.fsencode(vacuum)
    assert library.vacuum_intervals == 16
    assert result.unknowns == 24448
    assert result.negative_count == 7
    assert result.eigenvalue == pytest.approx(-4.870618368605979e-5)
    assert result.growth_rate == pytest.approx(8.315076745275407e-3)
    assert result.reference_eigenvalue == pytest.approx(-4.87061837182045e-5)
    assert result.reference_potential == pytest.approx(-9.949905417601407e-5)
    assert result.computed_potential == pytest.approx(-9.949905415935845e-5)
    assert result.reference_kinetic == pytest.approx(2.0428423370567863)
    assert result.computed_kinetic == pytest.approx(2.0428423370567965)
    assert result.reference_residual == pytest.approx(3.7920576355991595e-5)
    assert result.mode_overlap == pytest.approx(0.9999999992664245)
    assert gliss.TerpsichorePseudoplasmaResult is type(result)


@pytest.mark.parametrize(
    ("intervals", "exception", "match"),
    [
        (True, TypeError, "must be an integer"),
        (0, ValueError, "must be positive"),
        (2**31, ValueError, "signed 32-bit"),
    ],
)
def test_pseudoplasma_rejects_invalid_intervals(
    monkeypatch, tmp_path, intervals, exception, match
):
    monkeypatch.setattr("gliss.terpsichore._load_library", lambda: FakeLibrary())
    matrix = tmp_path / "fort.23"
    vacuum = tmp_path / "fort.24"
    matrix.touch()
    vacuum.touch()

    with pytest.raises(exception, match=match):
        gliss.solve_terpsichore_pseudoplasma(matrix, intervals, vacuum)


def test_pseudoplasma_rejects_missing_vacuum(tmp_path):
    matrix = tmp_path / "fort.23"
    matrix.touch()

    with pytest.raises(FileNotFoundError, match="FORT.24 does not exist"):
        gliss.solve_terpsichore_pseudoplasma(matrix, 16, tmp_path / "fort.24")


def test_pseudoplasma_validates_interval_before_paths():
    with pytest.raises(ValueError, match="must be positive"):
        gliss.solve_terpsichore_pseudoplasma("missing.23", 0, "missing.24")


def test_pseudoplasma_native_error(monkeypatch, tmp_path):
    library = FakeLibrary()

    def fail(*arguments):
        error = arguments[-2]
        error.value = b"TERPSICHORE FORT.24 IVAC does not match the request"
        return 1

    library.gliss_terpsichore_pseudoplasma = FakeFunction(fail)
    monkeypatch.setattr("gliss.terpsichore._load_library", lambda: library)
    matrix = tmp_path / "fort.23"
    vacuum = tmp_path / "fort.24"
    matrix.touch()
    vacuum.touch()

    with pytest.raises(gliss.GlissIOError, match="IVAC does not match"):
        gliss.solve_terpsichore_pseudoplasma(matrix, 16, vacuum)


def test_pseudoplasma_rejects_invalid_native_result(monkeypatch, tmp_path):
    library = FakeLibrary()

    def return_invalid(*arguments):
        result = arguments[-3]
        result._obj.unknowns = 1
        result._obj.negative_count = 1
        result._obj.eigenvalue = -1.0
        result._obj.certificate = 0.0
        result._obj.residual = 0.0
        result._obj.resolution = 0.0
        result._obj.growth_rate = 1.0
        result._obj.reference_kinetic = 1.0
        result._obj.computed_kinetic = 1.0
        result._obj.reference_residual = 0.0
        result._obj.mode_overlap = 1.01
        return 0

    library.gliss_terpsichore_pseudoplasma = FakeFunction(return_invalid)
    monkeypatch.setattr("gliss.terpsichore._load_library", lambda: library)
    matrix = tmp_path / "fort.23"
    vacuum = tmp_path / "fort.24"
    matrix.touch()
    vacuum.touch()

    with pytest.raises(gliss.GlissInternalError, match="invalid TERPSICHORE"):
        gliss.solve_terpsichore_pseudoplasma(matrix, 16, vacuum)


def test_pseudoplasma_binding_requires_symbol():
    with pytest.raises(OSError, match="gliss_terpsichore_pseudoplasma"):
        _bind_pseudoplasma(object())
