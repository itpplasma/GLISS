import ctypes

import numpy as np
import pytest

import gliss
from gliss.equilibrium import Equilibrium, GlissIOError, _bind


class FakeFunction:
    def __init__(self, callback):
        self.callback = callback
        self.argtypes = None
        self.restype = None

    def __call__(self, *args):
        return self.callback(*args)


class FakeLibrary:
    def __init__(self):
        self.creates = 0
        self.destroys = 0
        self.gliss_equilibrium_create = FakeFunction(self.create)
        self.gliss_equilibrium_destroy = FakeFunction(self.destroy)
        self.gliss_equilibrium_surface_count = FakeFunction(self.surface_count)
        self.gliss_equilibrium_schema_version = FakeFunction(self.schema_version)
        self.gliss_equilibrium_write = FakeFunction(self.write)
        self.gliss_mercier_profile_context = FakeFunction(self.mercier_profile)

    def create(self, path, length, handle, error, error_capacity):
        self.creates += 1
        handle._obj.value = self.creates
        error.value = b""
        return 0

    def destroy(self, handle, error, error_capacity):
        self.destroys += 1
        handle._obj.value = None
        if error is not None:
            error.value = b""
        return 0

    def surface_count(self, handle, surfaces, error, error_capacity):
        surfaces._obj.value = 3
        error.value = b""
        return 0

    def schema_version(self, handle, version, error, error_capacity):
        version._obj.value = 0
        error.value = b""
        return 0

    def write(self, handle, path, length, error, error_capacity):
        output = ctypes.string_at(path, length).decode()
        with open(output, "wb") as stream:
            stream.write(b"schema-1-export")
        error.value = b""
        return 0

    def mercier_profile(
        self,
        handle,
        n_theta,
        n_zeta,
        capacity,
        s_values,
        d_mercier,
        written,
        error,
        error_capacity,
    ):
        np.ctypeslib.as_array(s_values, shape=(capacity,))[:] = [0.1, 0.5, 0.9]
        np.ctypeslib.as_array(d_mercier, shape=(capacity,))[:] = [-1.0, 0.0, 2.0]
        written._obj.value = 3
        error.value = b""
        return 0


def test_equilibrium_context_lifecycle(monkeypatch, tmp_path):
    library = FakeLibrary()
    monkeypatch.setattr("gliss.equilibrium._load_library", lambda: library)
    export = tmp_path / "equilibrium.nc"
    export.touch()

    equilibrium = Equilibrium(export)
    assert repr(equilibrium) == (
        "<gliss.Equilibrium(path='equilibrium.nc', state='open')>"
    )
    assert str(tmp_path) not in repr(equilibrium)
    assert not equilibrium.closed
    assert equilibrium.schema_version == 0
    s_values, d_mercier = equilibrium.mercier_profile(n_theta=32, n_zeta=16)
    np.testing.assert_array_equal(s_values, [0.1, 0.5, 0.9])
    np.testing.assert_array_equal(d_mercier, [-1.0, 0.0, 2.0])
    equilibrium.close()
    equilibrium.close()
    assert equilibrium.closed
    assert library.creates == 1
    assert library.destroys == 1
    with pytest.raises(RuntimeError, match="closed"):
        equilibrium.mercier_profile()
    with pytest.raises(RuntimeError, match="closed"):
        _ = equilibrium.schema_version
    with pytest.raises(RuntimeError, match="closed"):
        equilibrium.write(tmp_path / "closed.nc")


def test_equilibrium_write_atomically_replaces_destination(monkeypatch, tmp_path):
    library = FakeLibrary()
    monkeypatch.setattr("gliss.equilibrium._load_library", lambda: library)
    export = tmp_path / "equilibrium.nc"
    export.touch()
    destination = tmp_path / "copy.nc"
    destination.write_bytes(b"old")

    with Equilibrium(export) as equilibrium:
        result = equilibrium.write(destination)

    assert result == destination
    assert destination.read_bytes() == b"schema-1-export"
    assert sorted(tmp_path.iterdir()) == sorted((export, destination))


def test_equilibrium_write_preserves_destination_on_native_error(monkeypatch, tmp_path):
    library = FakeLibrary()

    def fail_write(handle, path, length, error, error_capacity):
        error.value = b"write failed"
        return 1

    library.gliss_equilibrium_write = FakeFunction(fail_write)
    monkeypatch.setattr("gliss.equilibrium._load_library", lambda: library)
    export = tmp_path / "equilibrium.nc"
    export.touch()
    destination = tmp_path / "copy.nc"
    destination.write_bytes(b"old")

    with Equilibrium(export) as equilibrium:
        with pytest.raises(GlissIOError, match="write failed"):
            equilibrium.write(destination)

    assert destination.read_bytes() == b"old"
    assert sorted(tmp_path.iterdir()) == sorted((export, destination))


@pytest.mark.parametrize(
    ("output", "exception", "match"),
    [
        (b"bytes.nc", TypeError, "resolve to a string"),
        ("nul\0name.nc", ValueError, "null byte"),
    ],
)
def test_equilibrium_write_rejects_invalid_paths(
    monkeypatch, tmp_path, output, exception, match
):
    library = FakeLibrary()
    monkeypatch.setattr("gliss.equilibrium._load_library", lambda: library)
    export = tmp_path / "equilibrium.nc"
    export.touch()
    with Equilibrium(export) as equilibrium:
        with pytest.raises(exception, match=match):
            equilibrium.write(output)


def test_equilibrium_write_rejects_missing_directory_and_directory(
    monkeypatch, tmp_path
):
    library = FakeLibrary()
    monkeypatch.setattr("gliss.equilibrium._load_library", lambda: library)
    export = tmp_path / "equilibrium.nc"
    export.touch()
    with Equilibrium(export) as equilibrium:
        with pytest.raises(FileNotFoundError, match="output directory"):
            equilibrium.write(tmp_path / "missing" / "copy.nc")
        with pytest.raises(ValueError, match="output path is a directory"):
            equilibrium.write(tmp_path)


def test_equilibrium_context_manager_allows_independent_instances(
    monkeypatch, tmp_path
):
    library = FakeLibrary()
    monkeypatch.setattr("gliss.equilibrium._load_library", lambda: library)
    first_path = tmp_path / "first.nc"
    second_path = tmp_path / "second.nc"
    first_path.touch()
    second_path.touch()

    with Equilibrium(first_path) as first, Equilibrium(second_path) as second:
        assert first._handle.value != second._handle.value
        assert not first.closed
        assert not second.closed
    assert library.destroys == 2


def test_equilibrium_create_reports_native_error(monkeypatch, tmp_path):
    library = FakeLibrary()

    def fail_create(path, length, handle, error, error_capacity):
        error.value = b"failed to read equilibrium export"
        return 1

    library.gliss_equilibrium_create = FakeFunction(fail_create)
    monkeypatch.setattr("gliss.equilibrium._load_library", lambda: library)
    export = tmp_path / "invalid.nc"
    export.touch()

    with pytest.raises(GlissIOError, match="failed to read equilibrium export"):
        Equilibrium(export)


def test_equilibrium_failed_partial_create_is_cleaned_up(monkeypatch, tmp_path):
    library = FakeLibrary()

    def fail_create(path, length, handle, error, error_capacity):
        handle._obj.value = 1234
        error.value = b"failed after allocation"
        return 1

    library.gliss_equilibrium_create = FakeFunction(fail_create)
    monkeypatch.setattr("gliss.equilibrium._load_library", lambda: library)
    export = tmp_path / "invalid.nc"
    export.touch()

    with pytest.raises(GlissIOError, match="failed after allocation"):
        Equilibrium(export)
    assert library.destroys == 1


def test_equilibrium_rejects_source_changed_during_load(monkeypatch, tmp_path):
    library = FakeLibrary()

    def changing_create(path, length, handle, error, error_capacity):
        export = ctypes.string_at(path, length).decode()
        with open(export, "wb") as stream:
            stream.write(b"changed")
        handle._obj.value = 1234
        error.value = b""
        return 0

    library.gliss_equilibrium_create = FakeFunction(changing_create)
    monkeypatch.setattr("gliss.equilibrium._load_library", lambda: library)
    export = tmp_path / "equilibrium.nc"
    export.touch()

    with pytest.raises(GlissIOError, match="changed while loading"):
        Equilibrium(export)
    assert library.destroys == 1


def test_equilibrium_is_public():
    assert gliss.Equilibrium is Equilibrium


def test_equilibrium_bind_reports_missing_native_symbols():
    with pytest.raises(OSError, match="equilibrium context.*matching"):
        _bind(object())
