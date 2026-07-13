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
    assert repr(equilibrium) == (f"<gliss.Equilibrium(path={export!r}, state='open')>")
    assert not equilibrium.closed
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


def test_equilibrium_is_public():
    assert gliss.Equilibrium is Equilibrium


def test_equilibrium_bind_reports_missing_native_symbols():
    with pytest.raises(OSError, match="equilibrium context.*matching"):
        _bind(object())
