import ctypes
from pathlib import Path

import pytest

import gliss


@pytest.mark.native
def test_python_and_compiled_versions_match(native_library):
    assert gliss.version() == gliss.__version__


def test_bundled_library_is_preferred(monkeypatch, tmp_path):
    bundled = tmp_path / "libgliss_c.so"
    bundled.touch()
    monkeypatch.delenv("GLISS_LIB", raising=False)
    monkeypatch.setattr(gliss, "_bundled_library", lambda: bundled)
    monkeypatch.setattr(ctypes, "CDLL", lambda path: Path(path))
    assert gliss._open_library() == bundled


def test_incompatible_abi_is_rejected(monkeypatch):
    class Function:
        def __call__(self):
            return 1

    class Library:
        gliss_abi_version = Function()

    monkeypatch.setattr(gliss, "_open_library", lambda: Library())
    with pytest.raises(OSError, match="ABI version 1.*requires 4"):
        gliss._load_library()


@pytest.mark.parametrize("revision", [3, 5])
def test_incompatible_operator_revision_is_rejected(monkeypatch, revision):
    calls = []

    class Function:
        def __init__(self, name, value):
            self.name, self.value = name, value

        def __call__(self):
            calls.append(self.name)
            return self.value

    class Library:
        gliss_abi_version = Function("abi", 4)
        gliss_discretization_revision = Function("revision", revision)

    monkeypatch.setattr(gliss, "_open_library", lambda: Library())
    with pytest.raises(OSError, match=f"discretization revision {revision}.*requires 4"):
        gliss._load_library()
    assert calls == ["abi", "revision"]


def test_get_include_returns_bundled_header_directory(monkeypatch, tmp_path):
    include = tmp_path / "include"
    include.mkdir()
    (include / "gliss.h").touch()
    monkeypatch.setattr(gliss, "files", lambda package: tmp_path)
    assert Path(gliss.get_include()) == include
