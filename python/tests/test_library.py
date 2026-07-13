import ctypes
from pathlib import Path

import pytest

import gliss


def test_python_and_compiled_versions_match():
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
            return 2

    class Library:
        gliss_abi_version = Function()

    monkeypatch.setattr(gliss, "_open_library", lambda: Library())
    with pytest.raises(OSError, match="ABI version 2.*requires 1"):
        gliss._load_library()
