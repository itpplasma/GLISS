"""Shared fixtures that separate native evidence from contract tests.

Tests marked ``native`` load the real libgliss_c and check it against
independent oracles. They fail, rather than skip, when the library or the
repository fixtures are missing, so a native run cannot pass vacuously.
Contract tests substitute fakes for the library and check only the Python
side of the ABI.
"""

import os
from pathlib import Path

import pytest

_REPOSITORY = Path(__file__).resolve().parents[2]


def _build_library():
    # The CMake build tree of a checkout; the suffix is platform specific.
    for suffix in (".so", ".dylib"):
        candidate = _REPOSITORY / "build" / f"libgliss_c{suffix}"
        if candidate.is_file():
            return candidate
    return None


def pytest_configure(config):
    config.addinivalue_line(
        "markers", "native: loads libgliss_c and checks it against an oracle"
    )


@pytest.fixture(scope="session")
def native_library():
    library = _build_library()
    if not os.environ.get("GLISS_LIB") and library is not None:
        os.environ["GLISS_LIB"] = str(library)
    import gliss

    try:
        library = gliss._load_library()
    except OSError as error:
        pytest.fail(f"native tests require libgliss_c: {error}")
    assert gliss.version() == gliss.__version__
    return library


@pytest.fixture(scope="session")
def test_data():
    directory = Path(os.environ.get("GLISS_TEST_DATA", _REPOSITORY / "test" / "data"))
    if not (directory / "solovev_q1.035.nc").is_file():
        pytest.fail(f"native tests require the fixtures in {directory}")
    return directory
