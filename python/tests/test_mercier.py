"""Tests for gliss.mercier: ABI contract and golden-fixture regression."""
import csv
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))

import gliss  # noqa: E402

_DEFAULT_LIB = os.path.join(
    os.path.dirname(__file__), os.pardir, os.pardir, "build", "libgliss_c.so"
)


def _ensure_library():
    if os.environ.get("GLISS_LIB"):
        return
    if os.path.exists(_DEFAULT_LIB):
        os.environ["GLISS_LIB"] = _DEFAULT_LIB
    else:
        pytest.skip("libgliss_c.so not found; set GLISS_LIB to the built library")


def test_mercier_profile_nonexistent_path_raises():
    _ensure_library()
    with pytest.raises(RuntimeError):
        gliss.mercier_profile("/nonexistent/path.nc")


def _load_golden(path):
    with open(path, newline="") as handle:
        rows = list(csv.reader(handle))
    return np.array([[float(value) for value in row] for row in rows[1:]])


def test_mercier_profile_matches_golden():
    fixture = os.environ.get("GLISS_MERCIER_FIXTURE")
    golden_path = os.environ.get("GLISS_MERCIER_GOLDEN")
    if not fixture or not golden_path:
        pytest.skip("GLISS_MERCIER_FIXTURE and GLISS_MERCIER_GOLDEN not set")
    _ensure_library()

    golden = _load_golden(golden_path)
    s, d_mercier = gliss.mercier_profile(fixture)

    np.testing.assert_array_equal(s, golden[:, 0])
    np.testing.assert_allclose(d_mercier, golden[:, 5], rtol=1e-9)
