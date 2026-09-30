"""Tests for gliss.mercier: ABI contract and golden-fixture regression."""

import csv
import os
from pathlib import Path

import numpy as np
import pytest

import gliss


def test_mercier_profile_nonexistent_path_raises():
    with pytest.raises(FileNotFoundError, match="does not exist"):
        gliss.mercier_profile("/nonexistent/path.nc")


@pytest.mark.parametrize("resolution", [True, 0, -1, 1.5, "64"])
def test_mercier_profile_rejects_invalid_resolution(tmp_path, resolution):
    export = tmp_path / "equilibrium.nc"
    export.touch()
    with pytest.raises((TypeError, ValueError), match="n_theta"):
        gliss.mercier_profile(export, n_theta=resolution)


@pytest.mark.native
def test_mercier_profile_accepts_pathlike(native_library, tmp_path):
    export = Path(tmp_path, "invalid.nc")
    export.touch()
    with pytest.raises(RuntimeError, match="failed to read"):
        gliss.mercier_profile(export)


def test_mercier_profile_rejects_embedded_nul(tmp_path):
    with pytest.raises(ValueError, match="null byte"):
        gliss.mercier_profile(os.fspath(tmp_path) + "/bad\0name.nc")


def _load_golden(path):
    with open(path, newline="") as handle:
        rows = list(csv.reader(handle))
    return np.array([[float(value) for value in row] for row in rows[1:]])


def test_mercier_profile_matches_golden(request):
    fixture = os.environ.get("GLISS_MERCIER_FIXTURE")
    golden_path = os.environ.get("GLISS_MERCIER_GOLDEN")
    if not fixture or not golden_path:
        pytest.skip("GLISS_MERCIER_FIXTURE and GLISS_MERCIER_GOLDEN not set")
    request.getfixturevalue("native_library")

    golden = _load_golden(golden_path)
    s, d_mercier = gliss.mercier_profile(fixture)

    np.testing.assert_array_equal(s, golden[:, 0])
    np.testing.assert_allclose(d_mercier, golden[:, 5], rtol=1e-9)


def test_mercier_objective_is_instability_of_least_stable_surface(monkeypatch):
    profile = (np.array([0.25, 0.5, 0.75]), np.array([0.3, -0.2, 0.1]))
    monkeypatch.setattr(gliss.mercier, "mercier_profile", lambda *a, **k: profile)
    assert gliss.mercier.mercier_objective("unused.nc") == pytest.approx(0.2)
    stable = (profile[0], np.array([0.3, 0.2, 0.1]))
    monkeypatch.setattr(gliss.mercier, "mercier_profile", lambda *a, **k: stable)
    assert gliss.mercier.mercier_objective("unused.nc") == pytest.approx(-0.1)


def test_mercier_objective_rejects_nan_profile(monkeypatch):
    # A NaN surface must not reach an optimizer as a NaN objective.
    profile = (np.array([0.25, 0.5]), np.array([0.3, np.nan]))
    monkeypatch.setattr(gliss.mercier, "mercier_profile", lambda *a, **k: profile)
    with pytest.raises(RuntimeError, match="non-finite"):
        gliss.mercier.mercier_objective("unused.nc")
