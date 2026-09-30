"""Tests for gliss.mercier: ABI contract and golden-fixture regression."""

import csv
import os
from pathlib import Path
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
    with pytest.raises(FileNotFoundError, match="does not exist"):
        gliss.mercier_profile("/nonexistent/path.nc")


@pytest.mark.parametrize("resolution", [True, 0, -1, 1.5, "64"])
def test_mercier_profile_rejects_invalid_resolution(tmp_path, resolution):
    export = tmp_path / "equilibrium.nc"
    export.touch()
    with pytest.raises((TypeError, ValueError), match="n_theta"):
        gliss.mercier_profile(export, n_theta=resolution)


def test_mercier_profile_accepts_pathlike(tmp_path):
    _ensure_library()
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


_SOLOVEV = Path(__file__).resolve().parents[2] / "test" / "data" / "solovev_q1.035.nc"


def test_mercier_objective_is_instability_of_least_stable_surface(monkeypatch):
    profile = (np.array([0.25, 0.5, 0.75]), np.array([0.3, -0.2, 0.1]))
    monkeypatch.setattr(gliss.mercier, "mercier_profile", lambda *a, **k: profile)
    assert gliss.mercier.mercier_objective("unused.nc") == pytest.approx(0.2)
    stable = (profile[0], np.array([0.3, 0.2, 0.1]))
    monkeypatch.setattr(gliss.mercier, "mercier_profile", lambda *a, **k: stable)
    assert gliss.mercier.mercier_objective("unused.nc") == pytest.approx(-0.1)


def test_solovev_mercier_sign_matches_dcon():
    # Independent GPEC/DCON Newcomb runs of this analytic Solov'ev family near
    # q0=1.04 find the core Mercier-unstable (D_I > 0) for s below about 0.16
    # and stable outside.  Positive D_Mercier is stable in GLISS.
    _ensure_library()
    s, d_mercier = gliss.mercier_profile(_SOLOVEV)
    assert d_mercier[0] < 0.0
    assert np.all(d_mercier[s > 0.25] > 0.0)
    crossing = s[np.argmax(d_mercier > 0.0)]
    assert 0.1 < crossing < 0.3
    assert gliss.mercier_objective(_SOLOVEV) == pytest.approx(-d_mercier.min())
    assert gliss.mercier_objective(_SOLOVEV) > 0.0


def test_gvec_export_reports_left_handed_chart():
    _ensure_library()
    with gliss.Equilibrium(_SOLOVEV) as equilibrium:
        assert equilibrium.coordinate_handedness == "left-handed"
