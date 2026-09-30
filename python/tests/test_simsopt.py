"""Tests for gliss.simsopt: the SIMSOPT Optimizable wrapper for mercier_objective."""

import math
import os

import pytest

pytest.importorskip("simsopt")

from simsopt.objectives.least_squares import LeastSquaresProblem  # noqa: E402

from gliss.simsopt import MercierPenalty  # noqa: E402

W7X_EXPORT = os.environ.get("GLISS_MERCIER_W7X")
STABLE_EXPORT = os.environ.get("GLISS_MERCIER_STABLE")


@pytest.mark.skipif(not W7X_EXPORT, reason="GLISS_MERCIER_W7X not set")
def test_unstable_export_gives_large_positive_j(native_library):
    penalty = MercierPenalty(W7X_EXPORT)
    assert penalty.J() > 0.5


@pytest.mark.skipif(not STABLE_EXPORT, reason="GLISS_MERCIER_STABLE not set")
def test_stable_export_gives_near_zero_j(native_library):
    penalty = MercierPenalty(STABLE_EXPORT)
    assert penalty.J() < 1e-6


@pytest.mark.skipif(not W7X_EXPORT, reason="GLISS_MERCIER_W7X not set")
def test_composes_into_least_squares_problem(native_library):
    penalty = MercierPenalty(W7X_EXPORT)
    prob = LeastSquaresProblem.from_tuples([(penalty.J, 0.0, 1.0)])
    j = penalty.J()
    objective = prob.objective()
    assert math.isfinite(objective)
    assert objective == pytest.approx(j * j)
