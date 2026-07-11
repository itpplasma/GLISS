"""SIMSOPT integration for the GLISS Mercier stability objective.

``MercierPenalty`` wraps :func:`gliss.mercier.mercier_objective` as a
``simsopt`` ``Optimizable`` leaf node, so GLISS's validated worst-case
D_Mercier value can be used as a term in a SIMSOPT stellarator
optimization problem (e.g. composed into a
``simsopt.objectives.LeastSquaresProblem``).

Scope: this wrapper exposes the Mercier objective *value* computed from
a CAS3D/GVEC equilibrium export. It does not carry GVEC equilibrium
degrees of freedom and does not provide analytic derivatives -- both
come from the upstream GVEC equilibrium solve (the differentiable
chain that produces the export this class reads), which is a separate
follow-up. Between evaluations, an external GVEC solve step is expected
to write a new export and update ``export_path`` so the next ``J()``
call picks it up.
"""
from simsopt._core.optimizable import Optimizable

from .mercier import mercier_objective


class MercierPenalty(Optimizable):
    """Leaf SIMSOPT objective returning ``mercier_objective(export_path)``.

    Has no free DOFs: it is a value objective over an externally
    produced CAS3D/GVEC export, not a parameterized equilibrium model.
    """

    def __init__(self, export_path):
        self._export_path = str(export_path)
        super().__init__(x0=[])

    @property
    def export_path(self):
        """Path to the CAS3D/GVEC export read on the next ``J()`` call."""
        return self._export_path

    @export_path.setter
    def export_path(self, path):
        self._export_path = str(path)
        self.set_recompute_flag()

    def J(self):
        """Worst-case ``max(D_Mercier)`` for the current export path."""
        return mercier_objective(self._export_path)

    return_fn_map = {'J': J}


def mercier_penalty(export_path):
    """Return a :class:`MercierPenalty` wrapping ``export_path``."""
    return MercierPenalty(export_path)
