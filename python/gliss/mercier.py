"""Mercier stability diagnostics for GVEC/CAS3D equilibrium exports."""

import numpy as np

from .equilibrium import Equilibrium, PathLike, _resolution


def mercier_profile(
    path: PathLike, n_theta: int = 64, n_zeta: int = 64
) -> tuple[np.ndarray, np.ndarray]:
    """Return radial coordinates and Mercier discriminants.

    Parameters
    ----------
    path : str or os.PathLike
        GVEC/CAS3D NetCDF equilibrium export.
    n_theta, n_zeta : int, optional
        Positive angular quadrature sizes. Both default to 64.

    Returns
    -------
    s, d_mercier : numpy.ndarray
        One-dimensional ``float64`` arrays. Positive ``d_mercier`` is
        Mercier-unstable.
    """
    n_theta = _resolution(n_theta, "n_theta")
    n_zeta = _resolution(n_zeta, "n_zeta")
    with Equilibrium(path) as equilibrium:
        return equilibrium.mercier_profile(n_theta=n_theta, n_zeta=n_zeta)


def mercier_objective(path: PathLike, n_theta: int = 64, n_zeta: int = 64) -> float:
    """Return the most unstable Mercier discriminant in the profile."""
    _, d_mercier = mercier_profile(path, n_theta=n_theta, n_zeta=n_zeta)
    if d_mercier.size == 0:
        raise RuntimeError("GLISS returned an empty Mercier profile")
    return float(np.max(d_mercier))
