"""Mercier stability diagnostics for GVEC/CAS3D equilibrium exports."""

import ctypes
import operator
import os
from pathlib import Path

import numpy as np

from . import _load_library

_C_INT_MAX = 2 ** (ctypes.sizeof(ctypes.c_int) * 8 - 1) - 1
_STATUS_MESSAGES = {
    1: "failed to read equilibrium export file",
    2: "Mercier stability computation failed",
    3: "output capacity changed during evaluation",
    4: "invalid argument",
}


def _raise_for_status(status, path):
    reason = _STATUS_MESSAGES.get(status, "unknown GLISS status")
    raise RuntimeError(
        f"gliss_mercier_profile({os.fspath(path)!r}) failed with status "
        f"{status}: {reason}"
    )


def _resolution(value, name):
    if isinstance(value, (bool, np.bool_)):
        raise TypeError(f"{name} must be an integer")
    try:
        result = operator.index(value)
    except TypeError as error:
        raise TypeError(f"{name} must be an integer") from error
    if not 1 <= result <= _C_INT_MAX:
        raise ValueError(f"{name} must be between 1 and {_C_INT_MAX}")
    return result


def _export_path(path):
    try:
        value = os.fspath(path)
    except TypeError as error:
        raise TypeError("path must be a string or path-like object") from error
    if not isinstance(value, str):
        raise TypeError("path must resolve to a string")
    if "\0" in value:
        raise ValueError("path contains a null byte")
    export = Path(value)
    if not export.exists():
        raise FileNotFoundError(f"equilibrium export does not exist: {export}")
    if not export.is_file():
        raise ValueError(f"equilibrium export is not a file: {export}")
    encoded = os.fsencode(export)
    if len(encoded) > _C_INT_MAX:
        raise ValueError("encoded path is too long for the GLISS C API")
    return export, encoded


def _bind(library):
    library.gliss_mercier_profile.argtypes = (
        ctypes.c_char_p,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int),
    )
    library.gliss_mercier_profile.restype = None


def mercier_profile(path, n_theta=64, n_zeta=64):
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
    export, encoded = _export_path(path)
    library = _load_library()
    _bind(library)

    surfaces = ctypes.c_int()
    status = ctypes.c_int()
    library.gliss_mercier_profile(
        encoded,
        len(encoded),
        n_theta,
        n_zeta,
        0,
        ctypes.byref(surfaces),
        None,
        None,
        ctypes.byref(status),
    )
    if status.value not in (0, 3):
        _raise_for_status(status.value, export)
    count = surfaces.value
    if not 0 <= count <= _C_INT_MAX:
        raise RuntimeError(f"GLISS returned an invalid surface count: {count}")
    if count == 0:
        empty = np.empty(0, dtype=np.float64)
        return empty, empty.copy()

    s_values = (ctypes.c_double * count)()
    d_mercier = (ctypes.c_double * count)()
    library.gliss_mercier_profile(
        encoded,
        len(encoded),
        n_theta,
        n_zeta,
        count,
        ctypes.byref(surfaces),
        s_values,
        d_mercier,
        ctypes.byref(status),
    )
    if status.value != 0:
        _raise_for_status(status.value, export)
    if surfaces.value != count:
        raise RuntimeError("GLISS changed the surface count during evaluation")
    return (
        np.ctypeslib.as_array(s_values).copy(),
        np.ctypeslib.as_array(d_mercier).copy(),
    )


def mercier_objective(path, n_theta=64, n_zeta=64):
    """Return the most unstable Mercier discriminant in the profile."""
    _, d_mercier = mercier_profile(path, n_theta=n_theta, n_zeta=n_zeta)
    if d_mercier.size == 0:
        raise RuntimeError("GLISS returned an empty Mercier profile")
    return float(np.max(d_mercier))
