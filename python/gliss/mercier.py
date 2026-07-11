"""Mercier stability profile bindings (``gliss_mercier_profile`` C symbol)."""
import ctypes

import numpy as np

from . import _load_library

_STATUS_MESSAGES = {
    1: "failed to read equilibrium export file",
    2: "Mercier stability computation failed",
    4: "invalid argument (check path, n_theta, n_zeta)",
}


def _raise_for_status(status, path):
    reason = _STATUS_MESSAGES.get(status, "unknown gliss_mercier_profile status")
    raise RuntimeError(
        f"gliss_mercier_profile({path!r}) failed with status {status}: {reason}"
    )


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
    """Return ``(s, d_mercier)`` radial profiles from a GVEC/CAS3D export.

    ``d_mercier`` follows the D_Mercier sign convention: a positive value is
    Mercier-unstable.
    """
    library = _load_library()
    _bind(library)

    encoded_path = path.encode("utf-8")
    surfaces = ctypes.c_int()
    status = ctypes.c_int()
    library.gliss_mercier_profile(
        encoded_path, len(encoded_path), n_theta, n_zeta, 0,
        ctypes.byref(surfaces), None, None, ctypes.byref(status),
    )
    if status.value not in (0, 3):
        _raise_for_status(status.value, path)

    count = surfaces.value
    s_values = (ctypes.c_double * count)()
    d_mercier = (ctypes.c_double * count)()
    if status.value == 3:
        library.gliss_mercier_profile(
            encoded_path, len(encoded_path), n_theta, n_zeta, count,
            ctypes.byref(surfaces), s_values, d_mercier, ctypes.byref(status),
        )
        if status.value != 0:
            _raise_for_status(status.value, path)

    return (
        np.array(s_values, dtype=np.float64),
        np.array(d_mercier, dtype=np.float64),
    )


def mercier_objective(path, n_theta=64, n_zeta=64):
    """Return ``max(d_mercier)``: most positive value is most Mercier-unstable."""
    _, d_mercier = mercier_profile(path, n_theta=n_theta, n_zeta=n_zeta)
    return float(np.max(d_mercier))
