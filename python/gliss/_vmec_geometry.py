"""Boozer-coordinate geometry used by the VMEC importer."""

from dataclasses import dataclass
from typing import Any, Dict, Tuple

import numpy as np

_TWO_PI = 2.0 * np.pi
_POSITION_FRAME = "xhat,yhat rotated by winding*zeta_B"
_VMEC_WINDING = -1
_EVEN_FIELDS = {
    "mod_B",
    "xhat",
    "Jac",
    "g_tt",
    "g_tz",
    "g_zz",
    "II_tt",
    "II_tz",
    "II_zz",
    "B_contra_t",
    "B_contra_z",
}
_ODD_FIELDS = {"yhat", "zhat", "g_st", "g_sz"}


@dataclass(frozen=True)
class ConvertedGeometry:
    s: np.ndarray
    profiles: Dict[str, np.ndarray]
    harmonics: Dict[str, Tuple[np.ndarray, np.ndarray]]
    residuals: Dict[str, float]


def _position_frame(
    radius: np.ndarray, angle: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """Return the VMEC position in GLISS's field-period rotating frame."""
    return radius * np.cos(angle), -radius * np.sin(angle)


def _array(source: Any, name: str, shape: Tuple[int, ...] = ()) -> np.ndarray:
    values = np.asarray(getattr(source, name))
    if shape and values.shape != shape:
        raise ValueError(
            f"booz_xform returned {name} with shape {values.shape}; expected {shape}"
        )
    if not np.issubdtype(values.dtype, np.number) or not np.all(np.isfinite(values)):
        raise ValueError(f"booz_xform returned invalid values for {name}")
    return np.asarray(values, dtype=np.float64)


def _mode_tensor(
    cosine: np.ndarray,
    sine: np.ndarray,
    poloidal: np.ndarray,
    toroidal: np.ndarray,
    nfp: int,
) -> np.ndarray:
    if np.any(np.rint(toroidal / nfp) * nfp != toroidal):
        raise ValueError("booz_xform toroidal modes are not multiples of N_FP")
    n_modes = np.rint(toroidal / nfp).astype(np.int64)
    m_modes = np.rint(poloidal).astype(np.int64)
    if np.any(m_modes < 0):
        raise ValueError("booz_xform returned a negative poloidal mode")
    m_max = int(np.max(m_modes))
    n_max = int(np.max(np.abs(n_modes)))
    tensor = np.zeros((cosine.shape[1], m_max + 1, 2 * n_max + 1), dtype=np.complex128)
    for index, (mode_m, mode_n) in enumerate(zip(m_modes, n_modes)):
        tensor[:, mode_m, mode_n + n_max] = cosine[index] - 1j * sine[index]
    return tensor


def _radial_slope(coefficients: np.ndarray, s: np.ndarray, m: np.ndarray) -> np.ndarray:
    derivative = np.gradient(coefficients, s, axis=1, edge_order=2)
    power = 0.5 * m[:, None]
    regular = coefficients / s[None, :] ** power
    regular_slope = np.gradient(regular, s, axis=1, edge_order=2)
    regular_derivative = (
        s[None, :] ** power * regular_slope
        + power * s[None, :] ** (power - 1.0) * regular
    )
    use_regular = m[:, None] <= 8.0
    return np.where(use_regular, regular_derivative, derivative)


def _evaluate(
    coefficients: np.ndarray,
    theta: np.ndarray,
    zeta: np.ndarray,
    nfp: int,
    theta_order: int = 0,
    zeta_order: int = 0,
) -> np.ndarray:
    m = np.arange(coefficients.shape[1], dtype=np.float64)
    n_max = (coefficients.shape[2] - 1) // 2
    n = np.arange(-n_max, n_max + 1, dtype=np.float64)
    multiplier = (1j * m[None, :, None]) ** theta_order
    multiplier = multiplier * (-1j * nfp * n[None, None, :]) ** zeta_order
    theta_phase = np.exp(1j * np.outer(theta, m))
    zeta_phase = np.exp(-1j * nfp * np.outer(zeta, n))
    return np.einsum(
        "smn,tm,zn->stz",
        coefficients * multiplier,
        theta_phase,
        zeta_phase,
        optimize=True,
    ).real


def _series(source: Any, cosine_name: str, sine_name: str) -> Tuple[np.ndarray, ...]:
    modes = _array(source, "xm_b")
    toroidal = _array(source, "xn_b")
    surfaces = int(source.ns_b)
    expected = (modes.size, surfaces)
    cosine = _optional_component(source, cosine_name, expected)
    sine = _optional_component(source, sine_name, expected)
    return cosine, sine, modes, toroidal


def _optional_component(source: Any, name: str, shape: Tuple[int, int]) -> np.ndarray:
    values = _array(source, name)
    if values.shape == (0, 0):
        return np.zeros(shape, dtype=np.float64)
    if values.shape != shape:
        raise ValueError(
            f"booz_xform returned {name} with shape {values.shape}; expected {shape}"
        )
    return values


def _regular_tensor(
    source: Any,
    cosine_name: str,
    sine_name: str,
    s: np.ndarray,
    nfp: int,
) -> Tuple[np.ndarray, np.ndarray]:
    cosine, sine, modes, toroidal = _series(source, cosine_name, sine_name)
    values = _mode_tensor(cosine, sine, modes, toroidal, nfp)
    slope = _mode_tensor(
        _radial_slope(cosine, s, modes),
        _radial_slope(sine, s, modes),
        modes,
        toroidal,
        nfp,
    )
    return values, slope


def _angular_grid(source: Any, m_out: int, n_out: int) -> Tuple[np.ndarray, np.ndarray]:
    m_in = int(np.max(_array(source, "xm_b")))
    nfp = int(source.nfp)
    n_in = int(np.max(np.abs(_array(source, "xn_b"))) // nfp)
    n_theta = max(32, 2 * (m_in + m_out) + 2)
    n_zeta = max(32, 2 * (n_in + n_out) + 2)
    if n_theta > 512 or n_zeta > 512:
        raise ValueError("requested VMEC transform exceeds the 512-point angular limit")
    theta = _TWO_PI * np.arange(n_theta) / n_theta
    zeta = _TWO_PI * np.arange(n_zeta) / (nfp * n_zeta)
    return theta, zeta


def _project(
    values: np.ndarray,
    theta: np.ndarray,
    zeta: np.ndarray,
    nfp: int,
    m_max: int,
    n_max: int,
) -> Tuple[np.ndarray, np.ndarray]:
    m = np.arange(m_max + 1, dtype=np.float64)
    n = np.concatenate((np.arange(n_max + 1), np.arange(-n_max, 0))).astype(float)
    theta_factor = np.exp(1j * np.outer(theta, m))
    zeta_factor = np.exp(-1j * nfp * np.outer(zeta, n))
    coefficient = np.einsum(
        "stz,tm,zn->smn", values, theta_factor, zeta_factor, optimize=True
    ) / (theta.size * zeta.size)
    cosine = 2.0 * coefficient.real
    sine = -2.0 * coefficient.imag
    cosine[:, 0, 0] *= 0.5
    sine[:, 0, 0] = 0.0
    if n_max:
        cosine[:, 0, n_max + 1 :] = 0.0
        sine[:, 0, n_max + 1 :] = 0.0
    return cosine, sine


def _relative_max(actual: np.ndarray, expected: np.ndarray) -> float:
    scale = max(float(np.max(np.abs(expected))), np.finfo(np.float64).tiny)
    return float(np.max(np.abs(actual - expected)) / scale)


def _require_parity(name: str, cosine: np.ndarray, sine: np.ndarray) -> None:
    kept, rejected = (cosine, sine) if name in _EVEN_FIELDS else (sine, cosine)
    scale = max(float(np.max(np.abs(kept))), 1.0)
    if float(np.max(np.abs(rejected))) > 5.0e-10 * scale:
        raise ValueError(f"Boozer geometry violates stellarator symmetry in {name}")


def _reconstruct(
    pair: Tuple[np.ndarray, np.ndarray],
    theta: np.ndarray,
    zeta: np.ndarray,
    nfp: int,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    cosine, sine = pair
    coefficients = cosine - 1j * sine
    m = np.arange(cosine.shape[1], dtype=np.float64)
    n_max = (cosine.shape[2] - 1) // 2
    n = np.concatenate((np.arange(n_max + 1), np.arange(-n_max, 0))).astype(float)
    theta_phase = np.exp(-1j * np.outer(theta, m))
    zeta_phase = np.exp(1j * nfp * np.outer(zeta, n))

    def evaluate(multiplier: np.ndarray) -> np.ndarray:
        return np.einsum(
            "smn,tm,zn->stz",
            coefficients * multiplier,
            theta_phase,
            zeta_phase,
            optimize=True,
        ).real

    values = evaluate(np.ones((1, m.size, n.size)))
    theta_derivative = evaluate(1j * _TWO_PI * m[None, :, None])
    zeta_derivative = evaluate(-1j * _TWO_PI * n[None, None, :])
    return values, theta_derivative, zeta_derivative


def _force_balance_residual(
    harmonics: Dict[str, Tuple[np.ndarray, np.ndarray]],
    profiles: Dict[str, np.ndarray],
    s: np.ndarray,
    theta: np.ndarray,
    zeta: np.ndarray,
    nfp: int,
) -> float:
    reconstructed = {
        name: _reconstruct(pair, theta, zeta, nfp) for name, pair in harmonics.items()
    }
    jacobian = reconstructed["Jac"][0]
    contra_t, contra_t_t, contra_t_z = reconstructed["B_contra_t"]
    contra_z, contra_z_t, contra_z_z = reconstructed["B_contra_z"]
    g_st, g_st_t, g_st_z = reconstructed["g_st"]
    g_sz, g_sz_t, g_sz_z = reconstructed["g_sz"]
    beta_t = g_st_t * contra_t + g_st * contra_t_t
    beta_t += g_sz_t * contra_z + g_sz * contra_z_t
    beta_z = g_st_z * contra_t + g_st * contra_t_z
    beta_z += g_sz_z * contra_z + g_sz * contra_z_z
    slopes = {
        name: np.gradient(values, s, edge_order=2) for name, values in profiles.items()
    }
    toroidal_field = -slopes["Phi"][:, None, None]
    poloidal_field = -slopes["chi"][:, None, None] / nfp
    b_theta_slope = slopes["B_theta_avg"][:, None, None]
    b_zeta_slope = slopes["B_zeta_avg"][:, None, None]
    pressure_slope = slopes["p"][:, None, None]
    residual = 4.0e-7 * np.pi * pressure_slope * jacobian
    residual += toroidal_field * b_zeta_slope + poloidal_field * b_theta_slope
    residual -= toroidal_field * beta_z + poloidal_field * beta_t
    terms = np.abs(4.0e-7 * np.pi * pressure_slope * jacobian)
    terms += np.abs(toroidal_field * b_zeta_slope)
    terms += np.abs(poloidal_field * b_theta_slope)
    terms += np.abs(toroidal_field * beta_z)
    terms += np.abs(poloidal_field * beta_t)
    selected = (s >= 0.05) & (s <= 0.9)
    numerator = np.max(np.abs(residual[selected]), axis=(1, 2))
    denominator = np.maximum(
        np.max(terms[selected], axis=(1, 2)), np.finfo(np.float64).tiny
    )
    return float(np.max(numerator / denominator))


def convert_geometry(
    source: Any,
    beta_average: float,
    poloidal_max: int,
    toroidal_max: int,
) -> ConvertedGeometry:
    nfp = int(source.nfp)
    surfaces = int(source.ns_b)
    if nfp < 1 or surfaces < 5:
        raise ValueError("VMEC input must contain at least five half-grid surfaces")
    s = _array(source, "s_b", (surfaces,))
    indices = np.asarray(source.compute_surfs)
    available = int(source.ns_in)
    if (
        indices.shape != (surfaces,)
        or not np.issubdtype(indices.dtype, np.integer)
        or np.any(indices < 0)
        or np.any(indices >= available)
        or np.any(indices[1:] <= indices[:-1])
    ):
        raise ValueError("booz_xform returned invalid radial surface indices")
    indices = indices.astype(np.int64, copy=False)
    expected_s = (np.arange(surfaces) + 0.5) / surfaces
    if np.max(np.abs(s - expected_s)) > 1.0e-10:
        raise ValueError("VMEC half-grid is not uniformly spaced in normalized flux")
    theta, zeta = _angular_grid(source, poloidal_max, toroidal_max)
    radius, radius_s = _regular_tensor(source, "rmnc_b", "rmns_b", s, nfp)
    height, height_s = _regular_tensor(source, "zmnc_b", "zmns_b", s, nfp)
    angle, angle_s = _regular_tensor(source, "numnc_b", "numns_b", s, nfp)
    bmnc, bmns, modes, toroidal = _series(source, "bmnc_b", "bmns_b")
    gmnc, gmns, _, _ = _series(source, "gmnc_b", "gmns_b")
    mod_b = _mode_tensor(bmnc, bmns, modes, toroidal, nfp)
    gm = _mode_tensor(gmnc, gmns, modes, toroidal, nfp)

    def ev(coefficients: np.ndarray, dt: int = 0, dz: int = 0) -> np.ndarray:
        return _evaluate(coefficients, theta, zeta, nfp, dt, dz)

    r, rt, rz, rs = ev(radius), ev(radius, 1), ev(radius, 0, 1), ev(radius_s)
    z, zt, zz, zs = ev(height), ev(height, 1), ev(height, 0, 1), ev(height_s)
    p, pt, pz, ps = ev(angle), ev(angle, 1), ev(angle, 0, 1), ev(angle_s)
    rtt, rtz, rzz = ev(radius, 2), ev(radius, 1, 1), ev(radius, 0, 2)
    ztt, ztz, zzz = ev(height, 2), ev(height, 1, 1), ev(height, 0, 2)
    ptt, ptz, pzz = ev(angle, 2), ev(angle, 1, 1), ev(angle, 0, 2)
    phi_t, phi_z, phi_s = -pt, 1.0 - pz, -ps
    e_s = np.stack((rs, r * phi_s, zs), axis=-1)
    e_t = np.stack((rt, r * phi_t, zt), axis=-1)
    e_z = np.stack((rz, r * phi_z, zz), axis=-1)
    cross = np.cross(e_t, e_z)
    geometry_jacobian = np.einsum("...i,...i->...", e_s, cross)
    if np.any(r <= 0.0) or np.any(geometry_jacobian >= 0.0):
        raise ValueError("VMEC Boozer geometry is not a regular left-handed chart")
    normal = (
        cross
        * np.sign(geometry_jacobian)[..., None]
        / np.linalg.norm(cross, axis=-1)[..., None]
    )

    g_tt = np.einsum("...i,...i->...", e_t, e_t)
    g_tz = np.einsum("...i,...i->...", e_t, e_z)
    g_zz = np.einsum("...i,...i->...", e_z, e_z)
    g_st = np.einsum("...i,...i->...", e_s, e_t)
    g_sz = np.einsum("...i,...i->...", e_s, e_z)
    metric_determinant = g_tt * g_zz - g_tz**2
    if np.any(metric_determinant <= 0.0):
        raise ValueError("VMEC Boozer surface metric is singular")

    def second(
        r_ab: np.ndarray,
        z_ab: np.ndarray,
        phi_ab: np.ndarray,
        r_a: np.ndarray,
        r_b: np.ndarray,
        phi_a: np.ndarray,
        phi_b: np.ndarray,
    ) -> np.ndarray:
        vector = np.stack(
            (
                r_ab - r * phi_a * phi_b,
                r_a * phi_b + r_b * phi_a + r * phi_ab,
                z_ab,
            ),
            axis=-1,
        )
        return np.einsum("...i,...i->...", normal, vector)

    ii_tt = second(rtt, ztt, -ptt, rt, rt, phi_t, phi_t)
    ii_tz = second(rtz, ztz, -ptz, rt, rz, phi_t, phi_z)
    ii_zz = second(rzz, zzz, -pzz, rz, rz, phi_z, phi_z)
    phip = _array(source, "phip")
    if phip.shape != (available + 1,):
        raise ValueError("booz_xform returned phip on an unsupported radial grid")
    phip_half = 0.5 * (phip[indices] + phip[indices + 1])[:, None, None]
    gm_values = ev(gm)
    jacobian = phip_half * gm_values
    geometry_jacobian_residual = _relative_max(geometry_jacobian, jacobian)
    if geometry_jacobian_residual > 3.0e-2:
        raise ValueError(
            "VMEC geometry and Boozer Jacobian disagree: "
            f"{geometry_jacobian_residual:.6g}"
        )
    current_i = _array(source, "Boozer_I", (surfaces,))[:, None, None]
    current_g = _array(source, "Boozer_G", (surfaces,))[:, None, None]
    rotational_transform = _array(source, "iota", (available,))[indices]
    rotational_transform = rotational_transform[:, None, None]
    contra_t = rotational_transform / gm_values
    contra_z = 1.0 / gm_values
    field_strength = ev(mod_b)
    covariant_t = g_tt * contra_t + g_tz * contra_z
    covariant_z = g_tz * contra_t + g_zz * contra_z
    field_squared = covariant_t * contra_t + covariant_z * contra_z
    field_residual = _relative_max(field_squared, field_strength**2)
    current_scale = max(float(np.max(np.abs(current_g))), np.finfo(float).tiny)
    poloidal_current_residual = float(
        np.max(np.abs(covariant_t - current_i)) / current_scale
    )
    toroidal_current_residual = float(
        np.max(np.abs(covariant_z - current_g)) / current_scale
    )
    if field_residual > 5.0e-3:
        raise ValueError("VMEC Boozer metric and magnetic-field strength disagree")
    if poloidal_current_residual > 5.0e-3:
        raise ValueError(
            f"VMEC Boozer poloidal current identity failed: {poloidal_current_residual:.6g}"
        )
    if toroidal_current_residual > 5.0e-3:
        raise ValueError(
            f"VMEC Boozer toroidal current identity failed: {toroidal_current_residual:.6g}"
        )

    dtheta = -1.0 / _TWO_PI
    dzeta = -nfp / _TWO_PI
    xhat, yhat = _position_frame(r, p)
    pointwise = {
        "mod_B": field_strength,
        "xhat": xhat,
        "yhat": yhat,
        "zhat": z,
        "Jac": jacobian / (dtheta * dzeta),
        "g_tt": g_tt / dtheta**2,
        "g_tz": g_tz / (dtheta * dzeta),
        "g_zz": g_zz / dzeta**2,
        "g_st": g_st / dtheta,
        "g_sz": g_sz / dzeta,
        "II_tt": ii_tt / dtheta**2,
        "II_tz": ii_tz / (dtheta * dzeta),
        "II_zz": ii_zz / dzeta**2,
        "B_contra_t": dtheta * contra_t,
        "B_contra_z": dzeta * contra_z,
    }
    harmonics = {
        name: _project(values, theta, zeta, nfp, poloidal_max, toroidal_max)
        for name, values in pointwise.items()
    }
    for name, (cosine, sine) in harmonics.items():
        _require_parity(name, cosine, sine)

    phi = _array(source, "phi", (available + 1,))
    chi = _array(source, "chi", (available + 1,))
    pressure = _array(source, "pres", (available + 1,))[indices + 1]
    profiles = {
        "p": pressure,
        "B_theta_avg": -_TWO_PI * current_i[:, 0, 0],
        "B_zeta_avg": -_TWO_PI / nfp * current_g[:, 0, 0],
        "Phi": -0.5 * (phi[indices] + phi[indices + 1]),
        "chi": 0.5 * (chi[indices] + chi[indices + 1]),
        "iota": rotational_transform[:, 0, 0],
    }
    phi_slope = np.gradient(profiles["Phi"], s, edge_order=2)[:, None, None]
    chi_slope = np.gradient(profiles["chi"], s, edge_order=2)[:, None, None]
    residuals = {
        "toroidal_flux": _relative_max(
            -pointwise["Jac"] * pointwise["B_contra_z"], phi_slope
        ),
        "poloidal_flux": _relative_max(
            -nfp * pointwise["Jac"] * pointwise["B_contra_t"], chi_slope
        ),
        "boozer_jacobian": geometry_jacobian_residual,
    }
    residuals["force_balance"] = _force_balance_residual(
        harmonics, profiles, s, theta, zeta, nfp
    )
    if (
        residuals["toroidal_flux"] > 1.0e-2
        or residuals["poloidal_flux"] > 5.0e-2
        or residuals["boozer_jacobian"] > 3.0e-2
    ):
        raise ValueError(f"VMEC conversion failed field-identity checks: {residuals}")
    if not np.isfinite(beta_average):
        raise ValueError("VMEC volume-averaged beta is not finite")
    return ConvertedGeometry(s, profiles, harmonics, residuals)
