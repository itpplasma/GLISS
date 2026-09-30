"""Analytic GPEC Solov'ev equilibrium (GPEC ``equil/sol.f``, ``sol_run``).

psi(R, Z) = psio - psifac [R^2 Z^2 / e^2 + (R^2 - R0^2)^2 / 4] with
psio = e F a^2 / (2 q0 R0), psifac = psio / (a R0)^2, constant F and
mu0 p = pfac (1 - psi_n), pfac = 2 psio^2 (e^2 + 1) / (a R0 e)^2.
The plasma boundary psi_n = 1 is (R^2 - R0^2) = 2 a R0 cos t,
Z = e a R0 sin t / R.  q is proportional to q0; q(axis) = q0.
"""
import numpy as np
from scipy.optimize import brentq

R0, A, ELONGATION, F = 1.0, 0.33, 1.6, 1.0
MU0 = 4.0e-7 * np.pi


def flux_scales(q0):
    psio = ELONGATION * F * A**2 / (2.0 * q0 * R0)
    return psio, psio / (A * R0) ** 2


def surface(psi_n, points=4096):
    """Points of the surface psi_n in polar angle about the axis (R0, 0)."""
    angle = 2.0 * np.pi * np.arange(points) / points
    radius = np.empty(points)
    for i, t in enumerate(angle):
        def level(r):
            big_r = R0 + r * np.cos(t)
            z = r * np.sin(t)
            return ((big_r * z) ** 2 / ELONGATION**2
                    + (big_r**2 - R0**2) ** 2 / 4.0 - psi_n * (A * R0) ** 2)
        radius[i] = brentq(level, 1e-12, 0.9)
    return R0 + radius * np.cos(angle), radius * np.sin(angle)


def safety_factor_and_flux(psi_n, q0):
    """q and toroidal flux F * int dA / R inside the surface psi_n."""
    _, psifac = flux_scales(q0)
    big_r, z = surface(psi_n)
    wave = np.fft.fftfreq(big_r.size, 1.0 / big_r.size)
    dr = np.real(np.fft.ifft(1j * wave * np.fft.fft(big_r)))
    dz = np.real(np.fft.ifft(1j * wave * np.fft.fft(z)))
    step = 2.0 * np.pi / big_r.size
    grad_psi = psifac * np.hypot(2 * big_r * z**2 / ELONGATION**2
                                 + big_r * (big_r**2 - R0**2),
                                 2 * big_r**2 * z / ELONGATION**2)
    q = F / (2 * np.pi) * np.sum(np.hypot(dr, dz) * step / (big_r * grad_psi))
    flux = F * np.sum(np.log(big_r) * dz) * step
    return q, flux


def profiles(q0, samples=80):
    """Normalized toroidal flux s, iota, pressure in Pa and psi_n samples."""
    psi_n = np.linspace(0.0, 1.0, samples + 1) ** 2
    values = [safety_factor_and_flux(x, q0) for x in psi_n[1:]]
    q = np.r_[q0, [v[0] for v in values]]
    flux = np.r_[0.0, [v[1] for v in values]]
    psio, _ = flux_scales(q0)
    pfac = 2 * psio**2 * (ELONGATION**2 + 1) / (A * R0 * ELONGATION) ** 2
    return flux / flux[-1], 1.0 / q, pfac * (1.0 - psi_n) / MU0, psi_n, flux[-1]


def boundary_harmonics(modes, points=256):
    """Cosine R and sine Z coefficients of the boundary in the angle t."""
    t = 2.0 * np.pi * np.arange(points) / points
    big_r = np.sqrt(R0**2 + 2 * A * R0 * np.cos(t))
    z = ELONGATION * A * R0 * np.sin(t) / big_r
    r_cos = np.fft.rfft(big_r).real / points * 2
    r_cos[0] /= 2
    z_sin = -np.fft.rfft(z).imag / points * 2
    return r_cos[: modes + 1], z_sin[: modes + 1]
