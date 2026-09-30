#!/usr/bin/env python3
"""Exact vacuum energies frozen in test_scalar_potential_vacuum.f90.

On the torus eta = eta0 of toroidal coordinates (R0 = 3, a = 1) the field
grad Phi, Phi = sqrt(cosh eta - cos theta) P^1_{3/2}(cosh eta) cos(2 theta)
cos(phi), is harmonic and regular outside. Its normal flux per unit
(theta, phi) is g(theta) cos(phi); the script prints the cosine series of g
and the energy integral |grad Phi|^2 over (a) the whole exterior and (b) the
region between the torus and the coordinate torus of radius 1.6 (an ideal
wall, zero normal flux), whose field is a least-squares combination of P and
Q harmonics. Needs mpmath and NumPy.
"""
import mpmath as mp
import numpy as np

mp.mp.dps = 30
MAJOR, MINOR, WALL, M, N = 3, 1, mp.mpf("1.6"), 2, 1
A_T = mp.sqrt(MAJOR**2 - MINOR**2)
ETA0 = mp.asinh(A_T / MINOR)
ETA1 = mp.asinh(A_T / WALL)


def legendre(k, eta, kind):
    function = mp.legenp if kind == 0 else mp.legenq
    return mp.re(function(k - mp.mpf(1) / 2, N, mp.cosh(eta), type=3))


def harmonic(k, eta, theta, kind):
    # P grows and Q decays by e^(k eta); each is scaled to unit size on the
    # surface where it is largest so the least-squares columns stay O(1).
    scale = legendre(k, ETA0 if kind == 0 else ETA1, kind)
    return (mp.sqrt(mp.cosh(eta) - mp.cos(theta)) * legendre(k, eta, kind)
            / scale * mp.cos(k * theta))


def normal_derivative(k, eta, theta, kind):
    return mp.diff(lambda x: harmonic(k, x, theta, kind), eta)


def h_phi(eta, theta):
    return A_T * mp.sinh(eta) / (mp.cosh(eta) - mp.cos(theta))


samples = 128
thetas = [2 * mp.pi * (j + mp.mpf(1) / 2) / samples for j in range(samples)]
# Normal flux per d theta d phi through the torus, along -grad(eta).
flux = [-normal_derivative(M, ETA0, t, 0) * legendre(M, ETA0, 0)
        * h_phi(ETA0, t) for t in thetas]
coefficients = [2 * sum(f * mp.cos(k * t) for f, t in zip(flux, thetas)) / samples
                for k in range(24)]
coefficients[0] /= 2
exterior = abs(sum(harmonic(M, ETA0, t, 0) * legendre(M, ETA0, 0) * f for t, f in zip(thetas, flux))
               * 2 * mp.pi / samples * mp.pi)
basis = [(kind, k) for kind in (0, 1) for k in range(16)]
rows = np.array(
    [[float(normal_derivative(k, ETA0, t, kind)) for kind, k in basis] for t in thetas]
    + [[float(normal_derivative(k, ETA1, t, kind)) for kind, k in basis] for t in thetas])
target = np.array([float(-f / h_phi(ETA0, t)) for f, t in zip(flux, thetas)]
                  + [0.0] * samples)
weights, *_ = np.linalg.lstsq(rows, target, rcond=None)
residual = np.linalg.norm(rows @ weights - target) / np.linalg.norm(target)
annulus = abs(sum(sum(w * float(harmonic(k, ETA0, t, kind))
                      for w, (kind, k) in zip(weights, basis)) * float(f)
                  for t, f in zip(thetas, flux)) * float(2 * mp.pi / samples) * np.pi)
print("flux cosine coefficients:")
for k, c in enumerate(coefficients):
    print(f"    {k:2d} {float(c): .17e}")
print(f"exterior energy {float(exterior):.17e}")
print(f"annulus energy  {annulus:.17e}  (fit residual {residual:.1e})")
