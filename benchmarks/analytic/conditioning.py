"""Diagnose slow eigenvalue roundoff on dumped production FEEC pencils.

Run after:
  mkdir -p /absolute/output
  fo exec benchmark_theta_pinch /absolute/output/spectrum.csv dump /absolute/output
  python3 benchmarks/analytic/conditioning.py /absolute/output

Matrices and outputs stay in the supplied directory. No reference is fitted.
Requires NumPy, SciPy, mpmath; default 50 decimal digits for small-case checks.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

import mpmath as mp
import numpy as np
import scipy
import scipy.linalg as la

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('directory', type=Path)
parser.add_argument('--digits', type=int, default=50)
args = parser.parse_args()
if args.digits < 30:
    parser.error('Use at least 30 decimal digits')
mp.mp.dps = args.digits
va = 1 / (8 * mp.pi * mp.mpf('1e-7'))
cs = mp.mpf(250) / 3
cusp = va * cs / (9 * (va + cs))
root = mp.findroot(lambda z: (mp.besselj(2, z) - mp.besselj(4, z)) / 2, 4.2)
q = mp.mpf(1) / 9 + (2 * root)**2
summ = (va + cs) * q
fast = (summ + mp.sqrt(summ**2 - 4 * va * cs * q / 9)) / 2
first_slow = va * cs * q / (9 * fast)
report = dict(digits=args.digits, numpy=np.__version__, scipy=scipy.__version__,
              mpmath=mp.__version__, analytical_cusp=str(cusp),
              analytical_first_radial_slow=str(first_slow), cases=[])

for case in ('p1_n8', 'p4_n8', 'p4_n16', 'p4_n32'):
    path = args.directory / f'{case}.dat'
    lines = path.read_text().splitlines()
    n = int(lines[0])
    data = np.loadtxt(path, skiprows=1)
    k = data[:, 0].reshape((n, n), order='F')
    m = data[:, 1].reshape((n, n), order='F')
    if not (np.isfinite(k).all() and np.isfinite(m).all()):
        raise ValueError('Nonfinite matrix')
    entry = dict(case=case, dimension=n,
                 matrix_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                 stiffness_asymmetry=float(np.max(abs(k-k.T))),
                 mass_asymmetry=float(np.max(abs(m-m.T))))
    # Symmetrization explicitly excludes small assembly asymmetry as the cause.
    k = (k + k.T) / 2
    m = (m + m.T) / 2
    l = la.cholesky(m, lower=True)
    inv = la.solve_triangular(l, np.eye(n), lower=True)
    h = inv @ k @ inv.T
    h = (h + h.T) / 2
    d = 1 / np.sqrt(m.diagonal())
    entry.update(mass_condition=float(np.linalg.cond(m)),
                 scaled_mass_condition=float(np.linalg.cond(d[:, None]*m*d[None, :])),
                 stiffness_norm_over_cusp=float(la.norm(k, 2) / float(cusp)),
                 whitened_norm_over_cusp=float(la.norm(h, 2) / float(cusp)),
                 eps_whitened_norm=float(np.finfo(float).eps * la.norm(h, 2)),
                 variants=[])
    for scaled in (False, True):
        factor = d if scaled else np.ones(n)
        km = factor[:, None] * k * factor[None, :]
        mm = factor[:, None] * m * factor[None, :]
        for reverse in (False, True):
            kk = km[::-1, ::-1].copy() if reverse else km
            mass = mm[::-1, ::-1].copy() if reverse else mm
            chol = la.cholesky(mass, lower=True)
            for lower in (False, True):
                w, v = la.eigh(kk, mass, driver='gv', lower=lower)
                residual = kk @ v[:, 0] - w[0] * (mass @ v[:, 0])
                dual = la.solve_triangular(chol, residual, lower=True)
                entry['variants'].append(dict(scaled=scaled, reverse=reverse,
                    lower=lower, minimum=float(w[0]),
                    minimum_minus_cusp=float(w[0] - float(cusp)),
                    mass_dual_residual=float(la.norm(dual))))
    if case != 'p4_n32':
        start = time.time()
        # Parse decimal dumps directly: do not re-round through Python floats.
        high_k = mp.matrix(n)
        high_m = mp.matrix(n)
        for index, line in enumerate(lines[1:]):
            a, b = line.split()
            i, j = index % n, index // n
            high_k[i, j], high_m[i, j] = mp.mpf(a), mp.mpf(b)
        high_k = (high_k + high_k.T) / 2
        high_m = (high_m + high_m.T) / 2
        inverse = mp.cholesky(high_m)**-1
        white = inverse * high_k * inverse.T
        values = mp.eigsy((white + white.T) / 2, eigvals_only=True)
        entry.update(high_precision_minimum=str(values[0]),
                     high_precision_gap=str(values[0] - cusp),
                     high_precision_seconds=time.time() - start)
    report['cases'].append(entry)
    print(case, entry.get('high_precision_minimum', 'double-only large case'), flush=True)
    (args.directory / 'conditioning_report.json').write_text(
        json.dumps(report, indent=2) + '\n')
