"""Fresh runner output against independent dispersion relations; no fitting."""
import argparse
from pathlib import Path
import csv
import json
import numpy as np
from scipy.special import jnp_zeros
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("directory", type=Path)
parser.add_argument("--legacy", action="store_true")
a = parser.parse_args()
prefix = "legacy" if a.legacy else "theta"
mu0 = 4 * np.pi * 1e-7  # Shared legacy SI convention, physical_constants.f90
va2, cs2, k2 = 1 / (mu0 * 2), (5 / 3) * 100 / 2, 1 / 9
q2 = (jnp_zeros(3, 2) / 0.5) ** 2 + k2
summ = (va2 + cs2) * q2
fast = (summ + np.sqrt(summ**2 - 4 * va2 * cs2 * k2 * q2)) / 2
slow = va2 * cs2 * k2 / (va2 + cs2)
x = np.genfromtxt(a.directory / f"{prefix}_spectrum.csv", delimiter=",", names=True)
if not all(np.isfinite(x[name]).all() for name in x.dtype.names):
    raise ValueError("Nonfinite spectrum or diagnostic")
rows = []
for p in range(1, 5):
    for n in (8, 16, 32):
        for parity in (1, 2):
            v = x[(x['degree'] == p) & (x['surfaces'] == n) &
                  (x['parity'] == parity)]['omega2_s_minus2']
            # Fixed branch partition by independent Alfven point, no nearest fit.
            f = v[v > va2 * k2 * 1.01]
            for name, measured, reference in zip(
                    ('slow_minimum', 'fast_1', 'fast_2'),
                    (v.min(), f[0], f[1]), (slow, *fast)):
                rows.append(dict(degree=p, surfaces=n, parity=parity, branch=name,
                                 omega2=measured, reference=reference,
                                 relative_error=abs(measured/reference - 1)))
rates = []
for p in range(1, 5):
    for branch in ('fast_1', 'fast_2'):
        r = [r for r in rows if r['degree'] == p and r['branch'] == branch
             and r['parity'] == 1]
        for coarse, fine in zip(r, r[1:]):
            rates.append(dict(degree=p, branch=branch,
                              coarse=coarse['surfaces'], fine=fine['surfaces'],
                              measured_order=float(np.log2(coarse['relative_error'] /
                                                           fine['relative_error'])),
                              smooth_asymptotic_order=2*p))
(a.directory / f'{prefix}_rates.json').write_text(json.dumps(rates, indent=2)+'\n')
with (a.directory / f'{prefix}_branches.csv').open('w') as out:
    w = csv.DictWriter(out, fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
colors = ['#0072B2', '#E69F00', '#009E73', '#CC79A7']  # Okabe-Ito
markers = ['o', 's', '^', 'D']
fig, axes = plt.subplots(1, 3, figsize=(12, 4.4))
for ax, branch, label in zip(axes, ('slow_minimum', 'fast_1', 'fast_2'),
                            ('Slow minimum: upper dense solve', 'First fast branch', 'Second fast branch')):
    for p, color, marker in zip(range(1, 5), colors, markers):
        for parity, style in ((1, '-'), (2, '--')):
            r = [r for r in rows if r['degree'] == p and r['branch'] == branch
                 and r['parity'] == parity]
            ax.loglog([r['surfaces'] for r in r], [r['relative_error'] for r in r],
                      color=color, marker=marker, linestyle=style,
                      fillstyle='full' if parity == 1 else 'none', label=f'p={p}' if parity == 1 else None)
    ax.set(title=label, xlabel='Radial intervals (dimensionless)', ylabel='Absolute relative error')
    ax.set_xticks([8, 16, 32], ['8', '16', '32'])
    ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.grid(alpha=.2)
    if a.legacy and branch == 'slow_minimum':
        ax.set_yscale('linear'); ax.set_ylim(0, .8)
fig.legend(*axes[2].get_legend_handles_labels(), loc='upper center',
           bbox_to_anchor=(.5, .92), ncol=4)
fig.suptitle('Legacy fixture: invalid straight-cylinder comparison' if a.legacy else
             'Exact cylinder through production FEEC surface assembly: fast branches converge')
fig.text(.5, .01, 'B=1 T, p=100 Pa, density=2 kg/m³, a=0.5 m, L=6π m, (m,n)=(3,1); 8×4 angles (exact); 16×16 (legacy).\n'
         'Solid/filled: parity 1; dashed/open: parity 2. References: Bessel J₃′ roots and ideal-MHD dispersion. No filtering of slow minimum.',
         ha='center', fontsize=8)
fig.text(.5, .085, 'Slow drift: upper-triangle full-spectrum helper; ordering-sensitive roundoff, not the public certified solve.',
         ha='center', fontsize=8)
fig.tight_layout(rect=(0, .14, 1, .94))
for ext in ('pdf', 'png'): fig.savefig(a.directory / f'{prefix}_convergence.{ext}', dpi=180)

h = np.genfromtxt(a.directory / 'homogeneous.csv', delimiter=',', names=True)
t = h['angle_rad']; sound = (5 / 3) * 1e5 / 2
s = va2 + sound; root = np.sqrt(s*s - 4*va2*sound*np.cos(t)**2)
refs = [(s-root)/2, va2*np.cos(t)**2, (s+root)/2]
fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
errors = []
for i, (key, name) in enumerate(zip(h.dtype.names[1:], ('Slow', 'Alfvén', 'Fast'))):
    axes[0].plot(t*180/np.pi, refs[i], color=colors[i], label=f'{name}: analytical')
    axes[0].plot(t[::4]*180/np.pi, h[key][::4], linestyle='none', marker=markers[i],
                 fillstyle='none', color=colors[i])
    errors.append(float(np.max(np.abs(h[key]-refs[i]))/va2))
axes[0].set(xlabel='Angle between k and B (degrees)', ylabel='ω² (s⁻²)', title='Local homogeneous 3×3 kernel')
axes[0].legend(fontsize=8)
v = np.genfromtxt(a.directory / 'helical.csv', delimiter=',', names=True)
f = v['external_iota_fraction']; margin = v['normalized_vertical_margin']
axes[1].plot(f, 5*f-2, color=colors[0], label='Fu analytical limit: 5f−2')
axes[1].plot(f[::4], margin[::4], 's', fillstyle='none', color=colors[1], label='GLISS limit function')
axes[1].axhline(0, color='gray', linewidth=.7); axes[1].axvline(.4, color='gray', linestyle=':')
axes[1].set(xlabel='External rotational-transform fraction f', ylabel='Normalized vertical energy margin', title='Helical cylinder limit, elongation κ=2')
axes[1].legend(fontsize=8)
fig.suptitle('Local analytical controls — these do not validate the global FEEC solver')
fig.text(.5, .01, 'Homogeneous: |k|=1 m⁻¹, B=1 T, p=10⁵ Pa, density=2 kg/m³, γ=5/3; lines: formula, markers: fresh GLISS.\n'
         'Helical: Fu PPPL-3368 (1999), eqs. 15–16; negative/positive margins are unstable/stable controls.', ha='center', fontsize=8)
fig.tight_layout(rect=(0, .10, 1, .94))
for ext in ('pdf', 'png'): fig.savefig(a.directory / f'local_controls.{ext}', dpi=180)
assert max(errors) < 1e-12, errors
assert np.max(abs(margin-(5*f-2))) < 1e-14
(a.directory / f'{prefix}_metrics.json').write_text(json.dumps(dict(local_scaled_errors=errors,
    theta_max_error=max(r['relative_error'] for r in rows),
    maximum_relative_residual=float(x['relative_residual'].max()),
    maximum_mass_orthogonality=float(x['mass_orthogonality'].max()), theta_acceptance=('INVALID legacy geometry comparison' if a.legacy else
                      'Measured convergence; no new acceptance bound fitted')), indent=2)+'\n')
