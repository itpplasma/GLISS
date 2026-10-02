# Solov'ev n=1 comparison from public sources

This directory regenerates the GPEC analytic Solov'ev family (GPEC
`equil/sol.f`, `sol_run`: R0=1 m, a=0.33 m, elongation 1.6, F=1 T m,
edge toroidal flux 0.7071679553076 Wb) with public GVEC 1.5.0 and evaluates
the fixed-boundary n=1 stability count with `gliss_axisymmetric`.
It needs no restricted research data.

```sh
cmake -S . -B build -G Ninja && cmake --build build
benchmarks/solovev/public/run.sh /tmp/solovev 1.035 1.039062 1.039843 1.045
```

`run.sh` creates a virtual environment, installs `gvec==1.5.0`, applies
`gvec-1.5.0-cas3d.patch`, and for each q0 writes the GVEC input
(`make_gvec_input.py`), solves the equilibrium (M=20, 32 degree-5 elements,
force tolerance 1e-10), exports 64 half-grid surfaces with M=24, and runs GLISS
(m=0..8 with both sidebands, FEEC degree 2). Angular quadrature is rounded
upward to a power of two that satisfies both the cubic position-forcing and
displacement-product bandwidth rules: the default M=24 export uses 256 x 8.
The CLI's `--angular NTHETA NZETA` option and its CSV output record the actual
grid without changing the equilibrium export or displacement mode range.
Each q0 takes about a minute on
four threads. `make_fixture.py` shrinks a `--stellsym` export into the 24 kB
fixtures under `test/data` (16 surfaces, M=8) used by
`test_solovev_axis_regularity`.

## Result

The following numbers and timings are historical results on the previous
fixed 64 x 8 grid. That grid is no longer admitted for a declared M=24
position table. The revised runner preserves the physical truncations and
uses 256 x 8; these archived values have not been rerun at that grid and do
not qualify the current operator. See [the angular grid contract](../../../docs/angular_grids.md).

| q0 | DCON zero crossings | GLISS before #16: count, lowest eigenvalue | GLISS after #16: count, lowest eigenvalue | GLISS now: count, lowest omega^2 (s^-2) |
| --- | --- | --- | --- | --- |
| 1.035 | 1 | 1, -1.51e-3 | 1, -2.72e-4 | 1, -342.3 |
| 1.039062 | 1 | 1, -1.37e-3 | 1, -1.98e-5 | 1, -18.12 |
| 1.039843 | 0 | 1, -1.35e-3 | 0, +1.99e-5 | 0, +45.28 |
| 1.045 | 0 | 1, -1.21e-3 | 0, +2.32e-4 | 0, +472.8 |

These used the former `run.sh` defaults (64 surfaces, export M=24, m=0..8,
degree 2, angular quadrature 64 x 8).
The first two GLISS columns used the earlier compatible perpendicular-L2
norm, a pure number. The current eigenvalues are omega^2 in s^-2 with the
perpendicular kinetic form at unit mass density (1 kg m^-3) as the inertia
([#35](https://github.com/itpplasma/GLISS/issues/35)); only the sign and the
count are compared with DCON. One q0 takes about 5 s on one thread (about
three minutes before the matrix-product assembly). Use `THREADS=1` on a
loaded machine: OpenMP oversubscription slowed GVEC about a hundredfold here.

DCON (GPEC `f5595c06`, reproduced independently from its public regression
inputs) is converged at q0=1.03956-1.03959; edge truncation only lowers it.
Before the fix, the GLISS lowest eigenvalue at q0=1.039843 was
-1.35e-3, -7.14e-4 and -3.58e-4 on 64, 128 and 256 surfaces: a spurious mode
proportional to the mesh width. Its cause was a non-conforming |m|=1 trial
space: a regular displacement has eta ~ s^(-1/2) at the axis, but eta used an
unweighted L2 basis, so the compression term diverged like 1/s and the reduced
quadrature hid it. Both components now carry the regular-displacement axis
factor, and the compressional unknown is the regular part of mu.

## Upstream GVEC defects worked around by the patch

These are not GLISS defects; they are tracked in GLISS until reported to GVEC.

1. `gvec.fourier.fft2d` executes `c[0, -N:] = 0` for N=0 toroidal modes,
   which zeroes the entire m=0 row: every axisymmetric export with
   `--MN_out M 0` loses its mean values (for example `mod_B_mnc[:, 0, 0] = 0`).
2. `pygvec to-cas3d` rotates `xhat, yhat` by the logical toroidal angle, not
   the Boozer angle zeta_B used as the export coordinate, so positions cannot
   be reconstructed on the Boozer grid. It also writes no frame metadata.
3. With `which_hmap=1`, positions are `x = R cos(zeta)`, `y = -R sin(zeta)`;
   the default `--winding 1` therefore produces `xhat = R cos(2 zeta)`, a
   double-covered frame. The frame that is periodic in zeta_B requires
   `--winding -1`.
4. The export omits the radial chart metric `g_st`, `g_sz`, which the drive
   and force-balance diagnostics need.
5. `--stellsym` assigns every field except `yhat`, `zhat` to the cosine
   parity; after item 4 the odd `g_st`, `g_sz` must use the sine parity.

The patch fixes 1, 2, 4 and 5; `run.sh` passes `--winding -1` for 3.

## Radial convergence

`convergence.sh OUT 1.035` exports the same GVEC state at ns = 16, 32, 64 and
128 (M = 8) and solves the n = 1 family (poloidal_max 6) for FEEC degrees 1
to 4. The marginality mass is the perpendicular kinetic form at unit mass
density, so eigenvalues are omega^2 in s^-2:

| ns | degree 1 | degree 2 | degree 3 | degree 4 |
|---|---|---|---|---|
| 16 | +6.880e3 (stable) | +1.713e2 (stable) | −3.1348e2 | −3.4398e2 |
| 32 | +2.633e3 (stable) | −2.3985e2 | −3.4255e2 | −3.4470e2 |
| 64 | +8.448e2 (stable) | −3.2886e2 | −3.4537e2 | −3.4548e2 |
| 128 | +9.58e1 (stable) | −3.4389e2 | −3.45663e2 | −3.45659e2 |

The conforming axis space ([#35](https://github.com/itpplasma/GLISS/issues/35))
ties the leading |m|=1 coefficients of xi^s and eta and integrates the axis
element in sqrt(s), so every degree converges monotonically from above.
Richardson extrapolation of degree 3 gives -345.69; its successive
differences fall by 10.3 and 9.8 (order 3.3), degree 2 approaches third
order, and degree 1 converges at order 1.3 to 1.4 and resolves the
instability only beyond ns = 128. Degrees 3 and 4 agree to 1e-5 at ns = 128.
The earlier space left the |m|=1 leading coefficients independent (a
logarithmically divergent compression energy hidden by the Gauss rule in s)
and normalized with a coefficient norm that is unbounded for regular |m|=1
displacements; its eigenvalues converged only at first order and, with an
exact axis rule, showed a spurious high-m axis mode at degree 4.
