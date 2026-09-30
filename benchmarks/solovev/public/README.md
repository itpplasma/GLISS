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
(m=0..8 with both sidebands, FEEC degree 2). Each q0 takes about a minute on
four threads. `make_fixture.py` shrinks a `--stellsym` export into the 24 kB
fixtures under `test/data` (16 surfaces, M=8) used by
`test_solovev_axis_regularity`.

## Result

| q0 | DCON zero crossings | GLISS before #16: count, lowest omega^2 | GLISS now: count, lowest omega^2 |
| --- | --- | --- | --- |
| 1.035 | 1 | 1, -1.51e-3 | 1, -2.72e-4 |
| 1.039062 | 1 | 1, -1.37e-3 | 1, -1.98e-5 |
| 1.039843 | 0 | 1, -1.35e-3 | 0, +1.99e-5 |
| 1.045 | 0 | 1, -1.21e-3 | 0, +2.32e-4 |

These are `run.sh` defaults (64 surfaces, export M=24, m=0..8, degree 2) with
GLISS eigenvalues in the perpendicular-L2 marginality norm; one q0 takes about
three minutes on one thread. Use `THREADS=1` on a loaded machine: OpenMP
oversubscription slowed GVEC about a hundredfold here.

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

`convergence.sh OUT 1.035` exports the same GVEC state at ns = 16, 32 and 64
(M = 8) and solves the n = 1 family for FEEC degrees 1 to 4:

| ns | degree 1 | degree 2 | degree 3 | degree 4 |
|---|---|---|---|---|
| 16 | +6.98e-3 (stable) | −2.338e-4 | −3.867e-4 | −3.912e-4 |
| 32 | +1.83e-3 (stable) | −2.887e-4 | −3.056e-4 | −3.058e-4 |
| 64 | +3.47e-4 (stable) | −2.589e-4 | −2.638e-4 | −2.640e-4 |

Degree convergence at fixed ns is fast, but mesh convergence is first order,
and degree 1 misses the instability. The first radial element integrates
integrands that are smooth in sqrt(s) with a Gauss rule in s, and that
under-integration also hides an axis mode at degree 4
([#35](https://github.com/itpplasma/GLISS/issues/35)).

