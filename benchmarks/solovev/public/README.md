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

| q0 | DCON zero crossings | GLISS count before fix | GLISS count, lowest omega^2 |
| --- | --- | --- | --- |
| 1.035 | 1 | 1 | 1, -2.72e-4 |
| 1.039062 | 1 | 1 | 1, -1.97e-5 |
| 1.039843 | 0 | 1 | 0, +1.99e-5 |
| 1.045 | 0 | 1 | 0 |

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
