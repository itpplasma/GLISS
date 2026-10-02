# Solov'ev free-boundary kink against DCON

The n=1 external kink of the GPEC analytic Solov'ev family at q0=1.5
(`benchmarks/solovev/public`: R0=1 m, a=0.33 m, elongation 1.6) is unstable
without a wall and is stabilized by a conformal ideal wall close enough to
the edge. The critical wall distance is a single number that both codes
compute from public inputs, so it tests the whole free-boundary chain:
plasma energy, edge flux and the vacuum with a wall.

```sh
cmake -S . -B build -G Ninja && cmake --build build
benchmarks/solovev/free_boundary/run.sh /tmp/solovev_free
```

`run.sh` builds GPEC through `benchmarks/solovev/dcon/run.sh`, reruns its
Solov'ev regression example with `vac_flag=t` and a VACUUM `ishape=6`
conformal wall at `a` plasma half-widths for a list of `a` (`a=20` is no
wall), and records DCON's total energy. It then solves the same analytic
equilibrium with GVEC 1.5.0 (48 surfaces, export M=16, as in
`benchmarks/solovev/public`) and bisects the GLISS critical wall with
`scan.py --bisect 0.12 0.20`: parity class 1, modes (0, 1) and m=1..MMAX
with n=+-1, FEEC degree 2, angular quadrature 128 x 8, zero floor 1e-8, and
the conformal wall of `gliss.VacuumModel` at `a` times the half-width
(R_max - R_min)/2 = 0.3526573415939913 m along the outward edge normal.
The default half-width is derived from the public GPEC shape parameters:
`(sqrt(R0^2+2*a*R0)-sqrt(R0^2-2*a*R0))/2`, with R0=1 m and a=0.33 m.

## Result

The GLISS numbers below are historical results obtained with 64 x 8 angular
quadrature and the approximate half-width 0.35245 m, before the cubic MDE
forcing admission rule was enforced and that normalization was corrected. That
grid is rejected for the regenerated export's declared M=16 position table:
the current rule requires more than 96 poloidal points. The runner now uses
128 x 8. These historical values and timings have not been rerun with that
policy and corrected normalization and do not certify the current operator;
retain them for comparison
with a fresh run. See [the angular grid contract](../../../docs/angular_grids.md).

DCON (GPEC `f5595c06`, m=-12..18, `mthvac=960`) changes the sign of its
total energy between a=0.152422 (+3.10e-3) and a=0.1525 (-5.69e-4), so
a_crit = 0.15249 by linear interpolation. The value is unchanged with
`mthvac` 480 and 1920, with `delta_mlow=delta_mhigh=16` (m=-20..34) and with
ODE tolerances 1e-8.

GLISS brackets of a_crit (bisection to 6e-4):

| edge mesh (nu x nv) | m <= 4 | m <= 8 | m <= 12 |
| --- | --- | --- | --- |
| 32 x 16 | 0.1506-0.1513 | | |
| 48 x 24 | 0.1519-0.1525 | 0.1488-0.1494 | |
| 64 x 32 | 0.1525-0.1531 | 0.1500-0.1506 | |
| 128 x 64 | 0.1544-0.1550 | 0.1513-0.1519 | 0.1513-0.1519 |
| 256 x 64 | | 0.1519-0.1525 | |
| 128 x 128 | | 0.1519-0.1525 | |
| 256 x 128 | | | 0.1519-0.1525 |

The poloidal range is converged at m <= 8 (m <= 12 gives the same bracket)
and m <= 4 overestimates a_crit by 2%. The radial discretization is
converged: 96 GVEC surfaces or FEEC degree 3 leave the 32 x 16, m <= 4
bracket unchanged. Refining the edge mesh raises a_crit at second order in
the mesh spacing, poloidally and toroidally alike. Bisected to 6e-5 at
m <= 8 on square meshes:

| edge mesh | a_crit |
| --- | --- |
| 64 x 64 | 0.15051-0.15057 |
| 128 x 128 | 0.15197-0.15203 |
| 256 x 256 | 0.15238-0.15244 |

The differences 1.46e-3 and 0.41e-3 shrink by 3.6 per halving (order 1.8).
Richardson extrapolation gives a_crit = 0.15255 at second order and 0.15257
at the observed order, against DCON's 0.15249: agreement to 5e-4 relative.

## Cost

The vacuum of an axisymmetric edge is block circulant over its nv toroidal
sectors and the n=1 data excite one sector harmonic, so the vacuum costs
one assembly of the first sector's rows (far panels by a degree-5 rule)
and one complex solve of size 4 nu with the wall. Both parity classes of
the plasma come from one pass of angular products, and the lowest
eigenvalue is certified from the inverse iterate with two inertia probes.
One wall distance at m <= 8 on the 128 x 128 edge, end to end through
`scan.py` on this 4-vCPU Xeon (2.1 GHz), takes 9.2 s on one thread and
4.0 s on four. DCON takes 6.3 s for one wall distance (m=-12..18), with or
without more threads.
