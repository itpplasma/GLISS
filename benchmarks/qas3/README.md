# QAS3: VMEC, TERPSICHORE and GLISS from public sources

The TERPSICHORE 1.2 QAS3 benchmark family (nfp = 3, CURTOR = −75 kA) is used
to compare GLISS with public reference codes:

- VMEC's `DMerc*` Mercier terms against `gliss.mercier_profile` on an export
  written by `gliss.convert_vmec`.
- TERPSICHORE's fixed-boundary (IVAC = 0, MODELK = 0) lowest eigenvalue
  against the certified GLISS replay of its FORT.23 matrices.
- The TERPSICHORE inertia count against the independent GLISS FEEC
  discretization, with sign and count compared only, because the
  normalizations differ.

## Reproduce

```sh
./build.sh                             # STELLOPT 2f181f0d VMEC2000, TERPSICHORE 04dcf9a
export GLISS_LIB=$PWD/../../build/libgliss_c.so PYTHONPATH=$PWD/../../python
./run_case.sh qas3_p100 1.0 1.0        # base case: pressure and current factor 1
./run_case.sh qas3_c050 1.0 0.5        # half the net current
INDEPENDENT=1 ./run_case.sh qas3_c075 1.0 0.75
```

`build.sh` needs gfortran, OpenBLAS, ScaLAPACK/OpenMPI and NetCDF-Fortran
(Debian: `libopenblas-dev libscalapack-openmpi-dev libopenmpi-dev
libnetcdff-dev`). TERPSICHORE is fetched from the GitHub mirror of the public
EPFL GitLab release (`TERPSICHORE_URL` overrides it). Each VMEC run takes a
few minutes, TERPSICHORE about 15 s, and the replay seconds. The independent
GLISS path takes 15–35 min on four cores.

## Results (68-mode N = 1 family, m ≤ 8, |n| ≤ 5, ns = 65)

Mercier, median relative error over 0.1 < s < 0.9. The Mercier sign
convention agrees: positive is stable.

| Case | DShear | DCurr | DWell | DGeod | DMerc |
|---|---|---|---|---|---|
| DSHAPE tokamak | 0.00 % | 0.02 % | 0.02 % | 0.04 % | 0.01 % |
| W7-X β = 5 % | 0.1 % | 2.7 % | 0.1 % | 0.3 % | 2.2 % |
| QAS3 variants | ≤ 0.3 % | 1–1.4 % | ≤ 0.2 % | 0.4 % | 0.4 % |
| li383 | 0.1 % | 18 % | 1.8 % | 3.2 % | 17 % |

These are the medians before #33. After it, li383's minimum D_Mercier moved
from −7.39 to −0.118 (VMEC +0.022). W7-X and DSHAPE did not change.

Fixed-boundary stability against the current fraction:

| Current | TERPSICHORE λ | TERPSICHORE count | GLISS replay | GLISS independent |
|---|---|---|---|---|
| 1.00 | −7.03701e-7 | 5 | −7.0370098e-7, count 5 | unstable, count 1 |
| 0.75 | −8.7e-8 | 1 | agrees | stable (0) at ns = 65; unstable at ns = 257 (below) |
| 0.50 | +5.956134e-7 | 0 | +5.956091e-7, count 0 | stable (0) |

## Findings

- The replay reproduces TERPSICHORE's eigenvalue to 1e-8 relative, with
  mode overlap 0.99999999, in every unstable case, fixed and free boundary.
- Before [#32](https://github.com/itpplasma/GLISS/issues/32) the replay
  raised for stable cases and misreported IVAC > 0 input as a parity error.
- The independent count at base current rose from 1 to 12 with the widened
  magnetic-differential-equation table of #20 and returned to 1 with the
  cell-averaged resonant inverse of
  [#33](https://github.com/itpplasma/GLISS/issues/33), which also removed the
  QAS3 Mercier spike at ι = 2/3: minimum 0.1441 against VMEC's 0.1435.
- NaN Mercier terms on the zero-current case came from a truncated metric
  that is not positive definite
  ([#34](https://github.com/itpplasma/GLISS/issues/34)). They are now
  rejected at conversion and in Mercier; M = N = 16 gives a finite profile.
- GLISS Mercier spikes at rational ι
  ([#33](https://github.com/itpplasma/GLISS/issues/33)).
- TERPSICHORE's shifted inverse iteration does not converge when the shift
  AL0 is far from two close eigenvalues. At zero current and AL0 = −1e-5 it
  stops at NITMAX with nonconverged components and reports −3.24e-8, while
  the certified lowest is −3.54e-8 (overlap 0.58). `run_case.sh` warns when
  this happens. This is a usage limit, not a TERPSICHORE defect.
- The critical current fraction agrees once both radial resolutions are
  resolved. At 0.75 of the base current the lowest eigenvalue is close to
  marginal and depends on the VMEC radial resolution of the equilibrium as
  well as on the radial discretization of each code. GLISS independent
  (CAS3D marginality, 68 modes, 96 x 64 angular points), lowest eigenvalue
  by VMEC ns and FEEC radial cells:

  | VMEC ns | cells | degree 1 | degree 2 | TERPSICHORE |
  |---|---|---|---|---|
  | 33 | | | | +9.6e-8 |
  | 65 | 64 | +0.301 | +0.129 | −8.7e-8 |
  | 65 | 128 | +0.122 | | |
  | 65 | 256 | +0.061 | | |
  | 129 | 128 | +0.112 | +0.053 | −5.8e-8 |
  | 257 | 256 | −0.0037 | −0.149 | −2.7e-8 |

  The cell count is the `radial_cells` option, independent of the
  equilibrium surfaces; `compare.py refine EXPORT deck.data.modes 1 64 128
  256` reproduces a row on four threads in 2, 4 and 9 minutes.
  On the ns = 257 equilibrium both codes are unstable at 0.75. The earlier
  "GLISS stable" came from 64 cells on the ns = 65 equilibrium; refining
  only the cells there converges to a small positive value, so the
  equilibrium resolution decides the sign as much as the stability
  discretization. TERPSICHORE moves toward marginal from below with ns.
  Raising the VMEC angular resolution of the 0.75 equilibrium to
  mpol/ntor = 12/8 changes TERPSICHORE's eigenvalue by 1.6 % (it needs a
  Boozer table of M = 16, |N| = 12 there: a loop over the VMEC modes indexes
  the Boozer-table array, which overflows for mnmax_nyq = 326 with the
  benchmark's 14 x 10 table). The GLISS FORT.23 replay of the ns = 257 case
  (order 34,748, 287 MB FORT.23) assembles TERPSICHORE's matrices directly
  in block-tridiagonal storage and runs in 8.4 s and 434 MiB on four
  threads (dense storage would need 9.6 GB). TERPSICHORE's eigenvector has
  the quotient −2.82906511e-8 in GLISS, its own WP/WK to all printed
  digits; the certified lowest eigenvalue is −2.8365e-8 (overlap 0.9989),
  and the inertia counts two unstable modes, of which inverse iteration
  from one shift reports one.
  - The QAS D_curr difference below s = 0.4 (see li383 below for the
    mechanism) is not rechecked: its export needs Boozer resolution above
    M = N = 24 (Jacobian truncation 0.073), beyond this VM's memory on all
    127 surfaces.
- The li383 D_curr difference was VMEC's, not GLISS's. VMEC evaluates
  D_curr = −shear (⟨J·B/|∇φ|³⟩ − I′⟨B²/|∇φ|³⟩) with J·B from the curl of
  its Fourier-truncated field, while GLISS takes the parallel current from
  the magnetic differential equation of the force-balanced Boozer field.
  Raising the VMEC angular resolution of the same li383 input moves VMEC
  onto GLISS; GLISS is converged in Boozer harmonics (M = N = 16 and 24
  agree to four digits) and in radial surfaces (66 and 198):

  | VMEC mpol, ntor | VMEC D_curr at s = 0.25 | GLISS | median GLISS/VMEC D_curr | median D_Mercier |
  |---|---|---|---|---|
  | 9, 5 (input) | 0.0566 | 0.0685 | 1.102 | 1.098 |
  | 12, 8 | 0.0677 | 0.0672 | 0.985 | 1.018 |
  | 16, 10 | 0.0664 | 0.0679 | 0.999 | 1.025 |

  Close to the axis (s < 0.1) VMEC's D_curr also changes sign between
  ns = 99 and 199, while GLISS's does not. A 10 % local difference at
  s ≈ 0.66 persists at every resolution next to a low-order rational
  surface, where the ideal Pfirsch-Schlueter current is singular and each
  code regularizes it differently.
- `convert_vmec` used to reject converged VMEC 9.0 li383 and QAS outputs.
  Its |B| and current identities now share the Boozer Jacobian tolerance,
  and the averaged radial force balance is gated instead of the pointwise
  metric closure. li383 converts at M = N = 16; QAS needs more harmonics
  (the Jacobian truncation check reports 0.17 at 16).
