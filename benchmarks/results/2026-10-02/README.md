# GLISS validation measurements, October 2, 2026

These are fresh executions on a shared Ryzen 9 5950X Linux workstation with
GNU Fortran 16.2.1. Each artifact records its own source revision, tracked patch
digest and build provenance; the measurements are snapshots of those revisions.
Wall times include other workstation activity and are single observations.

## Independent DCON sign comparison

The reference runner built public GPEC DCON at
`f5595c0689c624834d720068ddc8b5a7e4028248` and evaluated five Solov'ev equilibria
with default and tighter tolerances. See `dcon_solovev.csv` and
`dcon_provenance.txt`. Default tolerances missed crossings at q0=1.030 and
1.035; tighter settings recovered them. This upstream defect is reported in
[GPEC #292](https://github.com/PrincetonUniversity/GPEC/issues/292).

The native Python refinement runner separately evaluated the committed public
GVEC exports at q0=1.035 and 1.045 using degree 3, poloidal_max=6, one OpenMP
thread and one BLAS thread. Both equilibria agree with tighter DCON's sign on
every tested mesh. Values below use GLISS's physical mass normalization.

| q0 | radial cells | negative inertia | lowest eigenvalue |
| --- | ---: | ---: | ---: |
| 1.035 | 32 | 1 | -343.559311 |
| 1.035 | 64 | 1 | -345.727881 |
| 1.035 | 128 | 1 | -345.843301 |
| 1.045 | 32 | 0 | 458.915732 |
| 1.045 | 64 | 0 | 453.726732 |
| 1.045 | 128 | 0 | 453.368212 |

`gliss-native-refinement-final.json` retains input/library hashes, certificates,
solver diagnostics, timing and peak process RSS. Sign agreement is the common
observable. DCON and GLISS raw eigenvalues have different normalizations, and
these data do not establish a precision interval or a converged marginal q0.

Reproduce the native comparison after building the shared library:

```sh
GLISS_LIB="$PWD/build/libgliss_c.so" PYTHONPATH=python \
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  python3 benchmarks/solovev/native_refinement.py refinement.json
bash benchmarks/solovev/dcon/run.sh "$PWD/dcon-results" \
  1.030 1.035 1.039062 1.039843 1.045
```

## Analytical families

`benchmarks/analytic/run.sh` evaluated homogeneous and helical local limits,
the retained legacy model, and the current theta-pinch family. Compact CSV,
metrics, rate estimates and complete build provenance are stored beside this
report. The three scaled local-limit errors range from 7.31e-17 to 1.02e-15.
The maximum spectrum residual is 1.16e-13 and maximum mass orthogonality defect
is 1.14e-14. Degree 1's first fast branch converges at order approximately 2;
degree 2's at order approximately 4. Higher-order rates encounter finite
precision limits, as shown in the retained rate tables. No new acceptance
tolerance was fitted to these observations.

## Remaining comparisons

The original CAS3D W7-X deck and its complete normalization are unavailable.
The existing L139 artifact measures sparse scaling on a different equilibrium;
it cannot establish same-deck agreement or superiority. Convention-complete
MISHKA/CASTOR mode fixtures, a qualified figure-eight equilibrium, a fresh
cross-code stability study on that equilibrium, and the complete force-balanced
GVEC design derivative remain outstanding. Issue #35's degree-1 detection
criterion at 64 surfaces also remains unmet. None of these gaps is counted as
a passed benchmark.
