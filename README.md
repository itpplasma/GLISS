# GLISS

GLISS (Global Linear Ideal Stability Solver) computes the linear ideal-MHD
stability of three-dimensional toroidal equilibria with nested flux
surfaces.  It solves the energy-principle eigenvalue problem
`K x = omega^2 M x` with Fourier harmonics in the angles and spline finite
elements in the radius, reads equilibria from the
[GVEC](https://gitlab.mpcdf.mpg.de/gvec-group/gvec) CAS3D export, and is
built for exact differentiability: assembly kernels carry
Enzyme-generated derivative actions so that eigenvalues, marginal points,
and local criteria expose exact gradients to optimization loops.

The present code covers the equilibrium interface, fixed-boundary
two-component and compressible global operators, physical kinetic norm,
variable-block eigensolver, inertia certificates, Mercier diagnostics, and
selected free-boundary operators. Validation includes analytic cylinders,
Solov'ev marginality and the complete 191-mode QAS3 TERPSICHORE comparison.
Free-boundary plasma-vacuum parity and the complete cross-code benchmark table
remain under construction.

## Python

The Python package is the primary user interface. Install it with
`python -m pip install gliss`; version 0.0.1 exposes the validated Mercier
profile and objective through NumPy. Development main also provides reusable
`Equilibrium` and fixed-boundary `StabilityProblem` contexts with typed,
certified lowest-eigenpair results, opt-in full spectra with per-pair
diagnostics, deterministic full-spectrum run containers, and atomic versioned
equilibrium export.
The development API also exposes the shared two-component marginality
operator through an explicit general 3-D mode table and the axisymmetric
family used for the pinned Solov'ev comparison with DCON.
See the [Python guide](python/README.md) for examples, conventions, input and
output contracts, direct VMEC conversion, and the optional SIMSOPT adapter.

## Build

Requires CMake, Ninja, and a Fortran compiler. LAPACK, PkgConfig, and the
NetCDF C library are also required.

```sh
cmake -S . -B build -G Ninja
cmake --build build
ctest --test-dir build --output-on-failure
```

The Enzyme gradient gate needs matching Flang, `opt`, `llvm-link`, and
LLVMEnzyme versions:

```sh
cmake -S . -B build-enzyme -G Ninja \
  -DCMAKE_Fortran_COMPILER=flang-new \
  -DGLISS_ENABLE_ENZYME=ON \
  -DENZYME_PLUGIN=/path/to/LLVMEnzyme-22.so
cmake --build build-enzyme
ctest --test-dir build-enzyme -L enzyme --output-on-failure
```

## Diagnostic spectrum counts

`gliss_compatible_marginality` can count the generalized eigenvalues below
fixed physical shifts without assigning level numbers to mesh-dependent
continuum samples.  For example:

```sh
build/gliss_compatible_marginality equilibrium.nc 4 32 8 1e-8 1 \
  --physical-density=1e-7 \
  --count-shifts=0,0.001,0.002,0.005 \
  0,0 1,0 2,0
```

The command writes CSV with columns `shift,eigenvalues_below_shift`.  A row at
shift `s` is the inertia count of eigenvalues satisfying `lambda < s` for the
assembled generalized pencil.  Shifts must be finite, nonnegative, and
strictly increasing.  The count mode is mutually exclusive with bracketing
and `--inertia-only`; malformed or conflicting input exits nonzero before an
equilibrium is opened.  This interface is intended for deterministic
cross-grid and cross-code spectral-distribution comparisons.

## Certified diagnostic profiles

The same executable can reconstruct one radial eigenfunction after an inertia
bracket has isolated it:

```sh
build/gliss_compatible_marginality equilibrium.nc 4 64 8 1e-10 1 \
  --physical-density=1 \
  --stored-powers=0.5,0.5,0 \
  --eta-stored-powers=0.5,0.5,0 \
  --eigenvalue-bracket=0.003,0.004,0.005 \
  --eigenprofile-index=1 --profile-points=201 \
  0,1 1,1 2,1
```

The bracket endpoint counts must be `INDEX-1` and `INDEX`; otherwise the
command fails.  Profile mode automatically refines the eigenvalue interval to
at most `min(TOLERANCE,1e-10)` relative width.  It then uses deterministic
shift-invert iteration at the outer-bracket midpoint and emits a profile only
if the isolated level is the closest eigenvalue and the mass-whitened residual
proves an eigenspace-angle upper bound no larger than `1e-3`.

The selected-eigenpair CSV row reports the outer and refined brackets, their
certified midpoint, the independently computed Rayleigh quotient, the raw and
diagonally equilibrated action-relative residuals, a Frobenius-norm backward
error, the mass-whitened absolute residual, the eigenspace-angle bound, and the
reciprocal condition estimate of the equilibrated mass matrix.  These are
different diagnostics: a low mode of a strongly cancelling energy pencil can
have a visibly larger action-relative residual while remaining backward stable
and having a small, explicitly bounded subspace error.  The following rows are
`normal` and `eta` field values at cell-centred coordinates; no coefficient
layout is exposed to downstream scripts.

## Formulation and provenance

The formulation follows the CAS3D energy-principle programme published by
Carolin Schwab, later Carolin Nuehrenberg (one author): the 1991
dissertation and the 1993 formulation paper appeared under her maiden
name, the capability papers from 1996 on under her married name.  Further
methods derive from Bernstein et al. (1958) for the energy principle,
Newcomb (1960) and Suydam (1958) for the cylindrical gates, Mercier
(1960) and Landreman and Jorge (2020) for the interchange criterion, and
Anderson et al. (1990) for eigenvalue counting by matrix inertia.
`PROVENANCE.md` maps each module to its sources.

## License

MIT.  See [LICENSE](LICENSE).
