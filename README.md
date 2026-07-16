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
family used for the pinned Solov'ev comparison with DCON. A separate
CAS3D2MN phase-envelope entry point translates the ordered carrier/envelope
table and calls the same production assembly and eigensolver. The default
physical-L2 norm canonicalizes coincident Fourier modes. The explicit Schwab
coefficient norm instead pulls that physical operator back to every labeled
envelope coefficient, retaining the exact redundant zero-stiffness directions
and evaluating inertia on the physical quotient.
Paired TERPSICHORE FORT.23/FORT.24 files from a `MODELK=0` pressureless-
pseudoplasma run can be solved through the same public Python package, with
the stored TERPSICHORE mode available for direct diagnostic comparison. This
dense same-basis compatibility path is a validation tool; it is not the
production physical plasma-vacuum interface.
See the [Python guide](python/README.md) for examples, conventions, input and
output contracts, direct VMEC conversion, and the optional SIMSOPT adapter.

## Build

Requires CMake, Ninja, a Fortran compiler, BLAS/LAPACK, PkgConfig, and the
NetCDF C library. A clean single-config build defaults to the optimized
`Release` configuration and prefers threaded OpenBLAS. If OpenBLAS is not
installed, CMake falls back to another available BLAS/LAPACK provider.
The provider controls its thread count; for example,
`OPENBLAS_NUM_THREADS=16` uses 16 OpenBLAS threads and
`OPENBLAS_NUM_THREADS=1` selects deterministic one-thread controls.

```sh
cmake -S . -B build -G Ninja
cmake --build build
ctest --test-dir build --output-on-failure
```

Before a release, audit the committed tree for compiler-generated Fortran
array temporaries and run the complete test suite under the audited `-O3`
build:

```sh
./ci/array_temporary_audit.sh
```

The script uses a detached temporary worktree and a private `fo` cache, so it
does not reconfigure the normal build tree.  Set `GLISS_AUDIT_TMPDIR` to place
the temporary build on a large or fast filesystem.

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
