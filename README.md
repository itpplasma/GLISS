# GLISS

GLISS (Global Linear Ideal Stability Solver) computes the linear ideal-MHD
stability of three-dimensional toroidal equilibria with nested flux
surfaces.  It solves the energy-principle eigenvalue problem
`K x = omega^2 M x` with Fourier harmonics in the angles and spline finite
elements in the radius, reads equilibria from the
[GVEC](https://gitlab.mpcdf.mpg.de/gvec-group/gvec) CAS3D export, and is
built for differentiability: verified assembly kernels carry
Enzyme-generated derivative actions. Fixed-boundary pressure-sample spectral
derivatives include the spline and stiffness response at fixed imported
geometry. The force-balanced equilibrium-to-spectrum design derivative chain
and optimization loop remain under construction.

Version 0.0.2 supports production fixed-boundary FEEC spectra and energies,
Mercier diagnostics, and symmetric GVEC or VMEC equilibrium input. Selected
free-boundary operators remain research components; they do not form a public
physical plasma-vacuum API. The TERPSICHORE FORT.23/24 entry points reproduce
that code's stored discretization for validation and are labeled compatibility
paths throughout the API and documentation.

## Documentation

The documentation is published at <https://itpplasma.github.io/GLISS/>:
`latest` follows the default branch and every release tag `vX.Y.Z` keeps its
own directory, selected with the version menu. `docs/` holds the Sphinx
sources (`make -C docs html`, with `sphinx` and `myst-parser` installed). The
`docs` workflow builds the site from the installed wheel after running every
quickstart block against it, checks links, and publishes it to the
`gh-pages` branch with `ci/publish_docs.py`; GitHub Pages serves that branch
once it is enabled in the repository settings.

## Python

The Python package is the primary user interface. Install it with
`python -m pip install gliss`; version 0.0.2 provides reusable
`Equilibrium` and fixed-boundary `StabilityProblem` contexts with typed,
certified lowest-eigenpair results, opt-in full spectra with per-pair
diagnostics, deterministic full-spectrum run containers, and atomic versioned
equilibrium export.
The API also exposes the shared two-component marginality
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

Release 0.0.2 provides a manylinux x86-64 wheel and a source distribution.
macOS wheels, the production free-boundary solve, and the complete
equilibrium-to-spectrum derivative chain are tracked as future work.
Asymmetric equilibria take the coupled parity operator, and precomputed
BOOZ_XFORM files convert with `convert_boozer`.

See [PLAN.md](PLAN.md) for the current audit findings, validation priorities,
and linked implementation issues.

## Build

Requires CMake, Ninja, a Fortran compiler, BLAS/LAPACK, PkgConfig, and the
NetCDF C library. A clean single-config build defaults to the optimized
`Release` configuration, enables OpenMP assembly, and prefers threaded
OpenBLAS. If OpenBLAS is not installed, CMake falls back to another available
BLAS/LAPACK provider. GLISS does not set a thread count: the OpenMP and BLAS
runtimes use their default thread counts.

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

## Diagnostic spectrum counts

`gliss_compatible_marginality` can count the generalized eigenvalues below
fixed physical shifts without assigning level numbers to mesh-dependent
continuum samples.  For example:

```sh
fo exec gliss_compatible_marginality equilibrium.nc 4 32 8 1e-8 1 \
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

## Diagnostic profiles

The same executable can reconstruct one radial eigenfunction after an inertia
bracket has isolated it:

```sh
fo exec gliss_compatible_marginality equilibrium.nc 4 64 8 1e-10 1 \
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
estimates an eigenspace-angle bound no larger than `1e-3`. These are
finite-pencil diagnostics with floating-point uncertainty.

The selected-eigenpair CSV row reports the outer and refined brackets, their
bracket midpoint, the independently computed Rayleigh quotient, the raw and
diagonally equilibrated action-relative residuals, a Frobenius-norm backward
error, the mass-whitened absolute residual, the eigenspace-angle bound, and the
reciprocal condition estimate of the equilibrated mass matrix.  These are
different diagnostics: a low mode of a strongly cancelling energy pencil can
have a visibly larger action-relative residual while remaining backward stable
and having a small, estimated subspace error.  The following rows are
`normal` and `eta` field values at cell-centred coordinates; no coefficient
layout is exposed to downstream scripts.

Counts and stand-alone bracket refinement use the production block storage. Profile,
external coefficient-energy and seeded-subspace diagnostics use the current
FEEC assembly with a conservative 2048-unknown dense limit. The default space
retains the Cartesian axis ties; an explicit eta-power table differing from
the normal table selects the historical unconstrained validation space.
See [mode diagnostics](docs/mode_diagnostics.md) for coefficient conventions,
branch integration and the remaining MISHKA/CASTOR acceptance requirements.

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
