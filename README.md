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

The present code covers the equilibrium interface (NetCDF reader,
coordinate adapter, Fourier reconstruction, integrals), the radial
B-spline basis, and a verified local plane-wave prototype with a dense
eigensolver and an opt-in Enzyme gradient gate. The global solver is under
construction.

## Python

The Python package is the primary user interface. Install it with
`python -m pip install gliss`; version 0.0.1 exposes the validated Mercier
profile and objective through NumPy. See the [Python guide](python/README.md)
for the API, sign convention, input contract, and optional SIMSOPT adapter.

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
