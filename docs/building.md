# Building from source

## Native library and tests

The build requires CMake 3.22 or newer, Ninja, a Fortran compiler, BLAS/LAPACK,
PkgConfig and the NetCDF C library.

```sh
cmake -S . -B build -G Ninja
cmake --build build
ctest --test-dir build --output-on-failure
```

A clean single-configuration build defaults to the optimized `Release`
configuration, enables OpenMP assembly and prefers threaded OpenBLAS, falling
back to another BLAS/LAPACK provider when OpenBLAS is not installed. GLISS does
not set thread counts; the OpenMP and BLAS runtimes use their defaults.

The Python package is built from the same sources with
`python -m pip install .` (add `".[vmec]"` for VMEC conversion); see the
[Python guide](https://github.com/itpplasma/GLISS/blob/main/python/README.md).

## Release audit

Before a release, check the committed tree for compiler-generated Fortran
array temporaries and run the full test suite under the audited `-O3` build:

```sh
./ci/array_temporary_audit.sh
```

The script uses a detached temporary worktree and a private `fo` cache, so it
does not reconfigure the normal build tree. Set `GLISS_AUDIT_TMPDIR` to place
the temporary build on a large or fast filesystem.

## Enzyme derivative checks

The Enzyme gradient tests need matching versions of Flang, `opt`, `llvm-link`
and the LLVMEnzyme plugin:

```sh
cmake -S . -B build-enzyme -G Ninja \
  -DCMAKE_Fortran_COMPILER=flang-new \
  -DGLISS_ENABLE_ENZYME=ON \
  -DENZYME_PLUGIN=/path/to/LLVMEnzyme-22.so
cmake --build build-enzyme
ctest --test-dir build-enzyme -L enzyme --output-on-failure
```

## Documentation site

The Sphinx sources are in `docs/`. With `sphinx` and `myst-parser` installed,
build them with

```sh
make -C docs html
```

The `docs` workflow installs the wheel, runs every quickstart block against
it, builds the site, checks links, and publishes it to the `gh-pages` branch
with `ci/publish_docs.py`. `latest` follows the default branch; each release
tag `vX.Y.Z` keeps its own directory, selected with the version menu.
