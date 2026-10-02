# macOS wheel validation

The `macOS wheels` workflow builds separate x86-64 and arm64 wheels on
GitHub-hosted macOS 15 runners. Each repaired wheel targets macOS 15.0 or
later and is installed into fresh CPython 3.9–3.14 environments. The workflow
uploads wheel artifacts for review; publication and platform support claims
require successful hosted runs of the complete matrix.

The build uses GNU Fortran 14.3.0, Clang 19.1.7, NetCDF C 4.9.3 and
OpenBLAS 0.3.30. The two `macos-*.lock` files fix every native package by URL
and checksum, including compiler and OpenMP runtimes. Micromamba 2.9.0 installs
those files without solving dependencies again. Build artifacts retain the
installed package list, compiler version, Xcode version and operating system
version.

Delocate 0.13.0 copies native dependencies into each wheel and removes host
run paths. Validation runs outside the checkout with library-path overrides
unset. `check_macos_wheel.py` checks every bundled library's architecture,
dependency paths and run paths, then compiles and runs the repository's C
consumer against the installed header and library. `check_installed.py` runs
the Python suite, including native Solov'ev stability and Mercier oracles,
finite-difference derivative checks, persistence and executable quickstart
examples. The existing Linux CI and manylinux release script run separately.

To update the native toolchain, solve both target platforms with micromamba
using these direct specifications and review all changed packages:

```text
python=3.11.14
gfortran_osx-{64,arm64}=14.3.0
clang_osx-{64,arm64}=19.1.7
sdkroot_env_osx-{64,arm64}=15.0
libnetcdf=4.9.3
openblas=0.3.30
cmake=4.1.2
ninja=1.13.1
pkg-config=0.29.2
```

Use `CONDA_OVERRIDE_OSX=15.0 micromamba create --dry-run --json --no-rc
--platform osx-64` or `osx-arm64`, with `-c conda-forge` and a scratch prefix.
Write `@EXPLICIT` followed by each `actions.LINK` package URL and its `#md5`
checksum into the corresponding lock file. Repeat both hosted build and
installation matrices before accepting new locks. Runner labels and their
architectures follow [GitHub's runner reference](https://docs.github.com/en/actions/reference/runners/github-hosted-runners).
