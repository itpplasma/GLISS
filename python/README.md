# gliss

`gliss` is the Python interface to the Global Linear Ideal Stability Solver.
Version 0.0.1 provides the Mercier stability diagnostic for GVEC/CAS3D NetCDF
equilibrium exports. The Linux wheel contains the compiled Fortran library and
uses a small, hand-written ISO C binding; no `f90wrap` dependency is required.

## Installation

```sh
python -m pip install gliss
```

NumPy is the only required Python dependency. The binary wheel also contains
the native libraries needed by GLISS. Installing from source requires CMake,
Ninja, and a Fortran compiler. LAPACK, PkgConfig, and the NetCDF C library are
also required.

## Mercier profile

```python
from pathlib import Path

import gliss

s, d_mercier = gliss.mercier_profile(
    Path("equilibrium_export.nc"),
    n_theta=64,
    n_zeta=64,
)
worst = gliss.mercier_objective("equilibrium_export.nc")
```

`s` and `d_mercier` are one-dimensional NumPy `float64` arrays with one entry
per retained radial surface. GLISS uses the `D_Mercier` convention: positive
values are unstable. `mercier_objective` returns `max(d_mercier)`, so a larger
positive result is less stable.

The input must be a regular file in the GVEC/CAS3D export format. Angular
quadrature sizes must be positive integers. Invalid Python arguments raise
`TypeError`, `ValueError`, or `FileNotFoundError`; invalid exports and solver
failures raise `RuntimeError` with the native status code.

### Reusing an equilibrium

The development API for version 0.0.2 can load an export once and reuse its
native data:

```python
from pathlib import Path

from gliss import Equilibrium

with Equilibrium(Path("equilibrium_export.nc")) as equilibrium:
    s_64, d_64 = equilibrium.mercier_profile(n_theta=64, n_zeta=64)
    s_128, d_128 = equilibrium.mercier_profile(n_theta=128, n_zeta=128)
```

`Equilibrium.close()` releases the native allocation and is safe to call more
than once. The context manager calls it on exit. Operations on a closed object
raise `RuntimeError`. Several contexts may coexist; calls using the same
context must not overlap. Concurrent context creation also requires a
thread-safe NetCDF C library.

The returned arrays are independent, writable NumPy arrays owned by Python.
No Fortran allocation crosses the ABI. Native failures map to subclasses of
`gliss.GlissError`: `GlissIOError`, `GlissComputationError`,
`GlissCapacityError`, `GlissArgumentError`, `GlissAllocationError`, and
`GlissInternalError`. The one-shot `mercier_profile()` function uses the same
context internally and remains convenient for a single evaluation.

## SIMSOPT

Install the optional dependency and import the adapter explicitly:

```sh
python -m pip install "gliss[simsopt]"
```

```python
from gliss.simsopt import MercierPenalty

penalty = MercierPenalty("equilibrium_export.nc")
value = penalty.J()
```

`MercierPenalty` is a leaf `Optimizable`. It has no equilibrium degrees of
freedom and no analytic derivative in version 0.0.1. Update `export_path` after
an external equilibrium solve before evaluating it again.

## Native library and ABI

`gliss.version()` reports the loaded native version. The Python package checks
ABI version 1 before every native call and rejects an incompatible library.
Platform wheels load their bundled library automatically. Developers can test
a local build explicitly:

```sh
cmake -S . -B build -G Ninja
cmake --build build --target gliss_c
GLISS_LIB=$PWD/build/libgliss_c.so python -c \
  'import gliss; print(gliss.version())'
```

`GLISS_LIB` deliberately overrides the bundled library. This is a development
and debugging facility, not required for normal use.

The installed C header is `gliss.h`. ABI version 1 defines fixed numeric status
values, opaque equilibrium handles, caller-owned output arrays, `size_t`
capacities, and caller-provided error buffers. An error buffer may be omitted
only by passing both a null pointer and zero capacity. Destroy accepts a null
handle and clears a live handle, which makes cleanup idempotent when callers
retain one authoritative handle variable. Existing ABI-v1 symbols and status
values remain unchanged when new functions are added; an incompatible change
requires a new ABI version.
