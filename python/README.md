# gliss

`gliss` is the Python interface to the Global Linear Ideal Stability Solver.
Version 0.0.1 provides the Mercier stability diagnostic for GVEC/CAS3D NetCDF
equilibrium exports. Development main for version 0.0.2 also provides reusable
equilibrium contexts and fixed-boundary stability problems. The Linux wheel
contains the compiled Fortran library and uses a small, hand-written ISO C
binding; no `f90wrap` dependency is required.

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

## Fixed-boundary spectrum

This example assembles the physical compressible problem for two Fourier
harmonics, then certifies the lowest eigenpair in each parity class:

```python
from pathlib import Path

from gliss import Equilibrium, StabilityProblem

with Equilibrium(Path("equilibrium_export.nc")) as equilibrium:
    with StabilityProblem(
        equilibrium,
        modes=[(1, 1), (2, 1)],
        adiabatic_index=5.0 / 3.0,
        density_kg_m3=1.0,
        zero_floor=1.0,
        radial_quadrature="midpoint",
    ) as problem:
        result = problem.solve()

for parity in result.classes:
    print(
        parity.parity_class,
        parity.lowest_eigenvalue,
        parity.certificate,
        parity.negative_count,
    )
```

`modes` contains explicit `(m, n)` integer pairs. GLISS uses the Fourier phase
`2*pi*(m*theta - n*zeta/N_T)`, where `N_T` is the number of field periods.
Poloidal mode `m` must be nonnegative; an axis mode with `m=0` also requires
`n>=0`. Duplicate modes are rejected. `adiabatic_index` is nonnegative,
`density_kg_m3` is a positive SI mass density, and `zero_floor` is a positive
`omega^2` threshold in `s^-2`. Radial quadrature is `"midpoint"` or
`"gauss2"`. Angular quadrature is currently fixed at 64 by 64.

The assembled problem uses the physical compressible stiffness and mass,
transformed one-period Fourier assembly, P1/P0/P0 radial spaces and exact
fixed-edge elimination. It owns the assembled matrices; the source
`Equilibrium` may be closed after construction. Repeated calls to `solve()`
reuse those matrices.

`StabilityResult.classes` contains parity classes 1 and 2. Each
`SpectrumResult` reports the lowest computed `omega^2` in `s^-2`, counts below
and inside the configured zero floor, and the certificate components:
eigenpair residual, floating-point resolution and inertia interval. Negative
`omega^2` is unstable. The returned eigenvector satisfies `x.T @ M @ x = 1`.
It is a read-only `float64` array in original dynamic-layout order: fixed-edge
normal coefficients, `eta`, then compressional `mu`. The `normal`, `eta` and
`mu` properties return the corresponding views. If the zero floor contains
the entire spectrum, `has_eigenvector` is false and the eigenvector is empty.

The current API returns the certified lowest eigenpair for each class, not the
complete spectrum. `StabilityProblem.close()` is idempotent. Calls on one
problem must not overlap, but independently constructed problems may coexist.

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

The installed C header is `gliss.h`. `gliss.get_include()` returns its directory
in a wheel installation; CMake source installs place it under the configured
include prefix. ABI version 1 defines fixed numeric status values, opaque
equilibrium and stability-problem handles, caller-owned output arrays, `size_t`
capacities, and caller-provided error buffers. An error buffer may be omitted
only by passing both a null pointer and zero capacity. Destroy accepts a null
handle and clears a live handle, which makes cleanup idempotent when callers
retain one authoritative handle variable.

Initialize `gliss_spectrum_summary.struct_size` with `sizeof` before a solve.
The struct reports conventions and certificate metadata; the caller owns the
eigenvector buffer. Problem construction copies and assembles everything it
needs from its equilibrium, so the two handles have independent lifetimes.
Existing ABI-v1 symbols and status values remain unchanged when functions are
added; an incompatible signature, layout or numeric status change requires a
new ABI version.
