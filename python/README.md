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

VMEC conversion is optional:

```sh
python -m pip install "gliss[vmec]"
```

This extra installs SciPy and `booz_xform`. They are used only while converting
a VMEC equilibrium; loading and solving the converted file still requires only
the base package.

## VMEC input

Convert a standard, converged VMEC `wout` file before constructing an
equilibrium:

```python
from pathlib import Path

import gliss

converted = gliss.convert_vmec(
    Path("wout_W7X.nc"),
    Path("W7X_gliss.nc"),
    poloidal_max=7,
    toroidal_max=7,
)
with gliss.Equilibrium(converted) as equilibrium:
    s, d_mercier = equilibrium.mercier_profile()
```

The conversion runs the maintained `booz_xform` implementation, then writes
the same left-handed, one-field-period Boozer convention used by pyGVEC's
CAS3D exporter. The output contains SI pressure, magnetic field and fluxes;
`s` is normalized toroidal flux. The Fourier phase is
`2*pi*(m*theta - n*zeta)`, with `zeta` normalized over one field period.

The importer currently accepts fixed-boundary, stellarator-symmetric VMEC
files with `signgs=-1` and a successful `ier_flag`. It rejects asymmetric and
reversed-field-pinch files instead of applying the symmetric GLISS operator to
an incompatible chart. Mode limits are integers from 0 through 64. The
default 7 by 7 truncation is deliberate: higher modes can amplify radial
noise in a surface-by-surface Boozer transform. Increase them only after a
convergence check.

The generated file supports GLISS's fixed-boundary operators. It does not
contain a vacuum-region mesh or conducting-wall model for free-boundary work.

Before writing, GLISS checks the chart orientation, surface metric, Boozer
Jacobian, magnetic-field strength, flux identities, currents, symmetry and
force balance. A conversion with an interior normalized force-balance
residual above 0.1 is rejected. This is a corruption and gross-resolution
gate, not a literature-validation tolerance. The residuals, source basename
and Boozer transform resolution are stored as NetCDF attributes so validation
runs can reproduce the exact preprocessing choice.

An existing output is never replaced unless `overwrite=True`. Writes use a
temporary file in the destination directory and an atomic replacement; a
failed conversion removes the temporary file. Install-time dependency errors,
invalid options, failed VMEC solves and numerical identity failures produce
distinct Python exceptions with the rejected condition in the message.

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

### Writing an equilibrium export

Legacy GVEC/CAS3D exports have schema version 0. GLISS writes schema version 1:

```python
from pathlib import Path

import gliss

source = Path("equilibrium_export.nc")
normalized = Path("equilibrium_gliss_v1.nc")
with gliss.Equilibrium(source) as equilibrium:
    print(equilibrium.schema_version)
    equilibrium.write(normalized)

with gliss.Equilibrium(normalized) as equilibrium:
    assert equilibrium.schema_version == 1
```

Version 1 names the `gvec-cas3d-export` schema in the NetCDF global
attributes. It stores the radial grid, mode order, symmetry, field periods,
winding, average beta, six profiles, thirteen required harmonic pairs and the
two optional chart-metric pairs. The writer preserves every value represented
by `Equilibrium`; it does not copy unrecognized variables or attributes from
the source file.

The Python writer creates a temporary file in the destination directory and
replaces the destination only after the native writer closes a complete
NetCDF file. A failed write leaves an existing destination unchanged and
removes temporary output. The low-level C function
`gliss_equilibrium_write()` uses exclusive creation and refuses to overwrite
an existing path.

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

The historical stopping controls are explicit and immutable:

```python
controls = gliss.SolverTolerances(
    eigenvalue_relative=1.0e-13,
    residual_relative=1.0e-12,
    negative_bracket_relative=1.0e-9,
    negative_bracket_floor=1.0e-3,
    inverse_iteration_limit=500,
    bracket_iteration_limit=200,
)

with gliss.Equilibrium("equilibrium_export.nc") as equilibrium:
    with gliss.StabilityProblem(
        equilibrium,
        modes=[(1, 1), (2, 1)],
        solver_tolerances=controls,
    ) as problem:
        result = problem.solve()
```

Inverse iteration stops when the eigenvalue change is no larger than
`max(eigenvalue_relative * max(1, abs(omega_squared)), resolution)` and the
backward residual is no larger than the analogous `residual_relative` bound.
The negative-eigenvalue inertia bracket stops at
`negative_bracket_relative * abs(midpoint) + negative_bracket_floor *
zero_floor`.  The iteration limits bound those two loops.  All tolerances must
be finite and positive; limits are signed 32-bit integers of at least one.
Changing them does not change the matrices, radial discretization, floor-band
classification or normalization.  Omitting them selects the exact constants
used before this API was exposed.  Each `SpectrumResult` reports the controls
that produced it.

Evaluate the physical terms before closing the assembled problem:

```python
with gliss.Equilibrium("equilibrium_export.nc") as equilibrium:
    with gliss.StabilityProblem(equilibrium, modes=[(1, 1), (2, 1)]) as problem:
        mode = problem.solve_class(parity_class=1)
        energy = problem.energy(1, mode.eigenvector)

print(energy.components)
print(energy.potential_energy, energy.kinetic_energy)
print(energy.rayleigh_quotient, mode.lowest_eigenvalue)
```

`EnergyTerms.components` contains field-line bending, magnetic shear,
magnetic compression, pressure drive and plasma compressibility, in that
order. Every term except pressure drive is a positive quadratic form. The
pressure-drive term can have either sign. GLISS reports the raw forms
`x.T @ K @ x` and `x.T @ M @ x`, without a factor of one half. Their quotient
is `omega^2`. For a solver eigenvector, the kinetic form is one and the
quotient agrees with `lowest_eigenvalue` within its numerical certificate.

`potential_energy` is evaluated from the assembled total operator rather than
constructed from the five reported terms. `closure_error` is their absolute
difference, and `closure_tolerance` is the floating-point acceptance bound.
Python rejects a result that does not close, has a nonpositive kinetic form,
contains a nonfinite value or makes a nominally positive term significantly
negative. Input vectors must be finite, one-dimensional and contain exactly
the parity class's number of unknowns. Noncontiguous array inputs are copied
to native `float64`; returned `EnergyTerms` objects own no native memory.

Request every eigenpair only when the dense cost is acceptable:

```python
with gliss.Equilibrium("equilibrium_export.nc") as equilibrium:
    with gliss.StabilityProblem(equilibrium, modes=[(1, 1), (2, 1)]) as problem:
        spectrum = problem.solve_full_spectrum_class(parity_class=1)

print(spectrum.eigenvalues)
print(spectrum.residuals)
print(spectrum.eigenvectors.shape)
```

`eigenvalues`, `rayleigh_quotients`, `residuals`, and `resolutions` are
read-only `float64` arrays with one entry per unknown. `eigenvectors` has shape
`(eigenpair, component)`. Its `normal`, `eta`, and `mu` views retain that first
axis. Eigenvalues, Rayleigh quotients, residuals, and resolutions use `s^-2`.
Each row satisfies the same `x.T @ M @ x = 1` normalization as the certified
pair.

The complete array includes modes inside `[-zero_floor, zero_floor]`. The
certified solver skips this numerical floor band when there is no unstable
direction. `spectrum.certified_index` identifies its independently certified
pair in the full array, or is `None` when the floor band contains the complete
spectrum. Inertia counts and the active eigenvalue are checked before Python
returns the result.

`solve_full_spectrum()` returns both parity classes. It runs the certified
block solve and then a dense LAPACK solve for each class. The dense stage costs
`O(unknowns^3)` time and `O(unknowns^2)` memory; use `solve()` or
`solve_class()` when only the stability margin is needed. NumPy and native
allocation failures raise typed GLISS exceptions. `StabilityProblem.close()`
is idempotent. Calls on one problem must not overlap, but independently
constructed problems may coexist.

### Configuration, results, and run manifests

`StabilityConfiguration` records the inputs independently of a native
context. Results and configurations use deterministic, versioned JSON:

```python
from pathlib import Path

import gliss

configuration = gliss.StabilityConfiguration(
    modes=[(1, 1), (2, 1)],
    adiabatic_index=5.0 / 3.0,
    density_kg_m3=1.0,
    zero_floor=1.0,
    radial_quadrature="midpoint",
)
configuration.write("configuration.json")

with gliss.Equilibrium(Path("equilibrium_export.nc")) as equilibrium:
    with configuration.create_problem(equilibrium) as problem:
        result = problem.solve()
        result.write("result.json")
        problem.write_manifest("run.json", result)

restored_configuration = gliss.StabilityConfiguration.read(
    "configuration.json"
)
restored_result = gliss.StabilityResult.read("result.json")
manifest = gliss.RunManifest.read("run.json")
manifest.verify_equilibrium("equilibrium_export.nc")
```

Configuration schema `gliss.stability.configuration`, version 1, records the
fixed boundary, mode pairs, physical scalars and radial quadrature. Result
schema `gliss.stability.result`, version 1, stores both parity classes with all
reported conventions, certificate terms and read-only eigenvectors. A round
trip preserves every binary64 value. Rewriting an unchanged object produces
the same bytes. This JSON schema stores the certified active pair only. It is
unchanged by the separate full-spectrum format.

Run schema `gliss.stability.run`, version 1, embeds the configuration and
result. It records the equilibrium export format, base filename, byte count
and SHA-256, including equilibrium schema 0 or 1. It also records the GLISS
Python/native versions and ABI, plus the NumPy and Python versions. Absolute
input paths are not stored. The manifest contains the complete run contract
and output without exposing a private directory; the NetCDF input remains a
separate checksummed artifact.

`StabilityProblem.write_manifest()` records the schema, size and checksum of
the file from which the problem was assembled. It refuses the manifest if the
file at that path has changed. The standalone `write_run_manifest()` collects
the same metadata from a stable file snapshot.

Complete spectra use a binary container because an eigenvector matrix has
quadratic size:

```python
with gliss.Equilibrium("equilibrium_export.nc") as equilibrium:
    with configuration.create_problem(equilibrium) as problem:
        full_result = problem.solve_full_spectrum()
        full_result.write("full-result.gliss")
        problem.write_full_manifest("full-run.gliss", full_result)

restored_full = gliss.FullStabilityResult.read("full-result.gliss")
full_manifest = gliss.FullRunManifest.read("full-run.gliss")
full_manifest.verify_equilibrium("equilibrium_export.nc")
```

Schema `gliss.stability.full-result` stores the certified result
metadata and every eigenvalue, Rayleigh quotient, residual, resolution and
eigenvector for both parity classes. Arrays are uncompressed little-endian
NumPy binary64 records in eigenpair-major order. Reading never enables NumPy
pickle support. Values and signed zeros round-trip exactly, returned arrays
are read-only, and rewriting an unchanged result produces the same container
bytes.

Schema `gliss.stability.full-run` adds the configuration,
equilibrium SHA-256 and software provenance to the same full-spectrum data.
It is self-contained except for the checksummed NetCDF equilibrium. Use
`write_full_run_manifest()` when no assembled problem is available; use
`StabilityProblem.write_full_manifest()` to reject an equilibrium file that
changed after assembly.

The historical controls write schema version 1 byte-for-byte.  Non-default
solver controls write version 2 and are embedded in the configuration and
each certified parity result.  Readers accept both versions and recover the
historical controls when reading version 1.  Full-spectrum readers require the
exact entry set for the declared version, stored without compression or
encryption. They reject invalid entry sets, malformed
metadata, incompatible versions, wrong array types or shapes, inconsistent
diagnostics, and truncated payloads.
Writes stream arrays directly into a temporary container in the destination
directory and replace the destination only after the complete file is closed
and synchronized. A failed write leaves an existing destination unchanged.

Writers use an atomic replacement in the destination directory. Readers
require UTF-8 JSON objects and reject duplicate or unknown fields, missing
fields, incompatible schema versions, nonfinite numbers, inconsistent parity
metadata, invalid vector lengths and truncated files. Errors name the invalid
field and its expected form. `verify_equilibrium()` rejects a file whose size
or SHA-256 differs from the manifest.

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
`gliss_stability_problem_full_spectrum()` accepts separate capacities for its
per-pair arrays and flattened eigenvectors. On a capacity error it reports
both required counts without modifying any data buffer. The header documents
the eigenpair-major layout and dense cost.
Existing ABI-v1 symbols and status values remain unchanged when functions are
added; an incompatible signature, layout or numeric status change requires a
new ABI version.

`gliss_equilibrium_schema_version()` returns 0 for legacy exports and 1 for
GLISS exports. `gliss_equilibrium_write()` writes schema 1 and never replaces
an existing path. Both functions validate null handles, output pointers, path
lengths, embedded null bytes and error-buffer ownership through the same typed
status contract as the other ABI calls.
