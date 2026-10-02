# Quickstart

The examples use the analytic Solov'ev equilibrium (R0 = 1 m, a = 0.33 m,
elongation 1.6) exported by `benchmarks/solovev/public`, with q0 = 1.035
(n = 1 unstable) and q0 = 1.045 (stable). GPEC/DCON puts the marginal point
at q0 = 1.0396. Every block below runs in the test suite, in order.

```python
import os
import tempfile
from pathlib import Path

import numpy as np

import gliss

data = Path(os.environ.get("GLISS_TEST_DATA", "test/data"))
unstable = data / "solovev_q1.035.nc"
stable = data / "solovev_q1.045.nc"
```

## Equilibria

An `Equilibrium` owns the native data. Use it as a context manager, or call
`close()`.

```python
with gliss.Equilibrium(unstable) as equilibrium:
    print(equilibrium.schema_version, equilibrium.coordinate_handedness)
    assert equilibrium.coordinate_handedness == "left-handed"
```

## Mercier stability

Positive `D_Mercier` is stable, as in VMEC. `mercier_objective` returns
`-min(D_Mercier)`, positive exactly when a surface is Mercier-unstable.

```python
s, d_mercier = gliss.mercier_profile(unstable)
assert d_mercier[0] < 0.0 and np.all(d_mercier[s > 0.25] > 0.0)
assert gliss.mercier_objective(unstable) > 0.0
```

## The n = 1 axisymmetric family

`solve_axisymmetric` counts negative directions of the two-component
marginality operator and certifies the lowest eigenpair. The trial space is
conforming, so its eigenvalues approach the converged ones from above; on
these 16-surface exports degree 3 resolves the near-marginal kink.

```python
for path, expected in ((unstable, 1), (stable, 0)):
    with gliss.Equilibrium(path) as equilibrium:
        result = gliss.solve_axisymmetric(equilibrium, poloidal_max=6, degree=3)
    assert result.negative_count == expected
    assert result.certificate < abs(result.lowest_eigenvalue)
```

## Fixed-boundary spectra

A `StabilityProblem` assembles the compressible fixed-boundary operator with
physical mass for both stellarator-symmetry parity classes. Eigenvalues are
omega^2 in s^-2.

The magnetic differential equation projects cubic position-harmonic products.
Its angular grid must exceed six times the largest absolute stored equilibrium
harmonic index in each direction. These Solov'ev exports contain poloidal
harmonics through m = 8, so the examples use 64 poloidal points; their
axisymmetric position has n = 0. The vacuum edge mesh has its own resolution
and is independent of this volume projection grid.

```python
modes = [(0, 1), (1, -1), (1, 1)]
with gliss.Equilibrium(stable) as equilibrium:
    with gliss.StabilityProblem(
        equilibrium, modes, degree=1, angular_theta=64, angular_zeta=8
    ) as problem:
        result = problem.solve()
        lowest = result.classes[0]
        energy = problem.energy(1, lowest.eigenvector)
        configuration = problem.configuration
assert lowest.negative_count == 0
assert abs(energy.rayleigh_quotient - lowest.lowest_eigenvalue) <= (
    lowest.certificate + 1e-9 * abs(lowest.lowest_eigenvalue)
)
```

## Pressure-sample derivatives

The pressure JVP follows the imported sample spline through the assembled
stiffness at fixed geometry. For an isolated eigenvalue, its contraction uses
the mass-normalized eigenvector; a complete cluster uses its basis-invariant
trace. The positive `gap` specifies an exterior separation in s^-2.

```python
with gliss.Equilibrium(stable) as equilibrium:
    with gliss.StabilityProblem(
        equilibrium, [(1, 1)], degree=1, angular_theta=64,
        angular_zeta=16, radial_cells=2,
    ) as problem:
        nodes, pressure_pa = problem.pressure_samples()
        spectrum = problem.solve_full_spectrum_class(1)
        stop = spectrum.eigenvalues.size
        direction_pa = 0.01 * pressure_pa
        jvp = problem.spectral_pressure_jvp(
            1, stop - 1, stop, direction_pa, gap=1e-4,
        )
        sensitivity = problem.spectral_pressure_sensitivity(
            1, stop - 1, stop, gap=1e-4,
        )
assert np.isclose(sensitivity.jvp(direction_pa), jvp, rtol=1e-9)
assert not sensitivity.gradient.flags.writeable
```

The gradient has units s^-2 Pa^-1 and remains usable after the problem closes.
Constructing it costs one tangent assembly per pressure sample; a direct JVP
uses one. Pressure must be strictly positive, including at the assembly
points, and zero-width resonances must admit the full sample-direction domain.
Geometry, magnetic profiles, density, gamma, modes, quadrature and resonance
topology are held fixed. A force-balanced external equilibrium response and
free-boundary pressure derivatives remain future work.

## Free-boundary spectra

A `VacuumModel` frees the plasma edge: its normal displacement drives a
current-free vacuum field, whose energy (a scalar-potential boundary integral
over the edge and wall) enters the stiffness. The model sets the full-torus edge mesh and an optional
ideal wall, here conformal at 3 cm. Without a wall this q0 = 1.045 Solov'ev
edge is kink unstable; the close wall stabilizes it.

```python
free_modes = [(0, 1), (1, -1), (1, 1), (2, -1), (2, 1)]
lowest_free = {}
with gliss.Equilibrium(stable) as equilibrium:
    for wall in (None, 0.03):
        with gliss.StabilityProblem(
            equilibrium, free_modes, degree=2, angular_theta=64,
            angular_zeta=8, vacuum=gliss.VacuumModel((24, 12), wall),
        ) as problem:
            lowest_free[wall] = problem.solve_class(1)
assert lowest_free[None].negative_count >= 1
assert lowest_free[0.03].negative_count == 0
assert lowest_free[None].boundary_condition == "free"
```

## Persistence

Configurations and results round-trip exactly through versioned JSON. A
configuration records the operator revision it was assembled with, so a
record from an earlier discretization is not silently replayed.
Fresh spectra also record the source equilibrium and canonical configuration
digests, including the radial mesh and exact wall model. Run-manifest writers
require the current operator revision, 4, and reject mismatched inputs.
Earlier revisions remain readable as
historical records; migrate the configuration explicitly with
`dataclasses.replace(configuration, discretization_revision=4)` and solve again.
An old spectrum cannot be relabelled or exported as a verified current run.
The Python loader checks both ABI 4 and native operator revision 4.
Historical manifests preserve their fingerprints. An old operator revision
reports `configuration_verified=False`; `equilibrium_verified` still reports
whether its recorded source digest is known. Archives without either digest
report both flags as false. Creating a new manifest requires a fresh solve.

```python
with tempfile.TemporaryDirectory() as directory:
    configuration.write(Path(directory) / "configuration.json")
    result.write(Path(directory) / "result.json")
    restored = gliss.StabilityConfiguration.read(
        Path(directory) / "configuration.json"
    )
    reread = gliss.StabilityResult.read(Path(directory) / "result.json")
assert restored == configuration
assert reread.classes[0].lowest_eigenvalue == lowest.lowest_eigenvalue
```
