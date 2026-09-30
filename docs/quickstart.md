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

```python
modes = [(0, 1), (1, -1), (1, 1)]
with gliss.Equilibrium(stable) as equilibrium:
    with gliss.StabilityProblem(
        equilibrium, modes, degree=1, angular_theta=32, angular_zeta=8
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

## Persistence

Configurations and results round-trip exactly through versioned JSON. A
configuration records the operator revision it was assembled with, so a
record from an earlier discretization is not silently replayed.

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
