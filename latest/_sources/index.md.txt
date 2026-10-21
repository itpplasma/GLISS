# GLISS

GLISS computes ideal-MHD stability of toroidal equilibria with the CAS3D
energy principle on compatible FEEC radial spaces.

```{toctree}
:maxdepth: 2

quickstart
vmec
compatibility
mode_diagnostics
source_inventory
conventions
angular_grids
free_boundary
sparse_stability
api
```

## Scope

Production
: Fixed- and free-boundary FEEC spectra (`StabilityProblem`, with a
  `VacuumModel` for the plasma-vacuum problem and an optional ideal wall),
  the axisymmetric and CAS3D marginality families,
  Mercier profiles, persistence of configurations and results, and the
  conversion of VMEC `wout` and precomputed BOOZ_XFORM `boozmn` files.
  Fixed-boundary pressure-sample spectral derivatives include the pressure
  spline and stiffness response at fixed imported geometry.

Compatibility
: The TERPSICHORE FORT.23/FORT.24 replays reproduce TERPSICHORE's matrices
  for validation. They are not a physical plasma-vacuum API.

Unfinished
: The force-balanced equilibrium-to-spectrum design derivative chain.
  Asymmetric equilibria are
  solved with the coupled parity operator (parity class 0).

Every code block in the quickstart is executed by the native test suite
(`python/tests/test_docs.py`) against the Solov'ev fixtures in `test/data`.

The version menu keeps this active `latest` manual separate from released
documentation. The v0.0.1 and v0.0.2 archives preserve each tag's Python
guide, with its historical APIs and conventions, and identify the exact
source commit in `source.json`. Their version and Mercier examples execute
against a wheel built from that same tag. Published release pages are
immutable; subsequent corrections belong to the active manual.
