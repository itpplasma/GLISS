# GLISS

GLISS computes ideal-MHD stability of toroidal equilibria with the CAS3D
energy principle on compatible FEEC radial spaces.

```{toctree}
:maxdepth: 2

quickstart
vmec
compatibility
conventions
api
```

## Scope

Production
: Fixed-boundary FEEC spectra of stellarator-symmetric equilibria
  (`StabilityProblem`), the axisymmetric and CAS3D marginality families,
  Mercier profiles, persistence of configurations and results, and the
  conversion of VMEC `wout` and precomputed BOOZ_XFORM `boozmn` files.

Compatibility
: The TERPSICHORE FORT.23/FORT.24 replays reproduce TERPSICHORE's matrices
  for validation. They are not a physical plasma-vacuum API.

Unfinished
: The free-boundary plasma-vacuum solve, the coupled operator for
  asymmetric equilibria, and the equilibrium-to-spectrum derivative chain.
  Asymmetric files convert and read, but the stability operator refuses
  them.

Every code block in the quickstart is executed by the native test suite
(`python/tests/test_docs.py`) against the Solov'ev fixtures in `test/data`.
