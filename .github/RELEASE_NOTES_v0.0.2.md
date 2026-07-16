# GLISS 0.0.2

GLISS 0.0.2 is a fixed-boundary physics and Python-package release.

## Supported scope

- Python-first `Equilibrium`, `StabilityConfiguration`, and
  `StabilityProblem` contexts
- certified lowest eigenpairs and opt-in complete spectra
- five-term energy decomposition and Mercier diagnostics
- deterministic result, run, and equilibrium persistence
- symmetric GVEC and direct symmetric VMEC conversion
- compatible FEEC degrees 1 through 4
- explicitly labeled fixed-boundary and pressureless-pseudoplasma
  TERPSICHORE compatibility replays

The validation set covers homogeneous and cylindrical manufactured cases,
Newcomb/Suydam behavior, the Solov'ev/DCON marginal sequence, a matched MISHKA
case, QAS3 native/replay/FEEC comparisons, and the W7-X CAS3D Figure 5 trend.
The normalized W7-X TERPSICHORE/FEEC shape difference remains below the frozen
`0.025` acceptance bound.

## Limits

- Production support is fixed-boundary. Selected free-boundary operators are
  not a public physical plasma-vacuum API.
- TERPSICHORE replay uses that code's stored discretization and normalization;
  it is a compatibility path, not an independent physical method.
- The equilibrium-to-spectrum and clustered-subspace derivative chain is not
  complete.
- VMEC conversion requires stellarator symmetry. Precomputed BOOZ_XFORM files
  and asymmetric equilibria are not supported.
- The unpublished W7-X CAS3D angular grid, form functions, and reference length
  prevent an exact same-deck reproduction. The comparison retains the audited
  digitization uncertainty and does not equate unlike positive norms.
- TERPSICHORE `LCURRF=9` is an unresolved sensitivity and is not a validated
  GLISS mode. The W7-X reference uses `LCURRF=1, DELTAJP=0.04`.
- This release provides a manylinux x86-64 wheel and source distribution.
  macOS wheels remain future work.

## Installation

```sh
python -m pip install gliss==0.0.2
```

VMEC conversion additionally requires:

```sh
python -m pip install "gliss[vmec]==0.0.2"
```
