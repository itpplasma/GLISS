# Mode diagnostics

A cross-code mode comparison needs both a coordinate map and a shared physical
inner product. The `gliss_compatible_marginality` diagnostic uses the current
compatible FEEC assembly to evaluate coefficient vectors and invariant
subspaces after that map has been established. Its CSV input represents GLISS
coefficients; it does not interpret an external code's variable layout.

The public production spectra remain on the block solver. Fixed-shift counts,
negative brackets and signed brackets in this diagnostic also assemble block
storage directly. External vector energies, seeded subspaces and radial profiles
use bounded dense matrices because they retain decomposed potential terms.
The conservative input limit is 2048 unknowns; larger diagnostic requests fail
before matrix assembly.

## Coefficients and axis regularity

For each active Fourier trial, normal coefficients use the retained H1 basis
and tangential coefficients use the L2 basis. The full ordering is H1 basis
index first, then active trial, followed by L2 basis index and active trial.
The current Cartesian axis constraint removes the first eta coefficient for
each regular |m|=1 trial. Input vectors use the resulting reduced ordering;
profile reconstruction restores each eliminated coefficient through the same
axis map used in assembly.

The perpendicular kinetic mass uses unit density unless `--physical-density`
sets another positive SI density. The default eta weights equal the normal stored-power table, preserving the
production regularity convention. An explicit different `--eta-stored-powers`
table selects the historical unconstrained validation space. Its results must
not be presented as production stability evidence.

`--evaluate-vector=FILE` expects consecutive `index,value` rows starting at one.
`--evaluate-subspace=FILE` expects a list of vector-file paths; relative paths
are resolved from the list's directory. The parser checks each listed vector
before opening the equilibrium. The seeded solve additionally requires both
`--subspace-shift` and `--subspace-iterations`. A nonpositive mass or exhausted
iteration budget is a failure; a returned Ritz vector alone does not certify
that the requested physical branch was selected.

## Numerical evidence

The vector reconstruction test uses explicit degree-one shape functions and
an independently prescribed axis relation. The seeded-solver tests use exact
diagonal pencils, an invariant seed with positive Gram matrix in an indefinite
full mass, a deliberately insufficient iteration budget, and rotated/scaled
bases spanning the same analytic low subspace. Cancellation tests use an exact
quadratic form whose terms cancel beyond ordinary summation precision.

Mass-whitened residuals estimate distance to the finite-pencil spectrum in
exact arithmetic. Floating-point diagnostics, inertia probes and subspace
angles require conditioning assessment. They do not certify equilibrium
quality, discretization convergence or agreement between unlike mass matrices.

## Integrated branches

The integration starts from `d0d1cfd988acb06ea21f5fb960d4b6b700c31127` and
includes `validation/mishka-low-subspace` at
`143dfc5`, which contains the complete `axis-regularity-trace` and
`validation/mishka-mode-energy` histories.

| Branch changes | Integration decision |
| --- | --- |
| Stored-power controls and eta weights | Preserve controls with current regularity defaults |
| Negative/signed brackets and fixed shifts | Use current sparse FEEC assembly and block inertia |
| Leading-cell trace selection | Add `--first-cells` beside existing explicit cell selection |
| Eigenprofiles | Reconstruct the current reduced axis space |
| Cancellation-aware energy closure | Retain compensated sums and forward-error diagnostics |
| External coefficient vectors and subspaces | Retain bounded validation diagnostics and path preflight |
| Legacy parity tests | Current coupled-parity and pairing tests cover the replaced backend |

## MISHKA/CASTOR acceptance

Issue #12 remains open until a convention-complete external adapter and
independent low- and high-beta fixtures are available. Each comparison must
record its coordinate, toroidal-sign, Fourier-layout and normalization
fingerprints; map the compatible constraints; and preserve mass, total
potential and decomposed energy within frozen uncertainty bounds. The same
invariant branch must survive at least three radial resolutions. Corrupted
layout/sign controls must fail, and continuum branches must be distinguished
from isolated or clustered eigenmodes.
