# GLISS verification and development plan

Updated September 9, 2026 after an independent Astra xhigh source and benchmark
audit. GLISS is a research-grade fixed-boundary ideal-MHD value evaluator.
The complete differentiable equilibrium-to-spectrum optimization workflow is
unfinished. A certificate for a discrete matrix eigenpair does not certify
equilibrium quality, discretization convergence, or agreement with another code.

## Audit baseline and evidence

The reviewed GLISS base is
[`b994f51aa9ebe1a5e03bad2c72c6fde730fcaa1b`](https://github.com/itpplasma/GLISS/commit/b994f51aa9ebe1a5e03bad2c72c6fde730fcaa1b)
(release 0.0.2), initially clean. The empty worktree patch has SHA-256
`e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`.
The research evidence was inspected at
[`gvec-stability` commit `822b2bb0296f48cedd246a67cc46ad69266a5cb7`](https://gitlab.tugraz.at/plasma/proj/stel/gvec-stability/-/tree/822b2bb0296f48cedd246a67cc46ad69266a5cb7).
That repository requires separate access; its data are not bundled with GLISS.

Read [PROVENANCE.md](PROVENANCE.md) for the formulation map and the research
repository's `ROADMAP.md`, `validation/pins.env`, acceptance manifests, and
benchmark READMEs before extending a physics claim. Preserve their distinct
solver pins: most archived numbers were produced before the released base.

The audit covers source, Python/C interfaces, tests, derivative gates, CI,
analytical evidence, and the archived cross-code comparisons. It does not
constitute an independent rederivation of every equation or a fresh execution
of every external solver.

## Corrections from this audit

- Scale symmetric 2-by-2 pivots before calculating determinant signs for
  inertia. The previous calculation underflowed or overflowed for finite
  matrices. The regression uses the exact eigenvalues `-s, 3s` of
  `s*[[1,2],[2,1]]` for `s=1e-200, 1, 1e200`; it failed before the fix.
- Apply the existing marginality Fourier bandwidth admission rule to the
  public compressible fixed-boundary constructor. Its 64-point grid previously
  accepted Nyquist modes and coincident sampled harmonics. The independent
  oracle is Fourier orthogonality: modes 31 and 33 have zero continuum cross
  mass, but the 64-point sum gives 0.5; the resolved 128-point sum gives zero.
  This conservative guard does not bound nonlinear geometry quadrature error.
- Reject nonfinite Enzyme gradients and finite-difference references before
  comparing errors. Previously a NaN could pass the inequality test. Preserve
  this requirement for every new derivative gate, including native and Python
  wrappers; the local eigenvalue-sensitivity test now checks finiteness too.
- Handle singular inertia probes when a bracket endpoint or midpoint equals
  an exact eigenvalue. Dense refinement and the positive-spectrum path must
  choose a successfully factorized shift. The independent regressions use
  exact diagonal pencils, repeated eigenvalues, zero modes, non-unit positive
  masses, and a rejected indefinite-mass control.
- Verify inverse-density scaling through the production constructor and solve:
  quadrupling density divides `omega^2` by four and preserves negative inertia.

## Priority 1: qualify the production operator

Related issues: [#13, higher-order FEEC certification](https://github.com/itpplasma/GLISS/issues/13),
[#10, asymmetric input](https://github.com/itpplasma/GLISS/issues/10).

- [ ] Complete production analytical coverage and public-API qualification.
  `benchmarks/analytic/run.sh` now runs an exact straight-cylinder case through
  current FEEC surface assembly for degrees 1–4, three meshes and both parities,
  with independent Bessel/dispersion references. Fast modes show second/fourth
  order for degrees 1/2; higher degrees reach numerical precision limits.
  The legacy cylinder fixture reconstructs as a double-covered torus in the
  primitive path and is retained only as an explicitly invalid comparison.
  The lower-level exact-cylinder runner does not test import or public admission.
  Newcomb/Suydam and complete public-API analytical coverage remain unfinished.
- [ ] Cover homogeneous Alfvén and magnetosonic branches, theta-pinch fast and
  slow branches, screw-pinch/Newcomb and Suydam thresholds, and the helical
  cylinder vertical threshold. Require analytical eigenvalues or independently
  derived marginal limits, both stable and unstable controls, density/field/length
  scaling, and units. Distinguish physical mass from artificial coefficient norms.
- [ ] For FEEC degrees 1 through 4, measure radial convergence on at least three
  meshes, verify commuting derivatives and traces, and test both parity classes.
  Use smooth manufactured solutions for optimal rates and state the separate
  behavior of singular or continuum solutions.
  The new cylinder slow-cluster diagnostic demonstrates ordering-dependent
  dense-solver roundoff: 50-digit solves of the same matrices recover the
  analytical cusp limit while a double-precision triangle/order choice drifts.
  `benchmarks/analytic/conditioning.py` reproduces this distinction. Small
  globally scaled residuals do not establish relative slow-mode accuracy.
- [x] Expose angular resolution through configuration, the C ABI, persistence,
  and Python. Defaults remain 64 by 64; C has an additive v2 constructor and
  configuration schema 4 migrates older files. Nonlinear metric convergence
  remains a separate qualification requirement from Fourier admission.
  Freeze mode topology and quadrature during differentiation.
- [x] Add a frame-aware physical symmetry admission check before solving separate
  parity classes. A full-storage export can be physically symmetric even when
  `stellarator_symmetry=False`; rejecting that flag alone is wrong. Conversely,
  genuine asymmetry couples parity classes and is not covered by today's solves.
  Verify symmetric full-storage files and an odd-harmonic cross-parity mass
  oracle. Coupled asymmetric assembly remains #10.
  The admitted operator is sampled at assembly nodes with dimensionally scaled
  tolerances; this does not certify continuum symmetry or angular convergence.
- [x] Reject incompatible nonconstant magnetic differential-equation resonances
  in `src/export_surface_geometry.f90`. The retained mean-projected equation now
  uses a RHS/projection roundoff bound and rejects unresolved nonzero forcing.
  Resolved denominators retain their exact inverse and matching flux derivatives.
  Independent resonance, near-resonance, scaling, finite-difference and injected
  NaN controls pass. Mercier failures propagate and invalidate legacy outputs.
  Mean force balance and Fourier truncation remain separate diagnostics.

Completion requires frozen acceptance bounds, independent reference formulas,
negative controls, and a clean-checkout runner. Do not tune a tolerance to a
new GLISS result or promote a stable sign from an unresolved trial space.

## Priority 2: establish production cross-code agreement

Related issues: [#12, MISHKA/CASTOR mode transfer](https://github.com/itpplasma/GLISS/issues/12),
[#11, sparse CAS3D L139 solve](https://github.com/itpplasma/GLISS/issues/11).

| Comparison | Evidence at the audited research pin | Remaining acceptance requirement |
| --- | --- | --- |
| TERPSICHORE FORT.23/24 | Useful same-discretization compatibility replay | Run the independent production FEEC operator on the same qualified equilibrium; match boundary condition, normalization, modes, and energy terms |
| QAS3 production FEEC | The 191-mode deck supplies a mode mask; ns64-to-ns128 lowest-eigenvalue drift is about 5.53%, with material force-balance residuals | Converge equilibrium and FEEC errors separately before claiming same-physics agreement |
| W7-X / Nuehrenberg 1996 | The documented coefficient-normalized L10 result is about -0.88701 versus digitized -0.37148 in the report's scaled units | Resolve normalization, radial form functions, reference length, and unavailable deck details; finish quotient-aware L139 scaling under #11 |
| MISHKA / CASTOR | Branch transfer is unresolved; the CASTOR low-beta stable control currently fails | Transfer compatible invariant subspaces across at least three meshes, distinguish continuum branches, compare mass and decomposed potential energy under #12 |
| DCON / Solov'ev | A pinned axisymmetric comparison exists | Reproduce it with current unrestricted 3-D assembly and both stable/unstable controls; retain convention and boundary-condition evidence |
| Moderate figure-8 | Research roadmap specifies a common VMEC reference, then GVEC reproduction | Qualify one canonical input, reproduce surfaces and profiles across representations, then compare converged stability and modes |

- [ ] Regenerate each retained comparison at an exact current GLISS commit;
  record external-code commit, input hashes, compiler/libraries, coordinates,
  normalization, uncertainty, and accepted resolution range.
- [ ] Compare sign and inertia first, then normalized eigenvalues, invariant
  subspaces, and energy decomposition. Individual eigenvectors inside a cluster
  and raw eigenvalues under unlike mass matrices are not valid comparisons.
- [ ] Retain failed controls and unexplained discrepancies in the report.
  Replaying another code's stored matrix does not validate the production
  physical plasma-vacuum solver.
- [ ] Obtain any missing original CAS3D decks rather than fitting undocumented
  settings. Helios data are request-only and not yet an available benchmark;
  linear ideal results must also be distinguished from nonlinear M3D-C1 evolution.

## Priority 3: complete the differentiable workflow

Tracked by [#9, equilibrium-to-spectrum and clustered-subspace derivatives](https://github.com/itpplasma/GLISS/issues/9).

Existing public Rayleigh JVP/VJP actions differentiate the displacement vector
with the assembled matrices held fixed. The eight Enzyme gates exercise local
kernels and compatibility maps. The SIMSOPT Mercier wrapper is value-only.
The public `spectral_parameter_sensitivity` now supplies isolated and cluster
trace derivatives for density and adiabatic index with the imported equilibrium
held fixed. Analytical contractions pass independent pencil, basis-rotation,
duality and production finite-difference plateau checks. Its residual-based
gap admission is a numerical diagnostic, not a rigorous spectral enclosure.
None of this establishes a GVEC-design-variable gradient. Present-tense claims
to that effect in the research
`docs/sections/differentiation.tex` need correction.

- [ ] Specify parameter ownership, units, constraints, fixed topology, and
  supported coordinate transformations. Define equilibrium sensitivities through
  force balance and the G-frame/Boozer map; do not silently treat external
  VMEC/BOOZ preprocessing as differentiable.
- [ ] Connect geometry/profile JVPs and VJPs to stiffness and physical mass
  assembly, then isolated generalized-eigenvalue sensitivities.
- [ ] Implement basis-invariant cluster objectives with a declared spectral
  gap and cluster-selection policy. Reject unresolved crossings and changing
  cluster dimensions; an ordered minimum need not be differentiable there.
- [ ] Require JVP/VJP dot-product duality, independent directional finite
  differences over a step-size plateau, refinement studies, and NaN/error
  controls at every layer and across complete GVEC-to-GLISS evaluations.
- [ ] Test any nonsmooth guards, resonant denominators, eigenvalue selection,
  and failed solver certificates explicitly. Admission checks and integer
  inertia are validity diagnostics, not differentiable optimization objectives.
- [ ] Run pinned derivative gates on changes to covered kernels and before
  promotion. The current GitHub Enzyme job is manual-only.
- [ ] Expose real optimization DOFs in SIMSOPT and verify improvement with a
  fresh converged evaluation and an independent derivative-free control.

## Priority 4: make evidence reproducible in CI

Related issue: [#15, executable versioned documentation](https://github.com/itpplasma/GLISS/issues/15).

- [ ] Publish the minimal redistributable fixtures and acceptance manifests
  needed for current production analytical and cross-code gates. Link restricted
  research evidence explicitly and report tests skipped for missing data.
- [x] Test an installed wheel with the physics suite, not just import/version
  smoke checks. `ci/check_installed.py` enforces bundled-library loading, checks
  production inverse-density scaling, and runs the portable Python suite from
  a temporary directory outside the checkout. The cylinder is a synthetic ABI
  control, not a qualified straight-cylinder analytical benchmark.
- [ ] Add controlled thread counts, pinned toolchain identities, and derivative
  gate artifacts to CI. Keep performance measurements separate from correctness.
  CI now fixes BLAS/OpenMP thread counts to one; automatic derivative gates and
  complete toolchain provenance remain unfinished.
- [ ] Run bounds/runtime checks and a second compiler, then the optimized
  array-temporary audit before release. Record failures and compiler warnings.
- [ ] Keep equation-to-source references and benchmark hashes current; distinguish
  analytical proof, independent solver comparison, compatibility replay, and
  deterministic regression in every reported result.

## Priority 5: extend the supported scope

- [ ] [#8: physical free-boundary plasma-vacuum Python solve](https://github.com/itpplasma/GLISS/issues/8):
  verify wall-distance/angular convergence, energy balance, analytical or
  independent-code controls, and derivatives before optimization use.
- [ ] [#10: asymmetric VMEC and precomputed BOOZ_XFORM](https://github.com/itpplasma/GLISS/issues/10):
  implement full parity coupling and convention-complete round trips after the
  immediate input-admission work in Priority 1.
- [ ] [#14: macOS wheels](https://github.com/itpplasma/GLISS/issues/14) and
  [#15: documentation site](https://github.com/itpplasma/GLISS/issues/15):
  complete clean-install physics checks on each advertised platform.

The GitHub audit found open issues #8 through #15 and closed
[#1](https://github.com/itpplasma/GLISS/issues/1), the earlier family/single-mode
assembly discrepancy. No issue was closed merely because this audit passed.

## Verification record

At the original baseline, controlled-thread `fo` passed build, tests, static
checks, and lint. The default unrestricted-thread attempt timed out; this is
not evidence of a physics failure. Python with the VMEC extra passed 237 tests
and skipped two optional cases. All eight original and hardened Enzyme gates
passed with Flang/LLVM 22 and LLVMEnzyme-22. A deliberately injected NaN
gradient failed the hardened drive gate as required.

The integrated corrections pass all 101 GNU Fortran CTest entries and all 109
Flang/Enzyme entries, with clean lint and changed-file formatting. An installed
wheel tested from a separate directory, explicitly loading its bundled native
library, passes 234 Python tests and skips SIMSOPT and the optional external
Mercier golden fixture. Three repository-registry checks require source files
and were excluded from that wheel run; they are structural checks, not physics
evidence. The source-tree Python run reports 237 passed and the same two skips.
The integrated GNU Fortran runtime-checking build also passes all stages with
`-fcheck=all -O1`; an earlier unoptimized runtime-checking snapshot passed too.

Before committing corrections, run the full `fo` pipeline, the native Python
suite, focused independent regressions, and the changed Enzyme gates. Freeze
uncommitted review inputs using the base commit and patch digest. Keep the
controller responsible for integration and promotion; reviewers return evidence.
