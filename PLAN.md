# GLISS verification and development plan

Updated October 2, 2026 during integration, independent regression review,
installed-package validation and fresh reference measurements. Historical
September 30 evidence includes public-source DCON/GVEC
reproduction of the Solov'ev benchmark (see the corrections section below);
previously updated September 9, 2026 after an independent Astra xhigh source and
benchmark audit, parallel fixes, and fresh analytical and Solov'ev runs.
GLISS is a research-grade ideal-MHD value evaluator. Its physical
free-boundary implementation still lacks the complete qualified acceptance case.
The complete differentiable equilibrium-to-spectrum optimization workflow is
unfinished. A certificate for a discrete matrix eigenpair does not certify
equilibrium quality, discretization convergence, or agreement with another code.

## Delivery status on October 2, 2026

Development stopped at the user's request. The completed integration is being
committed and promoted to `main`; further scientific work and issue closure
are deferred. The production source and native test snapshot is
[`6a27d11`](https://github.com/itpplasma/GLISS/commit/6a27d113d1cca29f641b68c83e2406a1b1ce1c68).
The subsequent
[`818826b`](https://github.com/itpplasma/GLISS/commit/818826b92e9a5b3150984526233f5b5b5d3fa451)
adds the bounded source inventory. This status update changes documentation.
The development package is 0.0.3, with native ABI 4 and discretization revision
4. Existing released tags 0.0.1 and 0.0.2 retain their original contents.

Completed work includes integration of the original open branch/PR history,
axis interpolation and regularity repairs, robust inertia and initial-rank
checks, exact-resonance derivative and zero-iota gauge repairs, pressure
JVP/VJP interfaces at fixed geometry, coupled asymmetric-wall support,
VMEC/BOOZ metadata admission, operator provenance and native-handle ownership.
Build, installed-wheel, documentation-example and strict macOS CI checks are
implemented. These bounded corrections do not complete all scientific issues.

At the source snapshot, GNU `fo` passes all stages and 175 native tests; the
Flang/LLVM 22 pipeline passes 183 tests, including eight Enzyme gates. The
optimized `-O3 -Warray-temporaries -Werror=array-temporaries` audit passes
175 tests. Source Python reports 329 passed and two optional skips. The
installed Linux wheel reports 324 passed and two optional skips, including
35 tests marked `native`; all ten quickstart/VMEC documentation blocks execute.
Sphinx HTML, link checking with warnings as errors, and `pip check` pass.
The skips require optional SIMSOPT or externally supplied Mercier golden data.
The installed wheel SHA-256 is
`beb4904769bf7da4794a09d560ffe5516734097f4a28b948a41de7ea5ce5363f`.

[Hosted CI run 37035917702](https://github.com/itpplasma/GLISS/actions/runs/37035917702)
passes at that exact source commit.
[macOS run 37035917761](https://github.com/itpplasma/GLISS/actions/runs/37035917761)
passes all 175 native tests on each of arm64 and x86-64, all twelve installed
CPython 3.9–3.14 wheel jobs, and the public C consumers. Each macOS Python
suite reports 313 passed and four optional skips. This strict-shell run
supersedes earlier runs that masked a failing native test; those earlier green
statuses are not qualification evidence. The final documentation handoff also
passes the full local `fo` pipeline (58.2 seconds); later hosted runs are separate.

Fresh analytical and DCON sign comparisons, input hashes, solver certificates,
timing/RSS observations, and failed free-boundary qualification controls are
retained in the [October 2 measurement record](benchmarks/results/2026-10-02/README.md).
They establish bounded analytical accuracy and stable/unstable sign agreement,
not a converged marginality interval or superiority over every competitor.

## Remaining acceptance work

The issue audit identified bounded corrections ready for closure review in
#16, #17, #20, #21, #23–#29, #32 and #34. The final macOS run also supplies the
remaining execution evidence for #14. These issues are not automatically
closed by this handoff. The substantive remaining requirements are:

- #8: qualify true-edge flux derivatives and equilibrium exterior matching,
  then run installed plasma-vacuum wall/angular refinement acceptance. The
  corrected M32 geometry passes, but native edge `chi'/Phi'` error is
  `3.1952e-7` against `1e-12`; flux-slope error is `1.1983e-14` against `1e-14`.
- #9: complete geometry/conversion, force-balance, mass and vacuum parameter
  derivatives and the force-balanced optimization chain with independent
  finite-difference plateaus and transpose checks.
- #10 and #22: add unsupported precomputed BOOZ-file version and corrupted
  metadata controls, and the requested installed LI383 conversion/native
  truncated-geometry acceptance case.
- #11 and #12: obtain the original CAS3D W7-X deck and normalization, and
  implement convention-complete MISHKA/CASTOR mode transfer with independent
  branch and invariant-subspace acceptance across radial resolutions.
- #13: certify the full higher-order manufactured FEEC gradient/curl sequence,
  dense element oracles and three-dimensional energy/spectrum convergence.
  Radial sequence controls and the passing optimized audit cover only parts.
- #15: verify an actual deployed documentation site and immutable pages for
  released tags. Installed examples and Sphinx checks pass; hosted publication
  and released-tag pages have not been qualified at this stopping point.
- #18 and #19: add the requested automated VMEC Mercier-sign control and the
  complete primitive exact-torus off-knot Jacobian sweep at M16/24/28/36.
- #30: file the retained GVEC export reproductions upstream when tracker access
  is available; the attempted tracker endpoint returned 404.
- #33: demonstrate resolved QAS/LI383 rational-surface convergence and expose
  affected-surface status. Finite-cell inverse regularization is a model choice.
- #35: meet the degree-1 q0=1.035 instability criterion by 64 radial cells and
  repeat the higher-degree sweep on qualified inputs.

The source inventory accounts for 311 files, including 148 required production
paths, 37 C prototypes and 50 file records carrying bounded evidence. It records zero
whole-component equation certifications. A complete independent source audit,
the qualified figure-eight equilibrium and common cross-code study remain open.

An unfinished released-documentation bootstrap helper was stopped before
execution or integration. It remains in the detached local worktree
`/home/ert/code/GLISS-release-docs-bootstrap` at `818826b`. A backup is
`/home/ert/code/GLISS-release-docs-bootstrap-818826b.tar.gz`, SHA-256
`a332a100d22935c1f9429475b3a133bd23511dca4f5f37461faf4a18b5389ab2`.
It has no test evidence and is excluded from the production integration.

## Audit baseline and evidence

The earlier reviewed GLISS base was
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

## Corrections from the September 30, 2026 audit

- The DCON Solov'ev discrepancy was a GLISS discretization error
  ([#16](https://github.com/itpplasma/GLISS/issues/16)). A regular displacement
  has `eta ~ s^(-1/2)` for |m|=1, but eta used an unweighted L2 basis; the
  compression term then diverged like `1/s`, hidden by reduced quadrature, and
  a spurious mode with eigenvalue proportional to the mesh width gave count one
  at every resolution. eta now carries the normal-component axis factor and the
  third unknown is the regular `nu = mu - (FP'/FT') sqrt(g) eta`. GLISS now
  reproduces the DCON bracket (1.039062 unstable, 1.039843 stable). A fresh
  GPEC `f5595c06` DCON build confirmed the archived bracket; its converged
  marginal q0 is 1.03956-1.03959.
- Positive-spectrum solves used the midpoint of an unrefined inertia bracket
  and could return a non-lowest eigenvalue
  ([#17](https://github.com/itpplasma/GLISS/issues/17)); both paths now bisect.
- `benchmarks/solovev/public` regenerates the benchmark from public GVEC 1.5.0
  with a patch for five export defects, tracked for upstream reporting in
  [#30](https://github.com/itpplasma/GLISS/issues/30).
- Fixed on `claude/zen-babbage-dgxz45`, one commit per issue:
  [#18](https://github.com/itpplasma/GLISS/issues/18) Python Mercier sign,
  [#19](https://github.com/itpplasma/GLISS/issues/19) axis-regular spline
  roundoff, [#20](https://github.com/itpplasma/GLISS/issues/20) truncated
  Pfirsch-Schlueter solve, [#21](https://github.com/itpplasma/GLISS/issues/21)
  handedness metadata, [#22](https://github.com/itpplasma/GLISS/issues/22) VMEC
  truncation checks, [#23](https://github.com/itpplasma/GLISS/issues/23) block
  inertia reliability, [#24](https://github.com/itpplasma/GLISS/issues/24)
  full-spectrum cost, [#25](https://github.com/itpplasma/GLISS/issues/25) path
  truncation, [#26](https://github.com/itpplasma/GLISS/issues/26) certificate
  semantics, [#27](https://github.com/itpplasma/GLISS/issues/27) scale-relative
  tolerances, [#28](https://github.com/itpplasma/GLISS/issues/28) native Python
  evidence, [#29](https://github.com/itpplasma/GLISS/issues/29) operator
  revision and handle finalizers,
  [#32](https://github.com/itpplasma/GLISS/issues/32) TERPSICHORE replay of
  stable and IVAC>0 files, [#33](https://github.com/itpplasma/GLISS/issues/33)
  resonant Pfirsch-Schlueter harmonics, and
  [#34](https://github.com/itpplasma/GLISS/issues/34) nonpositive truncated
  metric.
- Roadmap work on this branch: precomputed `boozmn` and asymmetric VMEC
  conversion (#10 parts A and B) and the coupled parity operator for
  equilibria without stellarator symmetry (#10 part C, checked against the
  class union of a shifted angle origin), radial FEEC convergence controls
  (complete certification remains #13), the native test harness without path
  hacks and macOS wheel workflow (#14), and documentation sources whose
  quickstart runs against the installed wheel. Site deployment is verified
  separately for #15.
- [#35](https://github.com/itpplasma/GLISS/issues/35) conforming axis space:
  the leading |m|=1 coefficients of xi^s and eta are tied as a smooth
  displacement requires, the axis element is integrated in sqrt(s), and the
  marginality mass is the perpendicular kinetic form (the coefficient norm
  was unbounded for regular |m|=1 fields). Degree 3 now converges at order
  3.3 on the Solov'ev sweep; degree 1 was observed to
  converge from above at order 1.3-1.4 and needs ns > 128 for q0 = 1.035.
  This observation does not prove upper bounds for the general discretization.
- Upstream trackers stay open until the upstream reports are filed:
  [#30](https://github.com/itpplasma/GLISS/issues/30) GVEC CAS3D export and
  [#31](https://github.com/itpplasma/GLISS/issues/31) DCON default tolerances.
  The DCON report is now filed as
  [GPEC #292](https://github.com/PrincetonUniversity/GPEC/issues/292), and #31
  is closed. GVEC's upstream tracker remains inaccessible to this session;
  #30 is open.
- `benchmarks/qas3` compares STELLOPT VMEC, TERPSICHORE 1.2 and GLISS on the
  QAS3 family. The TERPSICHORE replay agrees to 1e-8. The independent
  critical current agrees once resolved: on an ns = 257 equilibrium both
  codes are unstable at 0.75 of the base current, where GLISS on 64 cells
  of the ns = 65 equilibrium was stable. The li383 D_curr difference was
  VMEC's angular truncation of J.B: raising VMEC mpol/ntor from 9/5 to 16/10
  brings its median D_curr to 0.999 of GLISS's. The `convert_vmec` force-balance and |B|
  gates were resolved: the pointwise metric-B_s closure (0.1-0.26 at any M,
  N, ns and Boozer resolution, while the formula reproduces W7-X with unit
  coefficients) is a diagnostic, and the averaged balance is gated.
- Reduced quadrature of the compression terms keeps the scheme non-variational:
  discrete eigenvalues are not upper bounds (p=1 cylinder fast modes fall
  slightly below the exact values). Near a shearless rational surface the
  cell-averaged Pfirsch-Schlueter inverse (#33) still amplifies the
  equilibrium's own resonant force-balance error; QAS at M = 12 keeps a
  D_Mercier dip near iota = 2/3 at s = 0.89 that VMEC does not show.

## Priority 0: deep source audit

The next audit must trace the complete production calculation from imported
equilibrium data to reported spectra and derivatives. Earlier component reviews
and passing tests provide evidence, but leave gaps in equation derivations,
interpolation accuracy, global derivative propagation, and benchmark agreement.
Start from the delivered production source at
[`6a27d11`](https://github.com/itpplasma/GLISS/commit/6a27d113d1cca29f641b68c83e2406a1b1ce1c68),
then freeze the exact commit and any worktree patch digest for each review.

- [x] Build a source inventory covering every production Fortran module, Python
  wrapper, C ABI entry point, application, benchmark runner, and CI verifier.
  Record each component's callers, mathematical contract, units, assumptions,
  independent oracle, derivative coverage, reviewer, and unresolved findings.
  Mark unread or unverified components explicitly; test counts alone do not
  establish audit coverage. [source_inventory.jsonl](source_inventory.jsonl)
  and [the inventory guide](docs/source_inventory.md) provide this accounting;
  whole-component contracts and equation reviews remain explicitly unverified.
- [ ] Have Astra xhigh lead an independent source review, with bounded parallel
  workers for geometry, assembly, eigensolvers, derivatives, and interfaces.
  Workers return findings and reproducible evidence; one controller owns the
  inventory, integration, PLAN updates, commits, and promotion to main.
- [ ] Trace reader conventions through half-grid coordinates, Fourier signs,
  winding, handedness, flux and pressure profiles, axis regularity, and Cartesian
  reconstruction. Audit the `s^(-m/2)` interpolation conditioning, between-knot
  geometry and derivatives, force balance, and sampled admission rules first.
  Use independent analytic maps and profile formulas to isolate input errors
  before interpreting the DCON discrepancy.
- [ ] Rederive each retained stiffness and physical-mass term from its cited
  energy principle. Follow signs, Jacobians, normalization factors, boundary
  terms, parity coupling, FEEC basis/traces, split quadrature, and radial scatter
  into the actual production assembly. Check every alternate public path and
  distinguish compatibility replays from independently assembled operators.
- [ ] Audit dense and block eigensolvers, equilibration, inertia, singular
  shifts, nullspaces, residuals, spectral selection, and error certificates.
  Use exact pencils, high-precision references, reordered and rescaled controls,
  and ill-conditioned masses. Separate rigorous eigenvalue/gap bounds from
  backward-error diagnostics, especially for clustered and slow modes.
- [ ] Trace JVP/VJP actions through every supported parameter path and list each
  missing link to geometry, profiles, global matrices, eigenspaces, and GVEC
  force balance. Check constraints, gauges, resonance branches, fixed topology,
  nonfinite/error paths, and cluster-boundary changes. Require independent
  directional finite-difference plateaus and full-space transpose duality;
  preserve differentiability on each declared admissible domain.
- [ ] Audit Python/C/Fortran ownership, lifetimes, array layout, integer limits,
  optional arguments, failure outputs, configuration migration, and provenance.
  Review test oracles and benchmark scripts for shared implementation mistakes,
  skipped controls, stale solver pins, and claims stronger than their evidence.
- [ ] Give every finding a minimal reproducer, severity, affected paths,
  independent expected behavior, and a linked existing or new GitHub issue.
  Fix confirmed bugs in small reviewed increments; retain unresolved physics
  discrepancies with explicit acceptance requirements in Priorities 1–3.

Completion requires an accounted-for source inventory, equation-to-code
evidence, independent review of fixes, regressions that fail before each fix,
and the applicable compiler, installed-package, derivative, and optimized
checks. Publish the coverage gaps and remaining findings; a clean audit report
must not imply that untested physics or an unfinished derivative chain is valid.

## Priority 1: qualify the production operator

Related issues: [#13, higher-order FEEC certification](https://github.com/itpplasma/GLISS/issues/13),
[#10, asymmetric input](https://github.com/itpplasma/GLISS/issues/10).

- [x] Reject reconstructed volume folds: the signed Jacobian must retain either
  consistent handedness across angular and radial assembly nodes. Independent
  polynomial-map controls exercise folds, both signs, and zero determinants.
- [x] Repair high-mode axis-regular interpolation conditioning
  ([#19](https://github.com/itpplasma/GLISS/issues/19)). The spline divided
  harmonics by the full `s^(|m|/2)` and amplified node roundoff by about
  `s_1^(-|m|/2)` (5e9 at m=24, 2e12 at m=28, 1e17 at m=36), which folded the
  M36 Solov'ev export. The first interval now retains full harmonic
  `s^(|m|/2)` regularity with a bounded polynomial quotient and a C2 join.
  Interior interpolation divides out only bounded parity factors. Native
  m=4/10/24/36 tests cover axis limits, derivatives and refinement. The
  issue's complete off-knot primitive torus sweep remains to be measured.
- [x] Replace fixed left-handed result metadata with the actual chart
  orientation ([#21](https://github.com/itpplasma/GLISS/issues/21)).
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
- [ ] Complete degree 1 through 4 certification: radial convergence on at least
  three meshes, global commuting projections and traces, independent assembled
  energy/stiffness oracles, and both parity classes.
  `test_radial_feec_convergence` measures the optimal L2-projection rates p+1
  and p on three graded meshes with a monomial oracle and a corrupted-map
  control ([#13](https://github.com/itpplasma/GLISS/issues/13)); the exact
  cylinder covers both parity classes. The Solov'ev toroidal sweep converges
  at order 3.3 for degree 3 after the conforming axis space of
  [#35](https://github.com/itpplasma/GLISS/issues/35); before it, the axis
  element under-integrated sqrt(s)-smooth integrands, hid a non-conforming
  |m| = 1 combination and used a mass unbounded on regular fields.
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
  genuine asymmetry couples parity classes and requires class 0.
  Verify symmetric full-storage files and an odd-harmonic cross-parity mass
  oracle. The full coupled operator supports asymmetric equilibria and retains
  cross-parity vacuum terms from asymmetric walls. Real symmetric QA and
  asymmetric LSP VMEC/BOOZ round trips and native spectra are now tested.
  The admitted operator is sampled at assembly nodes with dimensionally scaled
  tolerances; this does not certify continuum symmetry or angular convergence.
- [x] Give nonconstant magnetic differential-equation resonances an explicit
  policy in `src/export_surface_geometry.f90`. Without cell spread, the
  mean-projected equation rejects unresolved nonzero forcing. With finite cell
  spread it uses `D/(D^2+w^2)`; this regularized model does not prove solvability
  of the original resonant equation. Exact-resonance derivative tests retain
  pre-inverse forcing rather than attempting to reconstruct it from zero output.
  Resolved zero-width denominators retain the exact inverse.
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
| W7-X / Nuehrenberg 1996 | The documented coefficient-normalized L10 result is about -0.88701 versus digitized -0.37148 in the report's scaled units. The sparse quotient-aware solve now runs the L139 table size (139 labels) in 523 MiB (`benchmarks/cas3d_envelope`) | Resolve normalization, radial form functions, reference length, and unavailable deck details; these are not public |
| MISHKA / CASTOR | Branch transfer is unresolved; the CASTOR low-beta stable control currently fails | Transfer compatible invariant subspaces across at least three meshes, distinguish continuum branches, compare mass and decomposed potential energy under #12 |
| DCON / Solov'ev | Resolved by [#16](https://github.com/itpplasma/GLISS/issues/16): GLISS now counts 1 at 1.039062 and 0 at 1.039843 (public GVEC, ns64, M24), matching archived and freshly rebuilt DCON; the earlier `(1.05,1.10)` interval was a spurious non-conforming |m|=1 mode. Free boundary: the q0=1.5 critical conformal wall is 0.1519-0.1525 half-widths on refined edge meshes (m<=8) against DCON's 0.15249 | Extend to a converged q0 bisection and higher n; compare normalized energies; `test_solovev_axis_regularity` guards the sign on coarse public fixtures |
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
  The new `pressure_surface_derivatives` module completes the fixed-geometry
  pressure-sample spline-to-surface fields/drive/`gamma*p` map and its analytic
  JVP/VJP. Full nodal basis tests on nonuniform nodes, independent production
  finite-difference plateaus, duality and resonance controls pass. These
  tangents now propagate through bilinear angular products, radial FEEC
  assembly, axis congruence and sparse scatter to public isolated eigenvalue
  and complete-cluster pressure JVPs/VJPs. Imported geometry and resonance
  topology are fixed, samples must remain positive, and the VJP currently
  performs one tangent assembly per sample. Force-balanced GVEC geometry,
  conversion and vacuum design derivatives remain incomplete.
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
- [x] Test an installed wheel against independent oracles, not just
  import/version smoke checks. `ci/check_installed.py` enforces bundled-library
  loading, checks production inverse-density scaling and both material
  derivatives against central differences, and runs the Python suite from a
  temporary directory outside the checkout. Most Python tests are contract
  tests against a fake library; only tests marked `native`
  (`python/tests/test_native.py` and a few others) load libgliss_c, and they
  check the Solov'ev DCON stability bracket and Mercier sign and exact
  Rayleigh-quotient identities. `pytest -m native` selects them, and they fail
  rather than skip without the library. The cylinder is a synthetic ABI
  control, not a qualified straight-cylinder analytical benchmark
  ([#28](https://github.com/itpplasma/GLISS/issues/28)).
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
  `StabilityProblem(vacuum=VacuumModel(...))` adds the scalar-potential
  vacuum energy (exact toroidal-harmonic oracle at second order, with and
  without a wall), and the energy decomposition closes. Asymmetric walls
  retain full parity coupling. The historical q0=1.5 wall comparison does
  not qualify its imported equilibrium. Fresh M16 and corrected M32 exports
  failed the frozen qualification gates; their hashes, independent analytic
  checks and native edge values are retained in
  [the measurement record](benchmarks/results/2026-10-02/README.md).
  The M32 geometry passes, but native boundary flux derivatives fail.
  Faithful true-edge export/ingestion, exterior matching and an installed
  quantitative plasma-vacuum acceptance case remain required. Derivatives
  of the vacuum block remain part of #9.
- [ ] [#10: asymmetric VMEC and precomputed BOOZ_XFORM](https://github.com/itpplasma/GLISS/issues/10):
  implement full parity coupling and convention-complete round trips after the
  immediate input-admission work in Priority 1.
- [x] [#14: macOS wheels](https://github.com/itpplasma/GLISS/issues/14):
  the strict final source run passes both native builds and all twelve
  installed-wheel jobs; issue closure review is deferred at this handoff.
- [ ] [#15: documentation site](https://github.com/itpplasma/GLISS/issues/15):
  qualify actual publication and released-tag version pages.

The earlier GitHub audit found open issues #8 through #15 and closed
[#1](https://github.com/itpplasma/GLISS/issues/1), the earlier family/single-mode
assembly discrepancy. No issue was closed merely because this audit passed.

## Verification record

At the original baseline, controlled-thread `fo` passed build, tests, static
checks, and lint. The default unrestricted-thread attempt timed out; this is
not evidence of a physics failure. Python with the VMEC extra passed 237 tests
and skipped two optional cases. All eight original and hardened Enzyme gates
passed with Flang/LLVM 22 and LLVMEnzyme-22. A deliberately injected NaN
gradient failed the hardened drive gate as required.

The initial integrated corrections passed all 101 GNU Fortran CTest entries and all 109
Flang/Enzyme entries, with clean lint and changed-file formatting. An installed
wheel tested from a separate directory, explicitly loading its bundled native
library, passes 234 Python tests and skips SIMSOPT and the optional external
Mercier golden fixture. Three repository-registry checks require source files
and were excluded from that wheel run; they are structural checks, not physics
evidence. The source-tree Python run reports 237 passed and the same two skips.
The integrated GNU Fortran runtime-checking build also passes all stages with
`-fcheck=all -O1`; an earlier unoptimized runtime-checking snapshot passed too.

The parallel follow-up through
[`66ab128`](https://github.com/itpplasma/GLISS/commit/66ab128f89a9e40cb7595b5e0208ba16913316d1)
passes 104 GNU Fortran tests and 112 Flang/LLVM 22 tests, including all eight
Enzyme gates. The runtime-checking build passes with `-fcheck=all -O1`.
The installed wheel passes 245 portable Python tests, plus real production
angular-grid, density-scaling, and material spectral-derivative checks; the
source suite passes 248. Both Python runs retain the two optional skips above.
Correction ([#28](https://github.com/itpplasma/GLISS/issues/28)): those Python
counts were almost entirely contract tests against a fake library and are not
native physics evidence; the fixed-boundary reference eigenvalues were
self-regression values that moved from -79.144 to +1.19184 when the spurious
axis mode was removed (#16).
Independent analytical benchmark reruns reproduce the reported spectra and
convergence rates. The integrated marginality path reproduces count one at
the DCON stable endpoint and zero at q0=1.1; this preserves the unresolved
comparison rather than establishing agreement.

The strict optimized audit and GitHub lint initially rejected temporary arrays
in the new pressure finite-difference test calls. Explicit reusable perturbation
storage fixes those warnings without changing the oracle. Focused Flang and
runtime tests pass after that change, and a fresh isolated optimized audit
passes all 104 tests with array-temporary warnings treated as errors.
The optional compressible Solov'ev check timed out after ten minutes and
provides no physical-mass comparison result. No open GitHub issue was closed.

Before committing corrections, run the full `fo` pipeline, the native Python
suite, focused independent regressions, and the changed Enzyme gates. Freeze
uncommitted review inputs using the base commit and patch digest. Keep the
controller responsible for integration and promotion; reviewers return evidence.
