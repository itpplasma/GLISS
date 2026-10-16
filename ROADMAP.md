# GLISS qualification and implementation roadmap

Updated October 16, 2026. This is the active execution roadmap. [PLAN.md](PLAN.md)
retains the earlier audit, measurements and scientific limitations; completed
checks there apply to their recorded revisions, not automatically to new code.

## Target and supported domain

Deliver a physically qualified ideal-MHD stability solver with reproducible
FortSym-generated local mathematics, validated design derivatives, measured
performance and qualified figure-eight geometry. Qualification applies to a
declared domain of smooth, embedded, nested flux surfaces, supported boundary
models and locally fixed topology. A projected figure-eight may be embedded in
three dimensions; a physically self-intersecting plasma is outside this model.

Discrete residuals, symbolic identities, physical convergence and independent
solver agreement are separate evidence. Neither a generated kernel nor a green
test count certifies the complete calculation. Freeze acceptance tolerances
before measuring a candidate; retain failed controls and missing-data limits.

## Baseline

The initial production baseline is `e89d6a5343c6827433e3bacb9d9f47fa852ee282`,
development version 0.0.3, native ABI 4 and discretization revision 4. Hosted
GNU and Flang each passed 175 tests at this commit; source Python passed 329
with two optional skips and the installed wheel passed 324 with two skips.
Enzyme was skipped by automatic CI; eight derivative gates have earlier
source-pinned execution evidence. The source inventory has 311 records and
zero completed whole-component reviews.

Existing strengths include Cartesian primitive geometry, compatible radial
spaces, coupled parity, physical mass, OpenMP assembly, BLAS contractions,
block-sparse inertia/lowest-mode solves and installed native Python checks.
Missing qualification includes the complete equation audit, physical free
boundary, force-balanced design derivatives and the canonical figure-eight.

## Parallel workflow and promotion

One controller owns this roadmap, integration, test evidence and publication.
Managed workers own disjoint source paths, return reviewable changes and
evidence, and never promote their own work. Independent mathematical review is
required when a new claim is load-bearing; consensus does not replace proof.

Use one resident `fo gremlin` integration lane and the project's CMake/CTest
authorities. Workers can build independent generators in private scratch
directories. Run focused affected checks while developing; before promotion,
run the applicable native, Python, derivative, generation and optimized checks.
Reuse results for unchanged revisions and record unavailable checks explicitly.

Commit and push coherent verified increments to `main` regularly. Stage explicit
paths, inspect the diff and preserve unrelated work. Update achieved gates and
the next unresolved dependency after each integrated wave. Do not defer all
publication until the complete scientific program is finished. Do not close an
issue whose broader requested acceptance remains unmet.

## Work packages

| ID | Work package | Dependencies | Current execution |
| --- | --- | --- | --- |
| K | Generated physical mathematics | Independent equations and oracle | Physical-mass generator and derivatives in progress |
| G | Geometry and figure-eight admission | Declared embedding contract | Rigid-motion-invariant boundary admission in progress |
| V | Discretization and physical qualification | Qualified input and K | Production off-knot torus regression in progress |
| D | Efficient and complete derivatives | K, G, V for complete claim | Contracted pressure reverse assembly in progress |
| S | Sparse backend selection | Reliable backend capabilities | FortSparse readiness assessment; retain block-LAPACK reference |
| P | Performance qualification | Accepted workload and accuracy | Existing architecture retained; baseline instrumentation pending |
| R | Reproducible release and evidence | Applicable gates above | Roadmap/tracker reconciliation and integration checks |

### K: FortSym-generated kernels

- [ ] Independently derive physical mass and divergence from the cited energy
  principle with explicit units, signs, orientation and admissible domains.
- [ ] Add consumer-owned `gen/` sources and committed `src/generated/` output.
  Pin FortSym/FortGen identities and add reproducible regeneration/diff checks;
  ordinary GLISS builds require no CAS and never hand-edit emitted code.
- [ ] Emit values, coefficients and JVP/VJP contractions together where useful.
  Retain exact identity evidence before floating-point lowering; randomized
  probes are numerical evidence, not algebraic proof.
- [ ] Preserve independent Cartesian/metric energy oracles; cover all active
  argument directions, both orientations, symmetry, scaling and transpose
  duality. Cross-check symbolic derivatives against Enzyme and step plateaus.
- [ ] Extend the same source of mathematics to bending, shear, compression,
  drive, basis jets, geometric factors and needed boundary/adjoint leaves.
  Generic missing capabilities belong upstream, with a focused reproducer and
  consumer recheck; application physics remains in GLISS.

### G: geometry and figure-eight compatibility

- [ ] Complete exact-torus off-knot primitive Jacobian/derivative sweeps at
  M16/24/28/36 with fixed error bounds, radial refinement and corrupt controls.
- [ ] Replace global-z projected winding as decisive vacuum admission. Check
  the actual sampled piecewise-linear boundary for degeneracy, covering and
  self-contact, with rigid-translation/rotation controls and retained
  double-cover rejection. State the distinction from smooth off-grid embedding.
- [ ] Freeze a redistributable canonical figure-eight equilibrium, profiles,
  producer/export commits, conventions and input hashes. Qualify a moderate
  VMEC-admissible case separately from general GVEC Cartesian/G-frame geometry.
- [ ] Compare positions, fields, fluxes, metric/Jacobian and derivatives at
  independent off-knot points; converge radial and both angular resolutions,
  axis regularity, parity coupling and rational-surface treatment.
- [ ] Compare independent production stability operators on the same accepted
  figure-eight: signs/inertia, normalized energies and invariant subspaces over
  at least three useful resolutions. Record unavailable external data explicitly.

### V: physical and discretization qualification

- [ ] Complete equation-to-source contracts for retained production paths,
  including quadrature, axis congruence, parity, traces and scatter.
- [ ] Complete higher-order manufactured gradient/curl sequences, commuting
  projections/traces, independent dense element oracles and 3D convergence.
- [ ] Complete public-path cylinder fast/slow, Newcomb/Suydam and helical-limit
  controls with physical scaling, both stable/unstable cases and valid inputs.
- [ ] Resolve #35's degree-1 q0=1.035 instability criterion at 64 cells; repeat
  higher-degree convergence on qualified equilibria. Reduced compression
  quadrature does not establish general eigenvalue upper bounds.
- [ ] Converge the Solov'ev marginal interval and resolved QAS/LI383 controls;
  separate equilibrium error, regularized resonance behavior and FEEC error.
- [ ] For physical free boundary (#8), repair true-edge export/ingestion and
  qualify zero edge pressure, tangential-field continuity and current-free
  exterior matching before wall/angular/edge-mesh refinement acceptance.
- [ ] Obtain missing original CAS3D decks/normalizations (#11) and complete
  MISHKA/CASTOR branch/subspace transfer (#12). Matrix replay alone does not
  qualify the independently assembled physical operator.

### D: derivatives and optimization

- [ ] Replace pressure VJP's per-sample tangent assemblies with a contracted
  reverse assembly and the exact spline transpose; retain tangent/cardinal
  actions as independent small-case oracles and measure gradient cost scaling.
- [ ] Define design DOFs, constraints, units, gauges and fixed-topology domains.
  Connect converged force-balance implicit tangents/adjoints through geometry,
  chart/Boozer/conversion, plasma/mass/vacuum assembly and spectral objectives.
- [ ] Preserve isolated-eigenvalue and complete-cluster semantics. Declare
  spectral-gap and selection policy; reject unresolved cluster-boundary changes.
- [ ] Require finite-difference step plateaus, full-space transpose duality,
  nonfinite/failure controls and refinement at every layer and end to end,
  including fixed/free boundary and the accepted figure-eight.
- [ ] Expose actual equilibrium design variables in SIMSOPT and verify an
  optimization improvement with a fresh converged solve and independent control.

### S: meaningful FortSparse use

FortSparse native mode is progressing rapidly and remains under development.
Use it only where a supported capability benefits GLISS. Keep the validated
block-LAPACK path as reference and fallback; do not add an unused dependency or
replace a reliable inertia method with pivot-perturbed LU.

- [ ] Confirm a concrete consumer benefit before adding a backend. Inertia
  bracketing requires symmetric-indefinite factorization with reliable inertia,
  explicit singular/unreliable status and truthful treatment of any perturbation.
- [ ] When that API exists, pin the tested revision and compare exact/scaled/
  reordered pencils, 1x1/2x2 pivots, singular shifts, solves and residuals against
  independent dense references and the production block path.
- [ ] Validate refactorization, multiple RHS and required transpose/adjoint
  operations. Preserve the original pencil and equilibration congruence.
- [ ] Make an immature native route explicit and opt-in. Report backend identity
  and failures; never silently count failed native qualification as fallback
  success. Promote defaults only after correctness and meaningful performance
  evidence. Coordinate reproducible generic defects with the upstream owner.

### P: performance at accepted accuracy

- [ ] Record stage timings, peak memory, input/code/toolchain/library hashes and
  repeated-run distributions for qualified fixed/free-boundary and figure-eight
  values and complete gradients.
- [ ] Sweep radial/Fourier resolution, degree and OpenMP x BLAS thread counts;
  measure strong scaling and avoid oversubscription. Compare time to the same
  accepted observable/error target on disclosed hardware.
- [ ] Optimize measured bottlenecks: workspaces, geometry/phase reuse, SIMD,
  batched contractions, sparse factorization and reverse actions. Preserve
  independent behavioral oracles and unchanged physics to accepted tolerance.
- [ ] Add CUDA execution only where measured workloads justify it; generated
  leaf emission alone supplies neither device residency nor a solver backend.
  GPU qualification requires actual hardware execution, not a skipped test.

### R: evidence, CI and release

- [ ] Automate affected derivative and regeneration gates with pinned tools.
  Keep bounds/runtime checks, GNU/Flang, native installed-package controls and
  the optimized array-temporary audit; report optional skips distinctly.
- [ ] Publish revision-bound acceptance manifests and minimal redistributable
  fixtures, preserving failed controls and unavailable external inputs.
- [ ] Reconcile #14/#16/#17 closure-ready corrections and retain successor work
  for broader acceptance left after #13/#18/#19/#22/#33 bounded closures.
- [ ] Finish #10 metadata/version controls and installed LI383 acceptance, #30
  upstream GVEC reporting, and #15 immutable released-version documentation.
- [ ] Maintain this roadmap and the bounded source inventory at integration;
  completed file accounting must not be presented as completed equation review.

## First integration wave

The active wave targets generated physical mass, exact pressure reverse actions,
production off-knot torus controls and discrete boundary embedding admission.
Each can be promoted independently once its affected checks pass. A canonical
force-balanced figure-eight, missing external decks, upstream export repairs
and full design derivatives remain separate acceptance dependencies.
