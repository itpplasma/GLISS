# GLISS qualification and implementation roadmap

Updated October 17, 2026. This is the active execution roadmap. [PLAN.md](PLAN.md)
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
| K | Generated physical mathematics | Independent equations and oracle | Mass, magnetic/fluid energy, pressure reverse and Cartesian primitive geometry leaves implemented |
| G | Geometry and figure-eight admission | Declared embedding contract | Strong GVEC N2-12 is now primary; QUASR is a weak-shape control |
| V | Discretization and physical qualification | Qualified input and K | Explicit-edge ingestion and optional radial meshes implemented; uniform degree-1 criterion unmet |
| D | Efficient and complete derivatives | K, G, V for complete claim | Contracted pressure reverse passes; generic equilibrium response unavailable |
| S | Sparse backend selection | Reliable backend capabilities | Native LDLT lacks public inertia reliability API; retain block-LAPACK |
| P | Performance qualification | Accepted workload and accuracy | Pressure reverse and Cartesian value benchmarks implemented and measured |
| R | Reproducible release and evidence | Applicable gates above | 178 native, 388 Python and nine Enzyme gates pass locally; new Linux generation run pending |

### K: FortSym-generated kernels

- [ ] Independently derive physical mass and divergence from the cited energy
  principle with explicit units, signs, orientation and admissible domains.
- [x] Add consumer-owned `gen/` sources and committed `src/generated/` output.
  Pin FortSym/FortGen identities and add reproducible regeneration/diff checks;
  ordinary GLISS builds require no CAS and never hand-edit emitted code.
- [x] Emit values, coefficients and JVP/VJP contractions for mass and local
  magnetic/fluid energies, plus pressure/spline reverse leaves.
  Retain exact identity evidence before floating-point lowering; randomized
  probes are numerical evidence, not algebraic proof.
- [x] Generate Cartesian metric, signed Jacobian/coordinate derivatives,
  magnetic components/magnitude and second-form values with local JVP/VJP.
  Preserve range-safe evaluation, both signed charts, frozen NFP/orientation
  and independent all-input derivative and extreme-scale regressions.
- [ ] Preserve independent Cartesian/metric energy oracles; cover all active
  argument directions, both orientations, symmetry, scaling and transpose
  duality. Cross-check symbolic derivatives against Enzyme and step plateaus.
- [ ] Extend the same source of mathematics to bending, shear, compression,
  drive, basis jets, geometric factors and needed boundary/adjoint leaves.
  Generic missing capabilities belong upstream, with a focused reproducer and
  consumer recheck; application physics remains in GLISS.

### G: geometry and figure-eight compatibility

- [x] Complete exact-torus off-knot primitive Jacobian/derivative sweeps at
  M16/24/28/36 with fixed error bounds, radial refinement and corrupt controls.
- [x] Replace global-z projected winding as decisive vacuum admission. Check
  the actual sampled piecewise-linear boundary for degeneracy, covering and
  self-contact, with rigid-translation/rotation controls and retained
  double-cover rejection. State the distinction from smooth off-grid embedding.
- [ ] Freeze a redistributable canonical figure-eight equilibrium, profiles,
  producer/export commits, conventions and input hashes. Qualify a moderate
  VMEC-admissible case separately from general GVEC Cartesian/G-frame geometry.
  The primary candidate is the published GVEC N2-12 G-frame, rather than the
  weakly folded QUASR tutorial case. Require fresh full-boundary initialization
  and an independent off-knot boundary comparison before equilibrium acceptance;
  restarting at larger Fourier bandwidth can retain the old truncated boundary.
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
- [x] Expose an optional validated native radial mesh, preserving uniform
  defaults and carrying the same mesh through pressure JVP/VJP. Radius grading
  improves degree-1 approximation but does not satisfy the original uniform
  criterion or resolve the narrow marginal interval at degree 1.
- [ ] Converge the Solov'ev marginal interval and resolved QAS/LI383 controls;
  separate equilibrium error, regularized resonance behavior and FEEC error.
- [ ] For physical free boundary (#8), repair true-edge export/ingestion and
  qualify zero edge pressure, tangential-field continuity and current-free
  exterior matching before wall/angular/edge-mesh refinement acceptance.
- [x] Add an optional strict true-edge schema and independently qualified
  Solov'ev fixture. Explicit boundary geometry and native flux derivatives
  drive vacuum-interface construction; legacy half-grid input remains supported.
- [ ] Obtain missing original CAS3D decks/normalizations (#11) and complete
  MISHKA/CASTOR branch/subspace transfer (#12). Matrix replay alone does not
  qualify the independently assembled physical operator.

### D: derivatives and optimization

- [x] Replace pressure VJP's per-sample tangent assemblies with a contracted
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
- [x] Close #16's tangential-space correction and #17's lowest-positive
  bracketing correction after current native, runtime-check and Flang acceptance.
  Broader degree-1 physical convergence remains tracked in #35.
- [x] Close #14 after both macOS native builds and all twelve installed-wheel
  Python 3.9–3.14/architecture combinations pass at `c70f6b8`.
  Broader scientific acceptance remains in the work packages above.
- [ ] Complete native pointwise physical qualification beyond #10's passing
  installed conversion/representation checks, resolve LI383 input accuracy,
  and finish #30 upstream GVEC reporting.
- [x] Reject malformed/fractional/inconsistent VMEC metadata and unsupported
  BOOZ versions before conversion; retain producer and runtime versions separately.
- [x] Build archived 0.0.1/0.0.2 guides from their immutable tags and execute
  their examples against the corresponding installed native wheels. Hosted
  deployment, public version URLs and exact tagged guide provenance pass at
  `c70f6b8`, completing #15's released-version documentation acceptance.
- [ ] Maintain this roadmap and the bounded source inventory at integration;
  completed file accounting must not be presented as completed equation review.

## Integrated evidence and active dependencies

The first code wave includes three reproducible generated modules. At FortSym
`c9b1ec1a0eaf20d5fcde47853498eeb3508715ba`, the strict suites prove 41 mass,
103 energy and 26 pressure/spline identities, with zero probe substitutes or
failures. All three registered reproduction tests pass byte for byte. The
installed FortSym consumer exposed missing deployment/runtime dependencies;
the upstream repair also passes a relocated-package expression, derivative
and emission oracle. Ordinary GLISS builds remain CAS-free.

The Cartesian primitive geometry wave pins published FortSym
`2f0571e39f36a0606b3ae150c8b3b3df92ce6750`. Its new generator proves 778 exact
identities, bringing the four suites to 948 with zero probes or failures.
All four registered reproduction tests pass, including both geometry outputs.
An opt-in upstream ordered-output schedule preserves field evaluation range;
composed generated normal/contraction leaves preserve curvature range. The
production builder retains its finite-input and regularity guards. A shared
generated primal include permits ordinary optimized compiler inlining without
link-time optimization; the public generated module retains value/JVP/VJP APIs.
These are pointwise derivatives, with chart orientation and NFP held fixed;
they do not supply the generic force-balanced equilibrium-to-spectrum chain.
The new Enzyme geometry gate independently checks 372 derivative contractions
against the generated JVP/VJP and Cartesian finite differences, covering every
continuous local input in both signed charts and NFP 1/2/5. All nine integrated
Enzyme gates pass with Flang/LLVM 22.1.8 and Enzyme 0.0.277. Changing only a
private generated include correctly rebuilds this gate and fails its independent
geometric oracle; restoring the include passes. Source accounting now includes
Fortran includes and their references: 343 records, still zero completed
whole-component equation reviews.

At `8e51c65`, all 178 native tests and 388 Python tests pass with two OpenMP
threads; two optional Python checks are skipped. The isolated optimized audit
also passes clean Fo lint, strict Release compilation with array temporaries
treated as errors, and all 178 default tests. These checks supplement the
four exact-reproduction and nine Enzyme gates above. New Linux execution was
unavailable because the permitted host was shutting down; earlier Linux
evidence applies to its recorded revisions. Hosted generation checks remain
required for this new wave.

The first wave passed 175 GNU 16.1 Release tests with
`-O3 -Werror=array-temporaries`, 332 Python tests with two optional skips, and
175 GNU bounds/runtime-check tests. Mandatory runtime-check CI also passes.
The integrated edge/mesh/metadata wave passes 177 native tests and 388 Python
tests with two optional skips, with two OpenMP threads and a shared compatible
OpenMP runtime. Native backend and pressure-action regressions explicitly
exercise two workers. All eight Enzyme gates
pass with LLVM/Flang 22.1.8 and Enzyme 0.0.277 after deferring scheduling
optimizations until after differentiation. New controls cover
off-knot torus jets through M36, sampled boundary contacts and rigid motions,
all local energy input directions, spline transpose duality, pressure trace
rotation invariance and the exact physical-dilation spectrum/energy law.
Dilation tests do not implement an arbitrary equilibrium-design derivative.
The VMEC round-trip oracle retains its tolerances and frozen-reference checks,
while strict serialization equivalence uses a transform produced on the same
platform to avoid amplified upstream roundoff differences.
Hosted GNU, Flang, GNU runtime checks, source Python, installed wheel/sdist,
lint and documentation deployment pass at `c70f6b8`. Historical guides execute
eight examples per tag and retain exact tagged-source and fixture provenance.
The official isolated optimized audit also passes at `78d1ee8`: clean Fo lint,
strict Release array-temporary errors and all 177 registered default tests.

The pressure benchmark in `benchmarks/derivatives/` records repeated timings,
inputs, library/toolchain hashes and thread settings. On this Apple Silicon
host, five repetitions with a shared LLVM runtime measured median
reverse/cardinal costs of 14.54/290.68 ms with one OpenMP thread and
9.02/178.86 ms with two. Component agreement passes; improvement is about
20.0x/19.8x for this case. The runtime repair removes a reproduced corruption
caused by mixing GNU libgomp with OpenBLAS's LLVM libomp. CMake admits the actual
GNU/LLVM combination using a concurrent BLAS and OpenMP feature oracle, resolves
the locked conda runtime from its declared loader paths, and rejects conflicting
overrides. `OPENBLAS_NUM_THREADS=1` does not control OpenMP-built OpenBLAS;
the benchmark now queries its actual configuration (one/two threads outside
assembly in these runs). These figures do not establish general optimality or
qualified figure-eight performance.

The reproducible `benchmarks/geometry/` control compares the generated builder
with the unchanged `f2da7f4` routine on 4096 regular nonorthogonal jets and both
Jacobian signs. On the same Apple M3 Ultra/GNU16 host, eleven alternating pairs
of 8,388,608 calls measure 0.192735/0.198593 seconds, a generated/original ratio
of 0.97050. All 22 outputs agree within 1.45e-15 and consumed checksums match.
This is comparable serial value-kernel performance on this workload, without
LTO or fast-math; it does not establish application throughput or optimality.

Active dependencies remain explicit:

- Fo's configured `--no-build` CMake-reply race is repaired upstream at
  `2f979d1`, exact frozen CMake cache metadata at `614d951`, proven CMake path
  aliases at `a0e8dea`, header-only native action handling at `01a2acb`, and
  FortFront nested SELECT RANK scanning at `07bfbb4`. The installed Fo retains
  original cache help/type/value metadata, including whitespace and literal
  characters, so automatic OpenMP discovery remains distinct from an explicit override.
  The old driver fails the focused behavioral regression and the repaired one
  passes; the original GLISS resident lane passes all twelve affected gates at
  `78d1ee8`. The separate mixed-runtime BLAS failure exposed by the lane is
  repaired in GLISS's runtime selection.
  The generated geometry include exposed an unused native action that the
  old driver attempted to archive as `.o`. The upstream correction passes
  the maintained native/reference CLI oracle, including static/executable
  header declarations and include-only rebuilds, and the full original GLISS
  optimized audit at `8e51c65`.
  The updated resident integration lane is healthy and quiescent with all
  thirteen required gates passing at `8e51c65`, including generated geometry.
  Its scope is 13/178 tests; the separate full audit supplies default-suite
  evidence rather than treating the focused lane as complete verification.
- The public QUASR 189686 NFP=2 boundary/coils are a hash-pinned weakly folded
  control (SHA256 `4bba5c53be41d29798d93faa4ee8bad91fe57f3a8ce094648bd09c883e52038e`).
  A 20,000-iteration GVEC `c0dc66f` attempt fails its frozen coefficient-force
  and current/iota goals; the independently evaluated force/current RMS are
  approximately 0.00495/0.00219 against goals of 0.0001. The producer shares
  the mixed-runtime limitation, so no fixture or convergence claim is promoted.
  Serial residual reevaluation confirms the physical failure. GVEC's tutorial
  calls this case figure-eight-like; independent axis geometry shows much weaker
  folding than its N2-12 reference. The primary qualification target is now the
  actual GVEC N2-12 G-frame and boundary, with fresh serial producer solves.
  Independent guiding-curve measurements give inclinations of 13.43 degrees
  for QUASR and 50.47 degrees for N2-12, with converged numerical writhe near
  0.0534 and 0.4152 respectively. These are frame-curve measurements, not
  magnetic-axis measurements from a qualified equilibrium. Both sampled
  boundaries retain monotone cylindrical angle, so neither certifies the
  angle-reversing domain. The stronger approximately 76-degree configurations
  in the [figure-eight paper](https://arxiv.org/html/2411.16411v2) remain a
  separate missing-input target.
  N2-12's pinned file is byte-identical to the CC-BY-4.0
  [G-Frame simulation archive](https://zenodo.org/records/14714598), SHA256
  `e3417b64285408f86180fc4b0f0c85d2e08e2c4c999e838123ba9f2add0361aa`.
  Its paper setup uses M2/N10 and two degree-5 radial cells; the current GVEC
  CI smoke setup uses N5 and truncates the prescribed boundary. The archive's
  available 400-iteration log explicitly reports nonconvergence and ships no
  native restart state, so it is provenance rather than an acceptance oracle.
  Fresh N10 initialization reproduces the symmetry-projected native boundary
  to 1.02e-14 m; restarting from the N5 state retains a 4.66e-4 m discrepancy.
  A fresh M4/N10, four-cell degree-5 serial solve reaches coefficient force
  3.50e-9 against its 1e-8 goal, but its independent physical force RMS is
  0.01652 against the frozen 1e-4 goal. Increasing the diagnostic angular
  quadrature also raises the coefficient force to 2.94e-7; radial quadrature
  doubling alone has negligible effect. Thus coefficient convergence is not
  yet angular-quadrature convergence or physical equilibrium qualification.
  A private true-edge export of this actual strong state passes GLISS native
  reader/spline reconstruction and off-knot geometry at nine radial positions,
  37 poloidal and 73 toroidal angles. Sampled piecewise-linear boundary
  admission passes at 25x64, 33x96 and 49x128 angular resolutions. This is
  executable geometry compatibility evidence; the failed physical state is
  not promoted as a fixture or a stability-spectrum reference.
  A further bounded M6/N20, six-cell degree-5 continuation from the correct
  boundary with 27x99 angular quadrature stopped at 540 seconds/14,500 steps.
  Its coefficient force 8.59e-9 and boundary error 1.02e-14 m pass, while physical
  RMS 0.00964139 still fails. The physical residual worsens from the 5,000-step
  checkpoint's 0.00840334 despite decreasing coefficient force. This does not
  establish continuous-force convergence; no final state or spectrum is promoted.
  Actual native rho=0 stored axes retain the strong fold: inclination 49.474
  degrees for M4/N10 and 49.493 degrees for M6/N20, with poloidal-angle position
  spread below 6e-17 m. The two branches at the bowtie projection crossing are
  2.146 m apart in 3D. `axis-evidence.json` pins these state measurements;
  cylindrical azimuth remains monotone and physical qualification still fails.
  On the identical final state, doubling angular quadrature to 55x199 changes
  coefficient force by only 1.02e-5 relative; both values pass 1e-8. This does
  not support remaining angular underintegration as the physical-error cause.
  Source inspection identifies the solver norm as constrained, preconditioned
  Fourier/B-spline energy-gradient coefficients. Four selected admissible
  physical virtual-work projections agree with raw native gradients within
  4.98e-13 absolute; the energy integral agrees within 5.92e-15 relative.
  These cancellation-sensitive projections establish no force-definition bug,
  full residual decomposition or continuous-force error bound. The benchmark's
  `force-norm-evidence.json` preserves contracts, values and scratch provenance.
  Frozen-field enrichment to M10/N30 preserves position, field and energy to
  roundoff but exposes added toroidal-mode energy derivatives of -1.894e-5 and
  -1.617e-5. Centered differences with extrapolation agree within 3e-8 relative.
  Radial enrichment also exposes weak residual. The state is not stationary
  in these larger spaces; coefficient norms still supply no physical-error
  bound. `trial-enrichment-evidence.json` preserves the transfer and FD checks.
  One 600-second M10/N30, twelve-cell continuation then reaches coefficient
  force 5.73e-8 and physical RMS 0.003643: both fixed goals still fail, despite
  a 2.65x physical improvement and 1.12e-10 relative diagnostic-grid agreement.
  Boundary error stays 1.02e-14 m and profiles are preserved. The final state
  is not exported or promoted; `enriched-continuation-evidence.json` records
  the bounded attempt and controls.
  Producer coefficient convergence must also pass independent physical force
  residuals and resolution sweeps before any fixture or spectrum is qualified.
- Explicit true-edge Solov'ev export now meets frozen flux/transform/geometry
  bounds: errors 1.99e-16 Wb, 1.20e-14 and 3.75e-10, with a corrupted geometry
  control rejected. Native ingestion, round trips and boundary reconstruction
  pass. Interior-volume-to-edge convergence, tangential-field continuity and
  a matched current-free exterior remain required for physical free boundary.
- Uniform-flux degree 1 at 64 cells gives positive lambda for q0=1.035;
  radius grading detects that instability but still misses the narrower
  marginal endpoint. No axis-quadrature defect was found: an independent
  Cartesian/logarithmic energy oracle agrees to 3.4e-16 relative. Fresh explicit
  exports at q0=1.039062/1.039843 preserve the unstable/stable bracket under
  serial execution: radius-graded degrees 2/3/4 give negative/positive lambda,
  as does the unchanged uniform degree-2 default. Radius-graded degree 1 gives
  positive lambda at both endpoints, so #35 remains open. The independent
  volume force-balance indicator is approximately 1.70e-5 and is separate from
  the producer's reported coefficient-force tolerance of 1e-10.
- GVEC exposes equilibrium energy/force but no raw constrained residual
  Jacobian, transpose or implicit design-response API. Those interfaces and
  differentiated export/reconstruction remain required for generic gradients.
- FortSparse native LDLT currently has no public reliable-inertia interface;
  adding it to GLISS now would not improve the required certified solves.
- The installed Linux wheel passes genuine asymmetric VMEC/BOOZ conversion,
  stored replay, native round trips and coupled-parity spectra. These establish
  representation consistency. Independent Cartesian curl calculations following
  the exported spline contract reduce sampled QA/asymmetric-LSP force residuals
  from 0.0343/0.0120 to 0.00266/0.00214 with angular and radial refinement.
  These use `Aminor*rms(curl(B) x B)/rms(|B|^2)`, a different normalization from
  the N2 protocol; they exclude the axis patch and edge and still need direct
  native pointwise comparison. LI383 Jacobian mismatch converges near 0.0858
  against its unchanged 0.03 gate, with independent radial alternatives still
  above 0.05. A resolved producer is required before assigning a converter fix.
- The external `gvec-stability` repository is now accessible at `822b2bb`;
  recovered Brezellator and historical TERPSICHORE assets supply bounded
  controls, not a qualified common figure-eight. Original CAS3D decks and the
  chosen common VMEC/GVEC figure-eight reference remain missing.
- Hosted macOS wheels pass for the edge/runtime wave at `c70f6b8`, including
  both native architectures and all twelve installed Python environments.
  Linux regeneration
  CI exposed a runtime-library discovery defect before any proof. The repaired
  workflow passes all 170 exact identities, byte-exact outputs and 3/3 CTests
  on an x86_64 Linux host with the locked GNU 15.3 Nix toolchain. Hosted
  reproduction now also passes at `c70f6b8`.
