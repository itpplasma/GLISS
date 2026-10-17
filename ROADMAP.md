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

## Highest priority: analytical and published-reference validation

The current priority is to validate GLISS's production physical operator
against analytical cases and published literature, explicitly including
CAS3D. Run the analytical and CAS3D tracks in parallel. Further kernel
generation, design derivatives, backend changes and performance work follow
these gates, except for repairs needed to execute or explain a comparison.
GVEC figure-eight convergence does not block the analytical GLISS track.

The active claim to establish is that the current production energy, mass,
assembly and solver reproduce independent physical observables on declared
test domains.
The verifier is the native CMake/CTest suite plus fresh, revision-bound runs
of the analytical and cross-code benchmarks below. Local formula checks,
stored-matrix replay, symbolic identities and small discrete residuals retain
their narrower meanings; none alone passes this claim.

| Gate | Existing evidence | Required work and acceptance |
| --- | --- | --- |
| A0: equation and normalization audit | Exact pinch raw energy/mass pass; matched W7-X trial exposes a raw/production energy gap on a force-defective retained equilibrium; zero completed whole-component reviews | Qualify equilibrium over the actual trial support, converge raw three-dimensional energy and P0 interface terms, and audit units, flux signs, axis constraints and boundaries |
| A1: analytical spectra and thresholds | Cylinder fast/slow/Alfvén branches, imported exact vacuum torus and 144 exact resonant screw-pinch negative-energy trials pass | Add global production Newcomb/Suydam and helical unstable/marginal limits; retain surface-assembly scope for the screw pinch and stable-domain scope for the imported vacuum torus |
| A2: axis, resonance and small eigenvalues | Higher-degree Solov'ev bracket controls pass; degree 1 still misses a tested marginal endpoint; dense slow-cluster ordering sensitivity is diagnosed | Resolve #35's frozen degree-1 criterion and retain the failed gate until it passes; test rational-surface domain and refinement; check slow-cluster accuracy using independent higher-precision pencils before attributing any failure to the production solver |
| A3: CAS3D formulation and numerical comparison | Independent dense quotient oracle reproduces floor inconsistency and zero-carrier component-topology defects; integrated repairs and ordered-CSV experiments pass their software gates | Recover and audit the full published mode mask; resolve provisional numerical discrepancies and compare matched signs/inertia, marginal parameters, energies and subspaces |
| A4: independent published-code comparisons | DCON sign bracket, TERPSICHORE replay and bounded independent QAS3 controls exist | Converge the Solov'ev marginal interval and higher-n controls; rerun current independent QAS3 assembly with qualified input; explain count/branch differences through energy and subspace comparisons; complete ideal-limit MISHKA/CASTOR transfer and the failed stable control under #12 |
| A5: physical free-boundary references | Analytic vacuum-energy oracle and true-edge ingestion pass | Qualify edge pressure, tangential-field continuity and exterior matching, then reproduce the DCON critical-wall interval and applicable published helical/wall limits with separate volume, edge and wall refinement |

Use the original [Bernstein energy principle](https://doi.org/10.1098/rspa.1958.0023),
[Newcomb pinch analysis](https://doi.org/10.1016/0003-4916(60)90023-3),
[Fu vertical-stability limit](https://bp-pub.pppl.gov/pub_report/2000/PPPL-3368.pdf),
[Schwab CAS3D formulation](https://doi.org/10.1063/1.860656),
[Nuehrenberg W7-X study](https://doi.org/10.1063/1.871924) and
[compressional CAS3D formulation](https://doi.org/10.1063/1.873268).
Freeze exact equation/figure references and their physical assumptions before
deriving each oracle. The 1996 coefficient norm and the 1999 physical kinetic
norm are distinct; raw eigenvalues from unlike norms are not comparison data.

The CAS3D dense oracle removes redundant labels by orthonormalizing the range
of the independently specified coefficient map. It checks the coefficient
pencil, not a different physical kinetic norm. The old returned zero floor
was about 90 times smaller than the floor used by the solved coefficient
pencil in a reproduced case. The repair uses one pencil and floor for inertia,
mode selection and reporting. A second reproduced defect concerns the absent
sin(0) trial: normal and eta components can have different active mode counts.
The sparse congruence now accepts their separate active rows. Exact map rank
is distinguished from stiffness rank, which can also contain physical neutral
modes. Keep near-zero positive/negative and redundant-label controls.
Origin-only envelopes now retain valid normal-only physical blocks without
empty trailing groups. Their positive eigenpairs agree with independent dense
oracles in both norms. The complementary exact symmetry nullspace returns
zero negative inertia and no resolved eigenpair; singular or indefinite mass
is still rejected.

The runner now accepts ordered CSV modes, explicit printed-n sign and L/L+1
prefix conventions, both supported single parities, separate assembly and
coefficient grids, reference length, native radial cells and quadrature. It
freezes assumptions and input/library identities before solving, preserves
failed hypotheses, and fixes Darwin RSS units and the unsupported `--surfaces`
usage. Coupled parity is not implemented for this coefficient norm.

Proceed with reasonable declared assumptions rather than waiting for every
historical setting. The current controlled W7-X experiments freeze carrier
(3,2), degree 1, 48 radial cells, assembly 96x64, coefficient grid 36x24 and
reference length 10, then vary both parities, both quadratures and both prefix
conventions. These are hypotheses, not recovered settings or fitted values.
The first eight envelope entries can be reconstructed exactly from a retained
native L8 control; their full Figure 6 transcription is not independently
verified. The pinned research server currently refuses connections, and the
primary PDF is inaccessible. Do not invent later mode entries.

All 24 provisional L3/L5/L7 baseline runs completed. Of 16 L3/L5 comparisons
with available digitized values, eight miss the frozen broad tolerance and
all miss the approximate digitization band. Caption/Gauss values in recorded
10^3 lambda units are -0.190629 and -0.261329 versus -0.14613 and -0.16725;
caption/midpoint gives -0.130592 and -0.187850. Both parities agree, but that
does not resolve the paper discrepancy. An independently calibrated Cartesian
calculation on five interior radii gives
rms(curl(B) cross B - mu0 grad(p))/rms(|B|^2) = 0.343274216 per metre.
Differentiation-step and angular refinements converge; native aggregate
geometry/field comparisons agree to 3.5e-15. This is an angle-sampled interior
diagnostic, with no axis/edge, volume weighting or native pointwise-curl check;
it does not qualify this frozen export. A fresh solve of the unchanged complete
VMEC input converges below its 1e-12 coefficient-force target and preserves the
boundary, pressure, flux and zero-current settings. The current converter
rejects m5/n4 at its unchanged Jacobian gate. Fresh m7/n7 and m10/n10 exports
pass that gate; their independently sampled characteristic Cartesian force
residuals are 0.0746074 and 0.0180241 using the producer's minor radius. Neither
is physically converged. Fixed BOOZ72 and retained m10/m14/m18 controls reduce
this norm to 0.00519838; that misses the frozen 0.005 gate. Keeping the retained
18/18 representation and using all 144 native radial surfaces reduces it to
0.00425274 and passes that bounded sampled gate. The residual is nevertheless
21.68% of the pressure-force RMS; neither gate covers the axis, edge or global
maximum. Resolution-only source solves at MPOL/NTOR 13/11 and 17/15 converge
below the unchanged native 1e-12 target and preserve the prescribed physical
profiles and boundary. At fixed BOOZ48, retained 18/18 and 144 surfaces their
characteristic residuals are 0.00150060 and 0.00156457, with pressure-relative
residuals 7.65% and 7.98%. Both fail the separately frozen 0.001 and 5% physical
limits. The final source-pair change passes 0.0005, but retained-representation
convergence is not yet isolated by those controls. A separately frozen 2x2
retained-22/26 study on both sources improves the final 17/15, retained-26
case to 0.000874519 characteristic force and 4.458% pressure-relative force.
It passes both strict sampled limits; source and retained-pair changes pass
0.0005. The other three individual cases still fail the physical limits.
This qualifies only the five specified interior radii and aggregate angle
norms, with no axis/edge or pointwise maximum claim. Vary equilibrium resolution
separately from the stability mesh. The original profile diagnostic's large
maximum is not itself a fractional force error.

All 26 fresh-export CAS3D controls complete with the original declared
assumptions. The m7/n7 to m10/n10 eigenvalue drift is at most 0.0283%, with
unchanged inertia; the fresh input shifts values by 2.15–3.03% from the archived
export. All 18 comparisons still miss the approximate digitization band.
For m10/n10, caption/Gauss L3/L5 give -0.185670 and -0.254027; midpoint gives
-0.127052 and -0.182200. Retain these failed comparisons while investigating
the independent mode mask, radial trial functions and normalization conventions.
Preserve the historical L10 discrepancy (-0.88701 versus -0.37148) separately.
Issue #11 remains open.

Fourteen fixed-BOOZ72 CAS3D controls retain six completed solves and eight
undersampled failures. The existing cubic-Jacobian MDE admission requires
85x85 angular points for retained 14/14 and 109x109 for 18/18. The public
validator now reports those requirements instead of a generic assembly error;
the original failed runs remain evidence. At 18/18, assembly refinement from
192x128 to 384x256 changes the lowest eigenvalue by less than 5e-8 relatively.

With the same 144-surface equilibrium, all sixteen degree-1 radial FE controls
at 48/96/192/384 cells complete. The final eigenvalue changes are below 0.203%
and Gauss/midpoint differences below 0.083%. The reported negative count drops
from two to one at 384 cells while the numerical zero floor grows; the frozen
mesh gate therefore fails. Independent sparse eigenpairs on six actual L3
pencils retain two disjoint negative residual intervals. At FE384 the second
mode is approximately -1.81e-5, inside the official 3.50e-5 marginal band;
the public count correctly excludes it. No floor or count policy is changed.
Its Gauss eigenvalue still changes 2.12% on the final mesh pair, so convergence
of that weaker mode remains unestablished. L5 is outside this independent
reassembly scope. All sixteen literature
comparisons still miss the cached digitization bands. The final L5/L3 ratio is
about 1.3487, with certificate-expanded interval [1.3221,1.3758], versus the
cached reference interval [1.0451,1.2542]. A common positive reference-length
or eigenvalue rescaling cannot align both values. These references remain
conditional on the unrecovered primary figure and full mask; do not fit them.

Eight unchanged CAS3D comparisons on the improved source13/11 and 17/15,
retained26, 144-surface inputs all complete. At FE192 and assembly192x192,
all four source-pair eigenvalue changes pass the frozen 1% gate (0.0974–0.1202%)
with negative count two throughout. The sampled-force-qualified 17/15 input
gives L3/L5 values -0.203176/-0.273623 for Gauss and -0.203844/-0.274375 for
midpoint, in recorded 10^3 lambda units. All eight conditional literature
bands still fail, including numerical certificates and reference rounding.
These source-angular controls do not establish force accuracy over the full
operator domain. Recover the mode/radial conventions and compare raw energy
on matched trials; no normalization is fitted. One FE mesh and one angular
grid for this improved representation do not establish full convergence.

The matched-trial investigation now takes priority over further CAS3D sweeps.
Check force balance throughout the actual assembly support, including the
near-axis and outer cells, and separate source, transform and radial
interpolation errors. Compare Cartesian curl current with reconstructed MDE
current, then converge raw Bernstein energy for the same physical trial,
including the one-sided P0 interface contributions. Retain the existing
five-radius protocol's bounded result; it does not certify this larger domain.
Repeat the frozen CAS3D comparison only after this operator/input consistency
check is qualified. Historical settings remain declared hypotheses.

[Matched-trial evidence](benchmarks/cas3d_envelope/raw_energy_comparison.json)
uses the actual lowest FE192/L3 midpoint vector on the final 17/15, retained26
input. Independent production replay agrees with its native matrix quadratic
to 3.96e-7 relatively. Full-torus raw/production 2 delta W is
-631.992/-811.071 J with midpoint integration and -627.226/-805.852 J with
two Gauss points per cell. The refined magnitude gap is 28.48% of the raw
value; angular refinement changes raw energy by 1.89e-10 relatively, while
the radial change is 0.754%. Sampled volume-weighted force over all 192
assembly midpoints is 0.006033 characteristic and 28.31% pressure-relative.
The first two physical vector-force RMS values are 0.287 and 0.214 T^2/m;
the second is outside the special axis patch. No pressure/flux unit defect
or missing Bernstein term was reproduced. Current, force substitution and
drive differences close the production/geometric gap algebraically, but
signed volume/interface corrections and continuum energy are not converged.
The trial has finite-energy axis behavior in inspected probes; Cartesian
smoothness and full operator-domain equilibrium qualification remain open.

[The source radial control](benchmarks/cas3d_envelope/source_radial_comparison.json)
changes only the final VMEC NS from 145 to 289, retaining MPOL17/NTOR15,
the physical deck, BOOZ48 and retained26; all 288 half-grid surfaces are
exported. The original 180-second timeout is preserved. One declared
600-second retry converges in 399 seconds with all three native force
residuals below 1e-12, and the unchanged public conversion gates pass.
On the frozen FE192 midpoint trial, sampled volume-weighted characteristic
force improves from 0.006033 to 0.004182 and pressure-relative force from
28.31% to 19.62%. Raw/production energy becomes -646.227/-692.575 J,
reducing the relative gap from 28.34% to 7.17%. All three preregistered
whole-support targets still fail: characteristic force <0.001,
pressure-relative force <5%, and energy gap <1%. No radial cutoff is used;
the first two midpoint characteristic forces remain 0.03244/0.02073.
The same physical x/eta coefficients are retained in normalized-flux
coordinates; Cartesian displacement changes with the equilibrium geometry.
This is one source-resolution control, not a new eigenpair or continuum
convergence result. The two conditional CAS3D solves are not executed.
Identify and repair the cause, then repeat these unchanged physical gates
before the frozen CAS3D pair.

[Matching-shell reconstruction diagnostics](benchmarks/cas3d_envelope/reconstruction_comparison.json)
compare four flux surfaces s=1/288, 2/288, 4/288 and 8/288. With native stored
current, Jacobian-weighted characteristic force is 1.29e-4 to 2.66e-4.
An independent continuous native-wout field gives 2.54e-4 to 6.31e-4;
the public flux/Cartesian-spline reconstruction gives
0.03343/0.03097/0.01922/0.0006424. Public/native-continuous ratios are
52.97/88.57/38.55/2.53. Native interpolation has nonzero divergence and
metric mismatches, and its 64x64 angular grid differs from the public 192x192
grid; this localizes a substantial extra reconstruction discrepancy without
qualifying the source or identifying a specific production defect.
The installed BOOZ reader's odd-mode sqrt(s) staggering is correct on both
frozen sources. Replacing only low-mode parity factors with full |m|/2
factors worsens force on all four shells by factors 2.33 to 125.93;
independent baseline values/fields match the native trace, so that proposed
repair is rejected. Next separate BOOZ geometry, Cartesian Fourier truncation
and continuous radial reconstruction, retaining the failed whole-support gates.

The old synthetic 70-envelope/139-label run measures scaling. A literal paper
L139 prefix would contain 139 or 140 envelopes and 277 or 279 labels before
collisions. Do not identify these different dimensions as the same benchmark.
Complete integrated A1/A3 checks and explain the controlled discrepancies
before expensive high-L or further figure-eight studies.

The [provisional CAS3D evidence](benchmarks/cas3d_envelope/provisional_comparison.json)
retains all 24 outcomes, frozen settings, input/library hashes and unavailable
data. Its [eight-entry control](benchmarks/cas3d_envelope/native_l8_envelope_control.csv)
is explicitly derived from the retained implementation test.
The separate [refinement evidence](benchmarks/cas3d_envelope/refinement_comparison.json)
retains the fixed-transform, radial-FE and equilibrium-source controls with
their unchanged acceptance limits and failed outcomes.

Every accepted case must record primary reference, independent oracle, input
and code identities, units/norms, mode and boundary conventions, frozen error
bounds, failed controls, and at least three useful refinement levels. Vary
equilibrium resolution separately from stability discretization. Use sign,
inertia and marginal intervals when frequency normalization is unmatched;
compare isolated modes or complete clusters under a shared physical inner
product. Fix reproducible GLISS errors before extending claims. Record
missing external data and untested behavior explicitly.

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
| A | Analytical and published-reference acceptance, including CAS3D | Independent equations, oracles and qualified inputs | Highest priority; execute A0–A5 above before broader optimization claims |
| K | Generated physical mathematics | Independent equations and oracle | Mass, magnetic/fluid energy, pressure reverse, primitive geometry and Cartesian field/curl leaves implemented; new field parameter AD remains open |
| G | Geometry and figure-eight admission | Declared embedding contract | Strong GVEC N2-12 is now primary; QUASR is a weak-shape control |
| V | Discretization and physical qualification | Qualified input and K | Explicit-edge ingestion and optional radial meshes implemented; uniform degree-1 criterion unmet |
| D | Efficient and complete derivatives | K, G, V for complete claim | Contracted pressure reverse passes; generic equilibrium response unavailable |
| S | Sparse backend selection | Reliable backend capabilities | Native LDLT lacks public inertia reliability API; retain block-LAPACK |
| P | Performance qualification | Accepted workload and accuracy | Pressure reverse and Cartesian value benchmarks implemented and measured |
| R | Reproducible release and evidence | Applicable gates above | 184 native, 407 source Python, 382 installed Python and five generation gates pass locally; prior unchanged AD kernels have nine Enzyme gates |

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
- [x] Generate Cartesian B, its spatial gradient, curl and divergence from
  primitive jets. Retain independent curved-chart, vacuum-torus, screw-pinch,
  signed-chart and finite-difference oracles and exact regeneration checks.
- [ ] Add parameter JVP/VJP and independent AD checks for the new Cartesian
  field/curl API; spatial differentiation alone does not satisfy this gate.
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

The analytical/CAS3D waves add six registered native tests. All 184 default
CTest cases at `ddaa8f6` pass in the strict GNU 16 Release build. All 407
Python tests pass with two optional skips and BOOZ_XFORM installed.
Source accounting contains
354 records and still zero completed whole-component equation reviews.
At published `ddaa8f6`, the isolated Fo audit also passes clean lint, strict
optimized compilation with array temporaries treated as errors, and all 184
then-registered tests. At published `fbb2430`, hosted GNU, Flang, runtime-check,
Python, lint, format, generation and docs checks pass, as do both macOS native
builds and all twelve installed-wheel Python/architecture combinations.
The resident integration
generation at `ddaa8f6` passes its eleven selected gates;
that selected scope is not a second full-suite claim. The immutable Fo driver
is `454c93dffd513b85bc8274c262a2d5364b20e4acf24c45f6ecd55de1769f2daf`.
The analytical cylinder identifies two fast/slow branches and the Alfvén
cluster for both parities, degrees 1–4 and 8/16/32 cells, with zero negative
inertia and independently bounded branch enclosures. Finest relative errors
are 9.76e-5, 3.05e-7, 1.79e-10 and 2.89e-14; six physical variants and exact
scaling laws also pass. This exercises production surface assembly and the
full-spectrum certifier, while retaining its prescribed-field scope.
Generic reference BLAS/LAPACK exposed a degree-4 branch whose seed satisfied
the solver's backward-error stopping rule but missed the stricter physical
mass-residual target. A local reference-LAPACK reproduction confirms that
indexed production refinement reduces the bound from 9.48e-4 to 1.99e-6,
below the unchanged 7.99e-4 target. The test reruns all spectrum certificates
after replacing a branch vector; accuracy budgets remain unchanged.
The installed-distribution verifier excludes the repository-only benchmark
runner test, which remains in the source suite. A rebuilt wheel passes all
382 installed tests plus production angular/scaling/derivative checks; a
fresh sdist builds and loads its bundled library outside the checkout.

The finite-pressure screw-pinch oracle starts from Bernstein et al. (1958),
printed page 23, equation 3.13, with independent cylindrical vector calculus
and polynomial integration. Twenty-four cases exercise current, pressure
gradient, fluid compression, absolute mass, both chart signs and both angular
assembly routes. The worst relative error is 7.51e-16 against the frozen 4e-12
limit. This checks surface assemblies, not a public primitive cylinder import
or global unstable threshold. It also reproduces and repairs division by zero
in the three-component surface variable nu=mu-(FP'/FT')*sqrt(g)*eta: zero FT'
now fails admission before mutation. Perpendicular and generated point kernels
retain their independently tested wider domain.

The resonant screw-pinch test adds an independently derived negative-energy
trial at k=-m A/Z. Direct cylindrical curl gives Q=-B div(xi); minimizing the
raw Bernstein density over that divergence gives -p_r^2 f^2/H, H=B^2/mu0+gamma p.
The explicit displacement is perpendicular, Cartesian regular at the axis and
zero at the wall. All 144 configurations pass the unchanged 4e-12 budget;
the largest total-energy error is 1.37e-12. Gamma=1 has a closed beta-function
integral, and independent 80-digit integrals provide energy decomposition and
mass for gamma=1 and 5/3. OpenBLAS and reference BLAS/LAPACK controls pass.
This supplies an unstable prescribed trial, not a global eigenvalue threshold.

The imported exact vacuum torus exercises the native writer/reader, Cartesian
primitive spline, public fixed-boundary assembly, lowest-mode and full-spectrum
solves. Positive Jacobian and GLISS's negative toroidal-flux slope give the
known field B=F/R e_phi with zero pressure/current. For xi^s=s(1-s), the absolute
quadratic magnetic-energy error is 3.56e-12 and the independent physical-mass
error is 1.76e-13. Three import refinements, both parities, coupled spectra,
neutral modes, B-squared/density scaling and invalid-input controls pass.
This is a stable finite-aspect-ratio reference, not an unstable threshold.
See [analytical benchmark scope](benchmarks/analytic/README.md).

The CAS3D algebra gate uses independent product-to-sum matrices and an
orthonormal dense quotient, including colliding and self-colliding labels,
both parities, independent normalization changes and positive/negative
near-zero selection. It verifies active map rank/nullity, mass normalization,
residuals and certificates after the two reproduced repairs above. Fresh
corrected-library L3/L5 checks retain the provisional literature discrepancy;
solver-reporting repairs do not turn it into agreement.

The Cartesian field/curl wave adds 116 exact identities at the same pinned
FortSym revision: five suites now prove 1,064 identities with zero probes or
failures, and all five registered reproduction tests pass byte for byte.
The native point and spline-grid API returns Cartesian B, its spatial
gradient, curl and divergence. Independent curved-chart, exact vacuum-torus,
finite-pressure screw-pinch, both-orientation, magnetic-null and range controls
pass. On 11,520 interior points per archived/m7/m10 export, native positions
and B agree with independent Python reconstruction to roundoff; three-step
Cartesian finite differences converge quadratically to the native gradient
and curl. Native divergence stays below 1.34e-15 T/m. These checks establish
the sampled representation mathematics, not force balance. Parameter JVP/VJP
and C/Python exposure for this new API remain open.
The improved source17/15 retained26 input passes the same native/Python
parity protocol; the finest curl relative RMS difference is 6.33e-6 and the
finer-step convergence order is two. Its outer sampled surfaces still exceed
the aggregate physical limits individually (16.82% pressure-relative residual
at s=0.85); the aggregate pass must not become an every-surface claim.

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
four exact-reproduction and nine Enzyme gates above. Direct execution on the
permitted Linux host was unavailable while it was shutting down. The hosted
[Linux generation run](https://github.com/itpplasma/GLISS/actions/runs/38973903951)
subsequently passes all four reproduction tests at `8e51c65`. Hosted CI,
documentation and macOS wheels also pass at `8e51c65` and the documentation-only
`fa72f97`. These software gates do not complete A0–A5 physical acceptance.

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
