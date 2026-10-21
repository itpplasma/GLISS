# CAS3D2MN coefficient envelopes

`run.py` solves the CAS3D2MN labeled phase envelope with the coefficient
normalization (`gliss.solve_cas3d_phase_envelope(...,
normalization="cas3d2mn_coefficient")`) for envelope tables of increasing
size, each in a fresh process, and writes a manifest of labeled and physical
dimensions, quotient rank, labeled nullity, widest block, inertia, lowest
eigenvalue, certificate, eigenpair residual, equilibrium force-balance
residual, numerical zero floor, peak resident memory and wall time.

```sh
GLISS_LIB=build/libgliss_c.so PYTHONPATH=python OMP_NUM_THREADS=2 \
    python benchmarks/cas3d_envelope/run.py cd_w7x.nc l139_manifest.json
```

The labeled pencil is the congruence of the physical block-tridiagonal
pencil by a block-diagonal map, applied block by block
(`apply_cas3d_phase_envelope_block_congruence`): no dense matrix is formed
at physical or labeled size. Coincident labels and inactive `sin(0)` trials
give an exact null space of the coefficient map (`labeled_nullity`);
`quotient_rank` measures the map rank, while physical neutral modes can further
reduce stiffness rank. Exact zero-shift inertia agrees on the quotient.
Finite marginal-band counts, mode selection and the reported zero floor all
use the coefficient pencil and its units. A floor from a different physical
mass pencil is not interchangeable with this one.
`test_cas3d_quotient_spectrum` checks small colliding, noncolliding and
zero-carrier cases against an independently specified dense quotient,
including positive/negative near-zero modes and norm scaling.
`test_cas3d_block_congruence` checks the blockwise congruence
entrywise against the retained dense congruence on a colliding table and
the quotient rank and nullity exactly.

## Ordered literature experiments

The default mode preserves the original synthetic envelope and parameters.
Its 70 envelope entries give 139 sideband labels; that size does not identify
the published mode mask. For a paper comparison, supply an ordered CSV with
integer `m,n` columns and `(0,0)` as its first row. Additional columns, such
as the printed cell number or transcription notes, are allowed. Rows are
never sorted, deduplicated or filled in. Duplicate pairs and invalid input
are rejected.

For example, after auditing a transcription against the original figure:

```sh
GLISS_LIB=build/libgliss_c.dylib PYTHONPATH=python OMP_NUM_THREADS=2 \
    python benchmarks/cas3d_envelope/run.py w7x.nc caption.json \
    --envelope-csv figure6_mode_mask.csv --envelope-n-sign opposite \
    --literature-l 3 5 7 --prefix-convention caption_l \
    --base 3 2 --parity 1 --angular 96 64 \
    --coefficient-angular 36 24 --reference-length 10 \
    --radial-cells 32 --radial-quadrature gauss5 \
    --assumption 'Original coefficient grid unavailable: use 36 x 24' \
    --assumption 'Original reference length unavailable: use 10'
```

These numerical choices illustrate declared hypotheses; this command does
not assert that the unavailable historical settings have been recovered.
Use the platform's actual shared-library suffix. The CSV path must point to
an independently audited ordered transcription. A table reconstructed from
an old GLISS test is implementation-derived evidence and must be recorded
as such; it is not an independent paper transcription.

`--envelope-n-sign opposite` negates the CSV's printed toroidal envelope
index to obtain GLISS's phase `2*pi*(m*theta-n*zeta)`; `gliss` preserves it.
The carrier always uses GLISS's `2*pi*(M*theta-N*zeta/N_T)` convention and
is set separately by `--base`. Both the CSV sign choice and the literature
prefix convention are explicit:

| prefix convention | requested literature L | envelope CSV rows | labeled sidebands |
| --- | --- | --- | --- |
| `caption_l` | L | first L | 2L-1 |
| `terminal_l_plus_one` | L | first L+1 | 2L+1 |

`--entries` instead selects an explicit number of envelope rows and does not
attach a paper L label. No conversion from a sideband count to a paper L is
inferred. `--parity 1` and `2` select the two parity classes. Coupled parity
is unsupported by the native coefficient-normalized solve. Use separate
manifests for the two supported parities or different prefix hypotheses.

`--angular THETA ZETA` sets stiffness assembly resolution.
`--coefficient-angular THETA ZETA` independently sets the coefficient norm
grid and defaults to the assembly grid for backward compatibility.
`--reference-length` defaults to 1; changing it by a factor `a` rescales
the coefficient eigenvalue by `a**-3`, leaving exact inertia unchanged.
`--radial-cells` uses the native radial FE mesh API (at least two cells),
without resampling the equilibrium export. Omitting it retains the
equilibrium radial grid. The supported radial quadratures are `gauss5`
(default) and `cas3d_midpoint`. Coefficient normalization uses degree 1 and
the implemented constant envelope form factor `f_l=1`; no arbitrary
historical radial form function is implied. The old unimplemented
`--surfaces` option has been removed from the usage text.

The manifest is written before the first solve, with the full ordered input,
planned prefixes, all parameters and repeatable `--assumption` notes. It
records hashes for the equilibrium, CSV, converted envelope sequence, loaded
native library, runner, Python executable and tracked production source
tree, together with the Git commit and changed tracked paths. An adjacent
CMake cache, compiler identity and build flags are recorded when available;
the cache is informational and does not prove which source built a binary.
Input, table or library changes during the experiment are rejected. Failed
native runs retain the frozen manifest and error instead of dropping the
unsuccessful hypothesis. Peak RSS is converted from Darwin bytes or Linux
KiB to MiB; other platforms require an explicit unit implementation.

The output is a code-specific coefficient Ritz value, not an SI frequency
or growth rate. Compare physical assumptions, equilibrium, signs, inertia
and marginal parameters before interpreting raw normalized eigenvalue
agreement. Do not choose a reference length to fit an observed eigenvalue.
The relevant primary references are
[Schwab's CAS3D formulation](https://doi.org/10.1063/1.860656) and
[Nuehrenberg's W7-X study](https://doi.org/10.1063/1.871924).

## Matched raw-energy diagnostic

[raw_energy_comparison.json](raw_energy_comparison.json) records an independent
Bernstein-energy evaluation of the actual lowest FE192/L3 midpoint trial on
the source17/15, retained26, 144-surface W7-X export. The production replay
matches the native matrix quadratic to `3.96e-7` relatively. With two Gauss
points per cell, full-torus `2*deltaW` is `-627.226 J` for raw Bernstein energy
and `-805.852 J` for the production completed square. The magnitude gap is
28.48% of the raw value. The raw radial-quadrature drift is 0.754%; angular
refinement changes it by `1.89e-10` relatively.

The physical force residual on the actual assembly support is much larger
than the five-radius aggregate suggests: sampled volume-weighted values are
0.006033 in characteristic units and 28.31% relative to pressure force.
The receipt separates actual Cartesian curl current from reconstructed MDE
current and records force/drive contributions. No missing energy term or
exact-equilibrium counterexample was found. P0 interface terms and signed
volume corrections require further radial convergence; the retained input
and this diagnostic do not establish continuum or literature agreement.
Improve and qualify the equilibrium representation across the trial support
before repeating the same comparison.

[source_radial_comparison.json](source_radial_comparison.json) preserves the
next control: only final VMEC NS145 becomes NS289, with unchanged physical
deck, source17/15, BOOZ48 and retained26. The native solve and public export
pass. Using the frozen FE192 midpoint coefficients, force improves to
0.004182 characteristic and 19.62% pressure-relative; raw/production energies
are `-646.227/-692.575 J`, with a 7.17% relative gap. All preregistered targets
still fail (`<0.001`, `<5%`, `<1%`, respectively). The first two cells remain
included. Identical normalized-flux-chart coefficients do not imply identical
Cartesian displacements across different geometries. This single source step
establishes neither continuum convergence nor a new eigenvalue. The preserved
timeout, extended retry and matched receipts are linked by SHA256; the two
conditional CAS3D comparisons remain unexecuted.

[reconstruction_comparison.json](reconstruction_comparison.json) compares
four matching flux shells using native stored current, an independent
continuous native-wout field and the public Cartesian reconstruction.
Public/continuous-native characteristic-force ratios are 52.97, 88.57,
38.55 and 2.53. Different angular grids, nonzero native-interpolant divergence
and metric mismatches prevent source qualification. A fixed low-mode radial
factor substitution worsens every shell and is rejected; no numerical code
is changed. BOOZ's installed odd-mode staggering agrees with its source.

[retention_comparison.json](retention_comparison.json) records a matched
retained26/retained40 pair on the same NS289 source and BOOZ48 data. Both use
256x256 assembly to satisfy the unchanged retained40 anti-alias guard; the
retained26 repeat reproduces the original 192x192 values. Force improves only
0.428489%, to 0.004164 characteristic and 19.54% pressure-relative, and the
energy gap remains 7.17%. All three physical targets fail. The converter's
projection grids differ, but common position coefficients agree within
2.33e-14 m; coefficient agreement alone does not bound derivative errors.

[raw_booz_comparison.json](raw_booz_comparison.json) bypasses Cartesian
projection using all 4608 raw cylindrical BOOZ harmonics and analytic jets.
Four-shell forces agree with public retained40 within 0.0751%; divergence
is at roundoff by construction. This locates most of the discrepancy before
Cartesian projection on these shells. It does not distinguish raw radial
interpolation from source-to-BOOZ transfer or certify angular convergence.
The earlier continuous-native control interpolates contravariant B and is
not exactly solenoidal.

[native_flux_comparison.json](native_flux_comparison.json) adds the
divergence-free native flux field from R/Z, lambda and shared total-flux/
pressure jets. Its characteristic forces are
0.137940/0.058494/0.018355/0.00065216 on the four native-angle shells, with
divergence below 1.94e-16 T/m. Field differences versus stored native B track
geometric/native Jacobian differences. A substantial discrepancy thus exists
in this chosen continuous native reconstruction before BOOZ. This does not
invalidate the discrete source solve. Isolate radial continuation/staggering,
repair the demonstrated cause and repeat the unchanged whole-support gates.
No new eigenpair, literature agreement or physical qualification follows
from these diagnostics.

[axis_continuation_control.json](axis_continuation_control.json) isolates
the native m1 quotient axis closure. Although q0=2q1-q2 forces zero first-node
quotient curvature in a uniform not-a-knot spline, removing this closure
worsens force by 7.40% to 20.49% on all four shells. Positive-node positions
and pure angular jets remain bitwise unchanged. The manufactured bias is
real; this tested substitution is unsuitable as a physical repair. The
BOOZ first-half interpolation is a separate hypothesis.

[firsthalf_force_control.json](firsthalf_force_control.json) tests that
separate hypothesis using only the first m1 R/Z input geometry, one BOOZ
surface transform and an unchanged retained26 public conversion. All other
BOOZ payloads are bitwise preserved. Full-support characteristic/pressure-
relative force worsens by 12.75% to 0.004715/22.12%; both physical gates fail.
The candidate is rejected without an energy or eigenvalue solve.

The [reference source audit](reference_source_audit.md) records accessible
VMEC++ and BOOZ contracts and their build-provenance limits. Reproducible
[interpolation](symbolic_interpolation/derivation.md) and
[metric-averaging](symbolic_metric/derivation.md) FortSym derivations prove
34 and 2 exact checks respectively, with zero probes. Their source/CMake
drivers and receipts are included. No reference implementation is treated
as ground truth, and these proofs establish neither physical causality nor
the correctness of unavailable CAS3D source.

## Historical synthetic scaling result

`l139_manifest.json`: base mode (1, 1), angular grid 64 x 64, 100 radial
surfaces, two OpenMP threads on a 4-vCPU Xeon at 2.1 GHz.

| entries | labels | quotient rank | nullity | widest block | negative | lowest | certificate | peak MiB | seconds |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 3 | 5 | 995 | 0 | 10 | 0 | 0.279558 | 8.8e-06 | 85 | 7 |
| 10 | 19 | 3781 | 0 | 38 | 0 | 0.271981 | 2.8e-05 | 91 | 21 |
| 20 | 39 | 7761 | 0 | 78 | 0 | 0.261682 | 5.5e-05 | 107 | 57 |
| 35 | 69 | 13731 | 0 | 138 | 0 | 0.258997 | 9.4e-05 | 169 | 132 |
| 70 | 139 | 27661 | 0 | 278 | 0 | 0.254767 | 1.9e-04 | 523 | 382 |

Memory grows with the widest block, not with the square of the unknown
count: a dense pencil at 139 labels would need 6.1 GB per matrix. The eigenvalue
decreases monotonically as sidebands are added, and the inertia stays
stable.

The equilibrium is a fixed-boundary W7-X VMEC solution (NFP 5, beta 4.5%,
101 surfaces; `wout` sha256 `f8abe4dc...c957e`) converted with
`gliss.convert_vmec` (BOOZ_XFORM 52 x 52, 100 surfaces); the manifest
records the export hash and the GLISS commit. It is not the deck of the
published CAS3D W7-X L139 study, whose angular grid, radial form functions
and reference length are not public, so the digitized CAS3D eigenvalue is
not a target for these numbers: the table measures scaling, quotient
dimensions and certified inertia, not same-deck agreement.

## Pinned Figure 6 transcription (2026-10-21)

[V] [figure6_numbered_cells.csv](figure6_numbered_cells.csv) preserves all
204 numbered cells in Figure 6 of
[Nuehrenberg (1996)](https://doi.org/10.1063/1.871924). The hashes, method,
phase and parity conventions are in [figure6_fixture.json](figure6_fixture.json).
Separate row-wise and column-wise readings from the two private PDFs agree
with the archived native-raster audit on every cell and with the retained
L8 control on its eight rows and 15 physical sidebands. The PDFs contain
byte-identical figure rasters; this is a second reading, not independent
image evidence. The private dossier retains both explicit readings and a
reproducible reconciliation script.

[V] Section III B states `L=203` and `2L+1=407`: 204 envelope rows including
the origin produce 407 labeled sideband contributions before quotienting.
The Figure 6 caption instead relates printed cell `l` and its `l-1`
predecessors to the plotted eigenvalue. [H] The full endpoint supports
`terminal_l_plus_one`; historical small-L alignment remains conditional.
Keep `caption_l` and `terminal_l_plus_one` explicit and preserve existing
failed comparisons. Equation (3) recovers even normal parity (`parity=1`)
for the paper's results except Figure 13. The printed n sign is opposite
to GLISS, with carrier `(3,2)` and sidebands `(3+m,2-5n)`, `(3-m,2+5n)`.
No reference-length or radial-form normalization is recovered by counting
cells, and none is fitted.

[V] The latest retained full-support NS289 control still fails the frozen
force and energy-gap gates in [source_radial_comparison.json](source_radial_comparison.json).
The frozen CAS3D pair and conditional literature bands remain unexecuted
until those gates pass. A P0 interface flux belongs to the cellwise
completed-square identity; it must not be added again to raw Bernstein
volume energy. Quadrature refinement must retain all near-axis and outer
cells and split at both FE interfaces and equilibrium-spline joins.
