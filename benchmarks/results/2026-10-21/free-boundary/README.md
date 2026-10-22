# Qualified-edge Solov'ev wall refinement

[V] `measurements.json` binds these fresh Python/native solves to source
6e9a3d3, input and native-library hashes. The source archive and raw CSVs,
settings, timings, guarded-memory logs and analysis script are retained under
`/home/ert/data/handoffs/20261021_cluster_triage/GL3`. The complete CTest suite
passed 187/187 on the same cluster build. DCON is the separately qualified,
archived public GPEC reference; it was not rebuilt for these measurements.

[V] All lengths below are in the exact analytic plasma half-width. Settings:
q0=1.5, n=1, parity 1, m<=8, FEEC degree 2, 64 radial cells,
volume256x8, physical zero floor 1e-8, native true-edge M=36 fixture.
The original frozen GLISS acceptance window is [0.1514,0.1535]; the archived
DCON sign interval is [0.152422,0.1525]. Neither was changed after these runs.

| Edge mesh | Stable–unstable wall bracket | Frozen GLISS window |
| --- | --- | --- |
| [V] 64x64 | 0.15042969–0.15046875 | fail |
| [V] 128x128 | 0.15187500–0.15191406 | pass |
| [V] 256x256 | 0.15230469–0.15234375 | pass |
| [V] 256x128 | 0.15226563–0.15230469 | pass |
| [V] 128x256 | 0.15199219–0.15203125 | pass |

[V] At fixed edge128², independently changing degree 2->3, radial cells
64->128, volume theta 256->512, or mmax 8->12 preserves the measured bracket
at its 3.90625e-5 width. The square-edge 128->256 midpoint change is 0.283%.
The finest midpoint differs from the archived DCON midpoint by -0.107%; the
intervals do not overlap. There is no fitted or rescaled reference value.
The largest recorded absolute energy-closure error is 1.242e-6 in the native
doubled-energy form; the native energy API's own closure check accepted all
samples. Raw closure values are not a continuum error estimate.

[V] The scalar-potential implementation collocates constant potentials on
flat triangles (`src/boundary/scalar_potential_vacuum.f90:9`, geometry at :179,
sector-row assembly at :211). Its independent toroidal-harmonic test verifies
absolute energies and mesh convergence. [H] The observed wall shifts point
to this vacuum-edge discretization as the leading remaining error. They do
not independently certify DCON or exclude a model discrepancy. The submitted
512² pilot has a 32 GiB guard based on the measured cubic sector-storage scaling.

[V] The [FortSym derivation](../../../../gen/free_boundary/README.md) proves
20 exact edge-flux and cylindrical-wall identities without numerical probes.
[H] Full physical acceptance still requires matched equilibrium exterior
fields, pressure/tangential-field qualification and convergence of the complete
operator to its analytic limits. Passing the broad window does not close #8/A5.
