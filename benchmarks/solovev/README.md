# Solov'ev comparison snapshot, September 9, 2026

This is a dated snapshot of fresh GLISS two-component marginality FEEC runs
against archived independent DCON Newcomb crossings. It records GLISS source
`3b45902177e64c00e3b35908225e62a7a92218cd`, not the behavior of a future checkout.
The current operator gave a coarse sign-change interval `(1.05, 1.1)`, outside
the frozen DCON bracket `(1.039062, 1.039843)`. This is a failed comparison.

Full-volume accuracy is unqualified. Although M24 exports pass the frozen
boundary bounds and sampled orientation checks, independent reconstruction
has implicit-surface errors `7.765e-7`, `5.504e-7`, `1.731e-7` on radial meshes
64, 128, 256, exceeding the equilibrium bound `1e-7`. The mismatch cannot yet
isolate the stability operator from its interpolated equilibrium. M20 fails
the boundary bound; M28 passes the boundary and orientation checks. The M36
folded-volume diagnostic is separate and is not plotted here.

Regenerate the figure with Python, Matplotlib, NumPy and Pillow installed:

```sh
python benchmarks/solovev/plot_snapshot.py --output /tmp/solovev-snapshot
```

`qscan.csv` contains fresh GLISS raw count-only CLI output; NaNs denote
unrequested eigenpair diagnostics. `refinement_controls.csv` retains the
radial, angular and export-bandwidth controls. `dcon_reference.csv` contains
factual archived DCON outputs. The plot explicitly restricts that reference
to the marginal-range q0 values, 1.025 through 1.1. Connecting lines guide the
eye; counts are not continuously interpolated. The CSVs are a historical
evidence snapshot, not an executable physics test of the current source.

The common observable is fixed-boundary n=1 stability sign. Negative inertia
and Newcomb crossings are dimensionless; equality of individual counts is
not a mode-transfer proof. Raw GLISS L2 eigenvalues and DCON response energies
are not equated. GLISS uses parity 1, regular axis, zero normal displacement
at the plasma edge, m=0..8 with both toroidal sidebands, FEEC degree 2, and
64-by-8 angular quadrature unless a control specifies otherwise. Coordinates
are left-handed Boozer, with normalized toroidal flux s and a poloidal flip.
This is not the compressible physical-mass frequency operator or a matrix replay.

The same analytic GPEC Solov'ev family has R0=1 m, a=0.33 m, elongation 1.6,
F=1 T m and edge toroidal flux 0.7071679553076464 Wb. Pressure and iota vary
with q0 according to the frozen analytic construction. Equilibrium files and
external solver code are not bundled in this small snapshot.

Provenance:

- Research repository: `822b2bb0296f48cedd246a67cc46ad69266a5cb7`.
- Fresh GVEC source: `26a2b9dbc47e0d935c836f814c4a961da4734234`.
- DCON artifact's recorded GPEC pin: `f5595c0689c624834d720068ddc8b5a7e4028248`.
  This output-specific pin differs from the global research pins and is
  authoritative for these archived rows. DCON was not rerun.
- Source files: `gvec-stability/benchmarks/dcon/solovev/output/solovev_marginal_scan.csv`
  and `benchmarks/gvec/solovev/acceptance.json` in that research repository.
- GLISS diagnostic-only worktree patch SHA256:
  `9cfe20cfe8f1ed8c190ba36d23542a67c1aa960bf6595fc96a5b7af7be989576`.
- GNU Fortran 16.2.1, OpenBLAS, FO_JOBS=2, OMP/BLAS threads=1.

The original full fo run passed build, tests, static checks and lint. That
software check does not establish cross-code agreement. The eigenpair at the
DCON stable endpoint has native L2 value `-1.3505662135e-3`, certificate
`1.9248e-8` and backward residual `7.6472e-12`; its sign discrepancy exceeds
the solver certificate but remains conditional on equilibrium quality.
