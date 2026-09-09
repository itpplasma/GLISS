# Current production Solov'ev comparison

`run_solovev_current.py inputs.csv output-directory` evaluates the current
production compatible FEEC operator through `fo exec gliss_axisymmetric`.
The input manifest has columns `q0,export` and names independently qualified
GVEC/CAS3D NetCDF exports of the GPEC analytic Solov'ev equilibrium.
The runner records input hashes, source revision, tracked patch digest, runner
hash, commands, full execution logs, and fresh CSV results. It never reads an
archived GLISS result.

Use the research repository's frozen Solov'ev equilibrium and boundary-export
acceptance limits before interpreting the result. The source manifests,
external DCON results, and precise pins live separately in
`gvec-stability/benchmarks/{gvec,dcon}/solovev`. These restricted research data
are not bundled here. An absent fixture is not a passed comparison.

The common observable is fixed-boundary stability for toroidal mode n=1:
DCON Newcomb zero crossings versus GLISS negative inertia. Zero means stable
within the admitted numerical model; a positive count means unstable. Raw
DCON response eigenvalues and GLISS perpendicular-L2 eigenvalues must not be
overlaid: their norms and meanings differ. The axisymmetric convenience call
uses the unrestricted production marginality assembly, with m=0 through
MMAX and both toroidal sidebands, regular axis, vanishing normal displacement
at the plasma edge, parity 1, and 64-by-8 angular quadrature. Record radial
surface count and FEEC degree separately. Repeat for increasing radial and
Fourier resolution; retain failed controls and unresolved drift.

The default count-only run is intended for sign comparison. `--eigenpair`
adds the native lowest-eigenpair diagnostic without changing its normalization.
This benchmark makes no claim about free-boundary stability, physical
oscillation frequencies, or TERPSICHORE matrix-replay validation.

`fo exec gliss_geometry_admission EXPORT` scans the actual primitive interpolant
between radial knots, at 256-by-8 angular points. Its status column reports
evaluation success, not physical admission. Inspect both signed-Jacobian extrema
and retain errors; a consistently signed sampled map does not prove injectivity
or volume accuracy. Its fixed toroidal grid is intended for this axisymmetric
case. Full-volume surface/profile/derivative qualification remains necessary
even when boundary and orientation checks pass.
