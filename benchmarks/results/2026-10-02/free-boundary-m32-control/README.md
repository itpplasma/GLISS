# Solov'ev M32 boundary control

This generated q0=1.5 equilibrium remains a failed control for physical
free-boundary acceptance in issue #8. The corrected analytic generator and
M32 export resolve the total-flux and geometry failures found in the M16
attempt, but the native extrapolated edge flux derivatives fail the frozen
boundary limits. No quantitative wall-window acceptance is claimed.

`manifest.json` pins the generated files, public generator sources,
equilibrium resolution and independently chosen wall window. The export
contains 48 interior half-grid surfaces, poloidal modes 0 through 32,
one field period and winding -1. It was generated with patched GVEC 1.5.0
from the public analytic GPEC Solov'ev family. These generated files are
released under the repository's MIT license.

`qualification.json` reports independently reconstructed analytic errors.
Its flux reference uses an 80-digit positive hypergeometric series.
Cylindrical radius is `hypot(xhat,yhat)`, including the Boozer frame shift.
The total-flux error is 9.86e-17 Wb, the generator's boundary implicit error
is 7.82e-10, and the exported half-grid implicit error is 4.17e-9. The
source-equivalent extrapolated edge has implicit error 2.27e-8 and maximum
parameter-matched shape distance 4.51e-9 m.
An independent corruption control shifts every exported constant Cartesian
x coefficient by 0.01 m. The implicit error rises to 0.07955, so the
analytic geometry oracle detects the defect; its result is recorded in
`corrupted_geometry_control.json`.

`native_edge.txt` and `native_edge.json` report the actual native quantities
used by the vacuum operator. The edge toroidal-flux derivative differs
from the analytic value by 1.1983e-14 Wb, exceeding 1e-14 Wb. The native
`chi'/Phi'` differs from the analytic edge transform by 3.1952e-7, exceeding
1e-12. The stored-iota spline is a separate diagnostic; the native field
uses the differentiated flux profiles. The reconstructed edge pressure is
4.77e-5 Pa. Exterior equilibrium-field matching remains unqualified.

The frozen boundary contract belongs to explicit schema-1 edge data at
s=1, with modes 0 through 36 and declared edge flux derivatives. This
ordinary half-grid export lacks that schema. Applying the corresponding
numeric derivative limits to the native extrapolation reveals the
remaining boundary requirement; this record does not substitute a stored
iota value or relax a limit. A faithful explicit edge export and native
ingestion, followed by qualification of the equilibrium boundary model,
are still needed before the complete plasma-vacuum comparison can pass.

The independent wall bounds remain 0.1514 and 0.1535 plasma half-widths.
The exact analytic half-width is 0.3526573415939913 m; the archived older
scan used the approximate 0.35245 m scale. The bounds were frozen before
fresh GLISS wall results were inspected and were preserved when the
qualification failure required M32 and a 256-by-8 volume grid. No wall
solve was used to select, widen or validate these bounds in this record.

From the repository root, reproduce the independent measurements with
Python 3.11 or newer and NumPy, SciPy and netCDF4 installed:

```sh
python benchmarks/solovev/free_boundary/qualify_export.py \
  benchmarks/results/2026-10-02/free-boundary-m32-control/parameter.toml \
  benchmarks/results/2026-10-02/free-boundary-m32-control/export.nc \
  --q0 1.5 --theta-points 32768
fo exec benchmark_boundary_profiles \
  benchmarks/results/2026-10-02/free-boundary-m32-control/export.nc
```

The native helper is built with `BUILD_TESTING=ON`. It prints diagnostics
and makes no acceptance decision. The analytic qualifier reports errors
without solving a stability spectrum or claiming exterior matching.
