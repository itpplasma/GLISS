# Free-boundary model

`StabilityProblem(vacuum=VacuumModel(...))` retains the plasma-edge normal
displacement and adds the magnetic energy of the current-free vacuum field
it drives. A scalar-potential boundary-integral solve supplies that vacuum
energy, with either decay at infinity or an ideal conducting wall.
`VacuumModel` sets a full-torus edge mesh and an optional conformal or
explicit Cartesian wall enclosing the plasma. The reported energy
decomposition includes the vacuum contribution.

The equilibrium model assumes:

- zero plasma pressure at the edge;
- continuity of the tangential equilibrium magnetic field across the edge,
  matching a current-free exterior;
- no equilibrium surface-current sheet requiring an additional surface term.

Callers must qualify these conditions for the imported equilibrium. The
constructor checks input shapes, geometry and numerical admission; it does
not solve the equilibrium exterior matching problem or certify these
physical assumptions. An eigenpair certificate concerns the assembled
matrix pencil. It does not certify equilibrium force balance, edge matching,
or convergence of the volume and vacuum discretizations.

GVEC exports may supply explicit interface data in addition to the ordinary
interior half-grid. Native ingestion accepts numeric
`boundary_schema_version=1`, scalar `boundary_s=1` with units `1`, and
`boundary_dPhi_ds` and `boundary_dchi_ds` with units `Wb`. The six
`boundary_{xhat,yhat,zhat}_{mnc,mns}` arrays use the existing `(m,n)` Fourier
grid and units `m`; stellarator symmetry permits the unused sine/cosine
complements to be omitted. Positions use the declared Boozer position frame.
The vacuum operator uses these positions and independently evaluated flux
derivatives directly. Interior geometry and profile splines retain their
half-grid contract. Malformed or unsupported declarations fail ingestion.

Without an explicit boundary declaration, positions retain the legacy
extrapolation to `s=1`, whose accuracy requires separate qualification. Flux
profiles are never silently extrapolated: a half-grid input requires independent
edge slopes via `VacuumModel(edge_flux_slopes=(dPhi_ds, dchi_ds))` (Wb), or
`plasma_vacuum_model_t%edge_flux_slopes` in Fortran. The slopes are recorded in
configuration documents and hashes. A profile grid containing `s=1` can instead
supply its endpoint derivative. Schema 1 owns its declared slopes and ignores
this schema-0 override. Schema 1 supplies no edge pressure or matched
exterior field, so its presence alone does not certify the physical assumptions
above.

The vacuum boundary mesh and volume angular quadrature are independent.
Use admitted [volume angular grids](angular_grids.md), then refine the
radial mesh, displacement mode range and vacuum edge mesh separately.
Record the wall geometry, equilibrium and solver identities for a
comparison. Parity-separated spectra require symmetry of the complete
plasma-vacuum operator, including the wall. A wall that breaks that symmetry
requires the coupled parity operator so its cross-parity stiffness is retained.

The native toroidal-harmonic oracle checks absolute vacuum energies with
and without a wall. The public Solov'ev benchmark compares a critical
conformal wall distance with DCON and studies discretization convergence.
[V] The archived 3.1952e-7 edge-transform error is due to differentiating
and extrapolating the interior half-grid flux spline, rather than an error in
the vacuum solve. Resampling the independent analytic flux on the same half-grid still gives
3.1805e-7; the remaining stored-versus-analytic spline difference is 1.47e-9.
Thus interpolation accounts for over 99% of this archived error. The native
legacy call is `primitive_boundary_geometry.f90:66`, with the derivative
formula at `radial_cubic_spline.f90:255`; schema 1 uses declared slopes at
`primitive_boundary_geometry.f90:88`. The schema-1 native-edge fixture reduces it to 1.20e-14
against the unchanged 1e-12 limit. The benchmark runner now requires this
explicit-edge qualification and uses an admitted 256 x 8 volume grid. The
[FortSym proof](https://github.com/itpplasma/GLISS/blob/main/gen/free_boundary/README.md) checks the edge-flux chain rule
and cylindrical wall limits exactly. Equilibrium exterior matching and a
fresh converged quantitative plasma-vacuum comparison remain unqualified.
Historical benchmark results obtained before the current angular admission
policy require a fresh run; qualitative stabilization and energy closure
alone do not establish quantitative agreement of the complete
plasma-vacuum problem.
