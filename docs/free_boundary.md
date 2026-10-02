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
Historical benchmark results obtained before the current angular admission
policy require a fresh run; qualitative stabilization and energy closure
alone do not establish quantitative agreement of the complete
plasma-vacuum problem.
