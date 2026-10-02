# Angular grids

The magnetic differential equation (MDE) for the Pfirsch-Schlueter response
forms its forcing on an angular grid before Fourier projection. Its Jacobian
is cubic in Cartesian position derivatives. If the imported position table
has maximum absolute indices `M_eq` and `N_eq`, the grid must satisfy

```text
angular_theta > 6 * M_eq
angular_zeta  > 6 * N_eq
```

The inequalities are strict. They retain every harmonic through `3*M_eq`
and `3*N_eq` below Nyquist. Truncating the projected table after forming an
undersampled forcing cannot remove aliases that have already folded into
lower harmonics. The rule uses the declared position table, including padded
zero coefficients, rather than the requested displacement mode table.
For example, an M=8, N=0 position table admits a 64 x 8 grid; M=16, N=0
requires at least 97 x 8, and a 128 x 8 grid satisfies that requirement.
The separate displacement-product bandwidth checks also apply.

The native axisymmetric convenience routine and CLI keep their default
64 x 8 grid. For higher-bandwidth exports, the CLI accepts
`--angular NTHETA NZETA` and reports the chosen grid in its CSV result.
The public Solov'ev regeneration runner derives an admitted grid from the
retained export and displacement truncations; its M=24 default uses 256 x 8.
Python callers needing an explicit axisymmetric comparison grid can pass
the same mode family to `solve_cas3d_marginality` with `angular_theta` and
`angular_zeta`.

Python constructors generate uniform grids from `angular_theta` and
`angular_zeta`. The direct native `solve_beta_derivatives_modes` interface
requires each angle array to cover one complete uniform period, expressed
in coordinates with period one. Origin shifts, wrapping and permutations
are allowed. Repeated points, nonuniform spacing and incomplete periods are
invalid inputs.

The low-level native routine reports
`mercier_angular_alias_error = 4` when the position bandwidth exceeds the
grid. Invalid angular point sets report `mercier_invalid_input = 1`.
Higher-level equilibrium and stability interfaces can map these failures to
their input or assembly error status; increase the relevant grid size before
retrying.

This admission rule prevents folding of the polynomial cubic forcing. Ratios,
inverse metrics and other nonlinear geometry quantities require a separate
angular refinement study. Compare spectra, energy terms and equilibrium
diagnostics on successively finer admitted grids, and report the observed
changes. Passing admission or obtaining a certified discrete eigenpair does
not establish angular convergence.

The vacuum boundary mesh is independent of this volume quadrature grid.
`VacuumModel.edge_resolution` must resolve the displacement mode table and
must also be refined when studying plasma-vacuum convergence.
