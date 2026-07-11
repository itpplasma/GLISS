# Provenance

Module-to-source map.  Equation-level traceability with executable
Wolfram gates is maintained in the research repository's dossier; this
file records which published work each module implements.

## Compatibility-source policy

TERPSICHORE revision `04dcf9af` is distributed under Apache-2.0 and carries
an EPFL NOTICE.  GLISS currently incorporates no TERPSICHORE source file or
translated source block.  The repository is used to establish observable
switch semantics, data conventions, formulas, and regression targets.
Implementations are normally derived into GLISS's own variables, module
boundaries, algorithms, and tests.  Literal reuse is allowed only when the
affected file is identified as Apache-2.0 material and retains the applicable
copyright, license header, NOTICE attribution, and modified-file notice.  Such
reuse must also be recorded here.

| Module | Implements | Sources |
|---|---|---|
| `src/gvec_cas3d_types.f90`, `src/gvec_cas3d_netcdf.f90`, `src/gvec_cas3d_reader.f90` | GVEC CAS3D NetCDF export schema, `(s,m,n)` storage, one-period phase convention | GVEC documentation (CAS3D interface); audited exporter behavior |
| `src/netcdf_c_bindings.f90`, `src/netcdf_c_api.f90` | NetCDF C ABI access | Unidata NetCDF-C |
| `src/gliss_capi.f90` | Public C ABI (`iso_c_binding`) entry points: library version string and integer ABI version, exported by the `gliss_c` shared library for the Python package | Fortran 2018 `iso_c_binding` |
| `src/gvec_cas3d_adapter.f90`, `src/gvec_cas3d_reconstruction.f90`, `src/gvec_cas3d_integrals.f90` | Coordinate and Fourier conventions, pointwise reconstruction, signed-volume and derivative checks | Schwab (1993) conventions; GVEC export contract |
| `src/radial_bspline_basis.f90` | Open clamped B-spline basis, Cox-de Boor recurrence, fixed-boundary elimination | de Boor, A Practical Guide to Splines (2001) |
| `src/local_mode_model.f90`, `src/symmetric_eigensolver.f90` | Local plane-wave ideal-MHD energy prototype, dense symmetric generalized eigensolver | Bernstein et al. (1958); Freidberg, Ideal MHD (2014); LAPACK |
| `src/helical_cylinder_limit.f90` | Large-aspect-ratio current-carrying stellarator vertical-mode threshold and normalized energy margin; elliptical-tokamak conducting-wall threshold | Fu, PPPL-3374 (1999), eq. 16; Fu et al., Physics of Plasmas 7 (2000), eq. 4 |
| `src/physical_constants.f90` | Shared SI vacuum permeability used by local and global physical-energy kernels | CODATA 2018 vacuum permeability |
| `src/three_component_kernel.f90` | Full compressible point energy: the established bending, shear, and field-compression components plus the product-rule divergence and physical `gamma*p*(div xi)^2` term | Nuehrenberg, Phys. Plasmas 6 (1999), eqs. 2-3 and A9-A12; Schwab (1991), eqs. 3.2.4 and 3.2.14-3.2.15 |
| `src/physical_mass_kernel.f90` | Pointwise three-component physical kinetic form in `(xi^s, eta, mu)`, including SI mass density. GLISS reverses the Appendix flux signs; its production variable is `eta=FT' xi^theta-FP' xi^zeta`, while the signed Jacobian changes the `mu` and `sigma_tilde` signs. | Nuehrenberg, Phys. Plasmas 6 (1999), eqs. A1 and A4-A8 |
| `src/perpendicular_kinetic_kernel.f90` | Pointwise physical kinetic form after projecting out the field-parallel square; this is not the TERPSICHORE reduced `C8/C10` norm | Nuehrenberg, Phys. Plasmas 6 (1999), eqs. A4-A8; TERPSICHORE revision `04dcf9af`, near-physical branch and `RHSMAT` parallel-component description |
| `src/mass_density_policy.f90`, `src/physical_mass_assembly.f90` | Positive SI radial mass-density profile and P1/P0/P0 physical-mass element in `(xi^s, eta, mu)`, with exact one-period or direct all-period phase averaging | Nuehrenberg, Phys. Plasmas 6 (1999), eqs. A4-A8 and radial discretization in section II |
| `src/dynamic_family_layout.f90`, `src/physical_mass_family_assembly.f90` | Endpoint-eliminated global P1/P0/P0 physical mass with all eta and mu cell unknowns retained | Nuehrenberg, Phys. Plasmas 6 (1999), section II; Schwab (1991), sections 3.3 and 4.2.1 |
| `test/enzyme_drive_gradient.f90`, `test/enzyme_helical_cylinder_gradient.f90`, `test/enzyme_phase_assembly_gradient.f90`, `test/enzyme_physical_mass_gradient.f90`, `test/enzyme_physical_mass_assembly_gradient.f90`, `test/enzyme_physical_mass_family_gradient.f90`, `test/enzyme_terpsichore_reduced_mass_gradient.f90`, `test/enzyme_terpsichore_reduced_mass_family_gradient.f90`, `test/enzyme_terpsichore_normalization_gradient.f90`, `test/enzyme_three_component_gradient.f90`, `cmake/EnzymeFortran.cmake` | Enzyme reverse-mode gradients of local and analytical benchmark energies, phase assembly, physical and reduced mass at their assembly levels, the TERPSICHORE normalization map, and the compressible point kernel; whole-module LLVM O2 canonicalization precedes AD | Moses and Churavy, Enzyme (NeurIPS 2020) |
| `src/mercier_diagnostic.f90`, `src/nonuniform_derivative.f90` | Mercier criterion terms from the CAS3D export, the spectral solve of the magnetic differential equation for the covariant radial field, and the kernel geometry builder with chart-metric consumption (pointwise `B_s`, `sigma_tilde`, drive chart term) | Mercier (1960); Landreman and Jorge, JPP 86 (2020), eqs. 4.16-4.20; Schwab (1991), eqs. 3.2.9, 3.2.10, 3.2.15; verified against the Suydam reduction and the gauge-shift gates in the research repository |
| `src/mode_topology.f90` | Perturbation Fourier grid, mode families, coupling and parity selection | Schwab (1991), sections 3.3-3.4; Schwab (1993) |
| `src/terpsichore_model_policy.f90` | Boundary adapter from legacy `MODELK` to independent potential-energy-model and kinetic-norm choices; the coupled integer is not retained internally | TERPSICHORE revision `04dcf9af`, `tprgl0_block_ap.f` switch descriptions and `tprgl5_block_ap.f` coefficient branches |
| `src/terpsichore_normalization.f90` | Explicit `flip=pol` maps from either raw GLISS export slopes or internal stability fluxes, plus unit-turn one-field-period Jacobian coordinates, to TERPSICHORE radian-coordinate `BJAC`, `FTP`, `FPP`, and `NI Delta(s)` | TERPSICHORE revision `04dcf9af`, `LGIKVM`, `VMTOBO`, `FOURIN`, and `RHSMAT`; audited GVEC CAS3D export and VMEC-to-TERPSICHORE converter |
| `src/terpsichore_reduced_mass.f90` | Independently expressed reduced kinetic element: `C8`, `C10`, zero `C11`, per-mode axis factors, and P1/P0 midpoint radial weights in GLISS ordering | TERPSICHORE revision `04dcf9af`, reduced `FOURIN` branch and `RHSMAT` |
| `src/terpsichore_reduced_mass_family_assembly.f90` | Caller-mapped global gather of the independently expressed reduced element; its separate fixed-boundary convenience entry eliminates normal displacement at the axis and plasma edge and retains no compressional unknown | TERPSICHORE revision `04dcf9af`, reduced `RHSMAT` block topology; Schwab (1991), sections 3.3 and 4.2.1 |
| `src/two_component_kernel.f90` | Incompressible two-component energy density: bending, shear, and compression scalar components with the interchange drive | Schwab (1991), eq. 3.2.14; Schwab (1993) |
| `src/newcomb_limit.f90` | Single-mode physical-radius cylinder assembly of the two-component functional; mode-dependent regular-axis powers and natural `|m|=1` traction; Newcomb line-bending coefficient and Suydam acceptance limit. GLISS uses `m_Freidberg=-m_GLISS` and `k=2*pi*n_GLISS/L`. The reported stiffness-only level uses an artificial identity norm and is restricted to signs and marginality. | Newcomb (1960); Suydam (1958); Freidberg, Ideal MHD (2014), eqs. 9.35-9.55; Schwab (1991) reduction |
| `src/family_assembly.f90` | Coupled Fourier-family assembly over the two-component kernel with both perturbation parities or a single symmetry class, static condensation of the tangential unknowns, and an artificial marginality norm. The default eliminates `xi^s` at a flux-coordinate axis. Physical-radius cylinder comparisons use harmonics that vanish there; the restricted `|m|=1` result is used only as a sufficient instability sign. | Schwab (1991), sections 3.3-3.5; Schwab (1993) |
| `src/eigenvalue_tracking.f90` | Certified lowest-eigenvalue tracking: Gershgorin start, inertia bisection, inverse iteration at the bracket top, inertia-window certificate | Anderson et al. (1990) inertia counting; Schwab (1991), sections 4.3-4.4; shift-tracking gate in the research repository |

Remaining layers follow the same discipline: the global compressional
operator and staggered physical-mass assembly per Nuehrenberg (1996, 1999),
free boundary per Nuehrenberg (2016), and subsonic flow per Nuehrenberg
(2021).
