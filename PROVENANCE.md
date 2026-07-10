# Provenance

Module-to-source map.  Equation-level traceability with executable
Wolfram gates is maintained in the research repository's dossier; this
file records which published work each module implements.

| Module | Implements | Sources |
|---|---|---|
| `src/gvec_cas3d_types.f90`, `src/gvec_cas3d_netcdf.f90`, `src/gvec_cas3d_reader.f90` | GVEC CAS3D NetCDF export schema, `(s,m,n)` storage, one-period phase convention | GVEC documentation (CAS3D interface); audited exporter behavior |
| `src/netcdf_c_bindings.f90`, `src/netcdf_c_api.f90` | NetCDF C ABI access | Unidata NetCDF-C |
| `src/gvec_cas3d_adapter.f90`, `src/gvec_cas3d_reconstruction.f90`, `src/gvec_cas3d_integrals.f90` | Coordinate and Fourier conventions, pointwise reconstruction, signed-volume and derivative checks | Schwab (1993) conventions; GVEC export contract |
| `src/radial_bspline_basis.f90` | Open clamped B-spline basis, Cox-de Boor recurrence, fixed-boundary elimination | de Boor, A Practical Guide to Splines (2001) |
| `src/local_mode_model.f90`, `src/symmetric_eigensolver.f90` | Local plane-wave ideal-MHD energy prototype, dense symmetric generalized eigensolver | Bernstein et al. (1958); Freidberg, Ideal MHD (2014); LAPACK |
| `test/enzyme_drive_gradient.f90`, `cmake/EnzymeFortran.cmake` | Enzyme reverse-mode gradient of the local energy kernel | Moses and Churavy, Enzyme (NeurIPS 2020) |

Planned layers follow, in the same discipline: basis topology and mode
families per Schwab (1991) sections 3.3 and 3.4; the static two-component
operator per Schwab (1991) Aufgabe 3.2.1 and Schwab (1993); the
compressional operator and physical norm per Nuehrenberg (1996, 1999);
free boundary per Nuehrenberg (2016); subsonic flow per Nuehrenberg
(2021); the Mercier criterion per Mercier (1960) in the form of
Landreman and Jorge (2020).
