# Changelog

## Unreleased

- Make the |m|=1 compatible FEEC trial space conforming at the magnetic axis:
  the tangential eta basis carries the regular-displacement axis factor of the
  normal component, and the third compressible unknown is the regular
  `nu = mu - (FP'/FT') sqrt(g) eta`. This removes a spurious unstable mode
  whose eigenvalue scaled with the radial mesh width; the Solov'ev n=1
  stability boundary now agrees with DCON. Eigenvalues and eigenvector
  coefficients of all |m|>=1 families change (#16).
- Refine positive-spectrum inertia brackets by bisection before inverse
  iteration; the previous midpoint shift could return a non-lowest eigenvalue
  with a certificate as wide as the coarse bracket (#17).
- Report `eigenpair_residual` as the rigorous eigenvalue distance bound
  `||K x - lambda M x||_{M^-1} / ||x||_M` on every fixed-boundary and
  marginality path, and include the inertia interval in the dense marginality
  certificate so dense and sparse certificates have the same meaning (#26).
- Add a public-source Solov'ev pipeline (GVEC 1.5.0 plus a documented export
  patch) and small toroidal and exact-displacement regression tests.

## 0.0.2 - 2026-07-16

This release supports production fixed-boundary FEEC calculations. The
TERPSICHORE replay entry points are compatibility tools for validation, not a
physical plasma-vacuum API. The complete equilibrium-to-spectrum derivative
chain, asymmetric and precomputed BOOZ_XFORM input, macOS wheels, and
production free-boundary coupling remain future work.

- Default clean single-config CMake builds to the optimized Release
  configuration, prefer threaded OpenBLAS, and retain a portable generic
  BLAS/LAPACK fallback.
- Remove the duplicate historical fixed-boundary solver and route explicit
  3-D modes, CAS3D2MN carrier envelopes, axisymmetric convenience calls, and
  compressible spectra through the shared degree-one through degree-four FEEC
  assembly.
- Canonicalize coincident carrier-envelope Fourier sidebands before assembly
  while reporting the original input-label count as provenance.
- Add reusable opaque equilibrium contexts to the C and Python interfaces.
- Return typed native status codes with caller-provided error buffers and
  caller-owned output arrays.
- Install the public `gliss.h` header, expose its wheel location through
  `gliss.get_include()`, and test it from C.
- Add context-managed fixed-boundary `StabilityProblem` objects and immutable
  certified lowest-eigenpair results for both parity classes.
- Return mass-normalized eigenvectors in documented dynamic component order
  through the C and Python interfaces.
- Add opt-in complete fixed-boundary spectra with independently recomputed
  Rayleigh quotients, backward residuals and roundoff resolutions for every
  eigenpair.
- Add five-term fixed-boundary energy decomposition through the C and Python
  interfaces, with independent total and kinetic forms, Rayleigh quotient and
  checked floating-point closure.
- Add deterministic version-1 full-spectrum result and run containers with
  exact binary64 arrays, strict archive validation and atomic writes.
- Add deterministic version-1 configuration and result files with strict
  validation, plus portable run manifests that checksum the equilibrium and
  record the Python/native software versions.
- Add schema-version queries and exact version-1 equilibrium export through
  the C and Python interfaces. Python writes replace destinations atomically;
  direct C writes refuse existing paths.
- Record the actual legacy or version-1 equilibrium schema in run manifests.
- Refuse problem manifests when the source equilibrium changed after assembly,
  and reject files modified while manifest metadata is collected.
- Evaluate variable-block residual norms with scaled compensated reductions and
  report their operation-count roundoff bound without the former naive
  dimension-linear norm estimate.
- Reconstruct radial metric and magnetic-field derivatives from primitive
  Cartesian jets and use the same checked kernel-field builder for exported
  and spline-evaluated surfaces.
- Assemble degree-one through degree-four compatible H1/L2 radial problems
  with five-point Gaussian quadrature while routing the historical P1/P0
  closure through the same angular kernel.
- Expose fixed-boundary and pressureless-pseudoplasma TERPSICHORE FORT.23/24
  compatibility solves through the hand-written C and Python interfaces,
  including inertia, eigenpair, energy, residual, and mode-overlap diagnostics.
- Convert symmetric standard VMEC equilibria through BOOZ_XFORM with explicit
  convention, mesh, provenance, corruption and force-balance gates.
- Assemble production marginality pencils directly into sparse radial blocks
  with deterministic default OpenMP ownership and no dense global temporary.
- Reject nonzero-winding equilibrium exports unless their Boozer rotating-frame
  provenance is explicit and verified.
- Require array-temporary-free optimized builds and package threaded OpenBLAS
  in the manylinux wheel.
- Publish a manylinux x86-64 wheel and source distribution.

## 0.0.1 - 2026-07-13

- Ship the hand-written ISO C and NumPy interface as the `gliss` package.
- Provide Mercier radial profiles and worst-surface objectives for GVEC/CAS3D
  equilibrium exports.
- Validate paths, quadrature sizes, native status codes, ABI compatibility, and
  package/native version agreement.
- Bundle the GLISS shared library in the Linux wheel and provide an optional
  SIMSOPT adapter.
