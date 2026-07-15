# Changelog

## Unreleased

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

## 0.0.1 - 2026-07-13

- Ship the hand-written ISO C and NumPy interface as the `gliss` package.
- Provide Mercier radial profiles and worst-surface objectives for GVEC/CAS3D
  equilibrium exports.
- Validate paths, quadrature sizes, native status codes, ABI compatibility, and
  package/native version agreement.
- Bundle the GLISS shared library in the Linux wheel and provide an optional
  SIMSOPT adapter.
