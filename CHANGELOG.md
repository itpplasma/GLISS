# Changelog

## Unreleased

- Add reusable opaque equilibrium contexts to the C and Python interfaces.
- Return typed native status codes with caller-provided error buffers and
  caller-owned output arrays.
- Install the public `gliss.h` header, expose its wheel location through
  `gliss.get_include()`, and test it from C and C++.
- Add context-managed fixed-boundary `StabilityProblem` objects and immutable
  certified lowest-eigenpair results for both parity classes.
- Return mass-normalized eigenvectors in documented dynamic component order
  through the C and Python interfaces.

## 0.0.1 - 2026-07-13

- Ship the hand-written ISO C and NumPy interface as the `gliss` package.
- Provide Mercier radial profiles and worst-surface objectives for GVEC/CAS3D
  equilibrium exports.
- Validate paths, quadrature sizes, native status codes, ABI compatibility, and
  package/native version agreement.
- Bundle the GLISS shared library in the Linux wheel and provide an optional
  SIMSOPT adapter.
