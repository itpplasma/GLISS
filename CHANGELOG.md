# Changelog

## Unreleased

- Solve equilibria without stellarator symmetry with a coupled parity
  operator (#10 part C). Every mode enters with both Fourier parities as one
  problem, parity class 0. A fixed-boundary problem couples when the
  reconstructed operator fails the parity admission test at an assembly
  point, whatever the file declares; `StabilityProblem.coupled`,
  `parity_classes` and `gliss_stability_problem_coupled` report it, and
  results, full spectra and their documents hold the single class 0.
  Marginality accepts `parity_class=0`, its decoupled classes refuse an
  asymmetric operator by name instead of solving it silently, and
  `solve_axisymmetric` couples an up-down asymmetric family. A shifted
  poloidal angle origin, which stores the same equilibrium with both
  parities, reproduces the union of the two symmetric class spectra on the
  Solov'ev export and on a three-dimensional QA export. Dense inertia
  certification now accounts for Sturm-count roundoff and re-sorts refined
  eigenpairs, which degenerate parity pairs otherwise break.
- Make the compatible FEEC space conforming at the axis for |m|=1 (#35). A
  smooth displacement ties the leading coefficients of `xi^s ~ s^(1/2)` and
  `eta ~ s^(-1/2)`; the independent pair had a logarithmically divergent
  compression energy, hidden by the Gauss rule in s. The first eta
  coefficient of each |m|=1 trial is now eliminated (`eta_unknowns` drops by
  one per |m|=1 trial), the axis element is integrated in `sqrt(s)`, and the
  marginality mass is the perpendicular kinetic form at unit mass density
  instead of a coefficient norm that is unbounded for regular |m|=1
  displacements. Marginality eigenvalues are now omega^2 in s^-2 for
  perpendicular inertia. On the Solov'ev q0 = 1.035 family degree 3
  converges at order 3.3 in ns instead of first order, degrees 3 and 4 agree
  to 1e-5 at ns = 128, and the degree-4 axis mode is gone; the DCON bracket
  holds at the benchmark resolution. Eigenvalues are upper bounds, so 16
  surfaces need degree 3 for this near-marginal kink. The CAS3D midpoint and
  coefficient replays keep their historical space. The operator revision is
  3.
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
- Make symmetry checks, eigenvalue convergence and bisection tolerances, and
  the marginality and TERPSICHORE zero floors relative to the pencil scale
  instead of absolute. Nonsymmetric input is rejected at every scale rather
  than silently symmetrized below unit magnitude. Marginality and
  axisymmetric results report `zero_floor`, which extends their C result
  structs (#27).
- Separate native Python tests (`pytest -m native`) from contract tests
  against fakes. Native tests load libgliss_c, fail rather than skip without
  it, and check the Solov'ev DCON stability bracket, the Mercier sign and
  Rayleigh-quotient identities. The closure registry now requires each
  evidence item to name a ctest-registered check or a collected test, and the
  installed-wheel check verifies parameter derivatives by central
  differences (#28).
- Configuration, result and run schemas are at version 5 and record the
  `discretization_revision` of the assembled operator. Older documents stay
  readable, but replaying a configuration recorded with a different operator
  (midpoint quadrature, or the pre-#16 FEEC space) raises "operator changed"
  instead of silently assembling another operator. `Equilibrium` and
  `StabilityProblem` release native memory through `weakref.finalize` when
  they are not closed (#29).
- The TERPSICHORE replay returns `negative_count = 0` with the certified
  lowest nonnegative eigenpair for stable files instead of raising, rejects
  IVAC>0 FORT.23 files passed to the fixed-boundary API with an explicit
  message instead of a misread parity, and PROVENANCE pins the fork revision
  that writes the FORT.24 schema the reader accepts (#32). Add the public
  QAS3 benchmark (`benchmarks/qas3`: STELLOPT VMEC, TERPSICHORE 1.2, GLISS).
- Solve the Pfirsch-Schlueter magnetic differential equation with the
  cell-averaged inverse D/(D^2 + w^2), where w is half the variation of
  m chi'/Phi' over the radial data cell, in the Mercier diagnostic and the
  FEEC drive. The exact inverse amplified residual resonant forcing near
  rational iota into Mercier spikes (QAS3 D_Mercier -6.9 at iota = 2/3 now
  0.144 against VMEC 0.1435) and spurious negative directions. Nonresonant
  results change by O((w/D)^2) (#33).
- Report a nonpositive exported surface metric as an error instead of
  clipping det(g) to zero, which made every Mercier flux-surface integral
  0/0 and returned NaN profiles without an error. `convert_vmec` rejects
  metric harmonics that are not positive definite, and `mercier_objective`
  raises on a non-finite profile (#34).
- Certify the radial FEEC complex for degrees 1 to 4 with manufactured
  solutions: L2 projections converge at the optimal rates p+1 (H1) and p (L2
  and the mapped derivative) on three graded meshes, monomials up to the
  space degree are reproduced to roundoff, and a corrupted derivative-map
  entry is detected (#13).
- `convert_vmec` gates on the flux-surface-averaged radial force balance
  (the solvability condition of the Pfirsch-Schlueter equation GLISS
  solves, 2e-4 for converged QAS3 and 1e-3 for li383) instead of the
  pointwise metric-B_s closure, which is 0.1-0.25 for VMEC stellarators at
  any resolution and is now reported as `force_balance_pointwise`. The
  Boozer Jacobian, |B| and covariant-current identities share one 3e-2
  tolerance because all three measure the same booz_xform-versus-geometry
  Jacobian mismatch; converged VMEC 9.0 li383 now converts at M = N = 16.
- Add `convert_boozer` for precomputed BOOZ_XFORM `boozmn` files (the export
  equals the direct `convert_vmec` result bit for bit for the same
  transform), convert asymmetric VMEC files with both parities and
  `stellarator_symmetry="False"`, record `vmec_signgs` and
  `booz_xform_source`, and refuse asymmetric equilibria in the parity-class
  operator with a named `GlissArgumentError`. A restored `boozmn` file
  carries a zeroed axis `phip`, which is now extrapolated before the
  full-grid midpoint is taken (#10, parts A and B).
- Add Sphinx documentation sources in `docs/` (`make -C docs html`) that
  separate production, compatibility and unfinished scope. Every quickstart
  example runs in the native test suite and the installed-wheel check
  against the Solov'ev fixtures (#15; hosting is left to the release
  workflow).
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
