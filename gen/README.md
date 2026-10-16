# FortSym kernel generation

Ordinary builds and wheels compile the checked-in `src/generated/` modules and
do not need FortSym or another CAS. Generation is an explicit developer option.
Install the FortSym commit recorded in `fortsym-revision.txt` with the same
Fortran compiler used for GLISS, then configure its installed CMake package:

```sh
cmake -S . -B build-codegen -G Ninja -DGLISS_ENABLE_CODEGEN=ON \
    -DCMAKE_PREFIX_PATH=/path/to/fortsym/install
cmake --build build-codegen --target gliss_check_generated
ctest --test-dir build-codegen -L codegen --output-on-failure
```

`gliss_regenerate_kernels` explicitly replaces the committed sources after each
generator succeeds. `gliss_check_generated` and the `codegen` CTest tests instead
regenerate into the build tree and compare every byte. A failed or undecided
strict identity, missing output/report, or source difference fails the check.
Identity reports remain under `build-codegen/generated-kernels/`. The automatic
generated-kernels workflow uses FortSym's pinned Nix dependency/compiler closure,
builds FortSym from the recorded commit, and retains provenance and proof reports.
The installed package exports no compiler/revision metadata: developers must
select the matching installation; compilation and exact reproduction check it.

## Physical mass

`gliss_physical_mass_gen.f90` defines the pointwise physical kinetic energy in
GLISS's `(xi^s, eta, mu)` variables from the orthonormal displacement components
of Nuehrenberg, *Physics of Plasmas* 6 (1999), Appendix A1 and A4–A8. The flux,
eta and signed-Jacobian conventions are recorded in `PROVENANCE.md`. This is a
source-authored derivation; no personal Mathematica results or Wolfram oracle
are inputs to the generator.

The energy's displacement Hessian supplies the symmetric physical mass matrix.
Native FortSym differentiation supplies its thirteen argument partials, from
which the generator forms forward and reverse contractions. Each leaf emits
its value and derivative product together to share subexpressions. The emitted
Fortran module is committed in `src/generated/`; ordinary builds require no CAS.
The existing `physical_mass_matrix` and `physical_mass_energy` interfaces call
these generated leaves.

The additive Fortran JVP/VJP interfaces accept thirteen packed inputs in this
order: `flux_t_slope`, `flux_p_slope`, `current_i`, `current_j`, `signed_sqrtg`,
`bmag`, `grad_s2`, `signed_sigma_tilde`, `beta_tilde`, `density_kg_m3`, `xi_s`,
`eta`, `mu`. Tangents use the same order and units; a VJP returns the cotangent
times the energy gradient. The domain is `signed_sqrtg /= 0`, `bmag > 0`,
`grad_s2 > 0`, `flux_t_slope**2 + flux_p_slope**2 > 0`, and positive density.
Both nonzero Jacobian orientation branches are supported. These pointwise
leaves assume admissible inputs, as did the original mass kernel.

Build the optional `gliss_physical_mass_gen` executable against the FortSym
revision in `fortsym-revision.txt`, using the same compiler as that FortSym
installation. The generator accepts an existing output directory and an
optional identity-report path and FortSym revision:

```sh
gliss_physical_mass_gen src/generated physical-mass-identities.json \
    "$(cat gen/fortsym-revision.txt)"
```

Before output is written, strict symbolic checks verify symmetry, the
energy/matrix quadratic identity, displacement derivatives, density scaling,
JVP/VJP duality and every simplification accepted for emission. All checks
must be decided symbolically; probe-only evidence fails the generator. The
identity report records the proof/probe distinction. This algebraic evidence
does not certify floating-point rounding or physical/discretization convergence.

`check_generated.cmake` runs the generator and compares its result to the
committed source. The native `physical_mass_kernel`
test independently reconstructs contravariant displacement and metric kinetic
energy for two positive-definite metrics and both orientations. It also checks
all thirteen coordinate partials and mixed directions over three finite
difference steps, reverse/forward duality and the Jacobian derivative sign.

The generated contractions use native scalar `diff`. FortSym's contracted
products currently need separate qualification for `abs` and opaque
multivariate heads; this generator does not claim support for those broader
APIs. It also does not supply equilibrium, assembled spectral or vacuum design
derivatives.

## Fixed-geometry pressure reverse leaves

`gliss_pressure_adjoint_gen.f90` derives pressure-energy contractions and the
cubic-spline evaluation, right-hand-side and endpoint pullbacks with native
FortSym differentiation. Each derived expression and accepted simplification
must pass a strict symbolic identity. Its command-line interface is:

```sh
gliss_pressure_adjoint_gen src/generated/gliss_pressure_adjoint_generated.f90 \
    "$(cat gen/fortsym-revision.txt)" pressure-adjoint-identities.json
```

The native pressure-surface and pressure-stiffness tests retain independent
cardinal-gradient, finite-difference and transpose-duality checks. These leaves
support a fixed-geometry pressure adjoint; they do not establish force-balanced
equilibrium or vacuum design derivatives.

## Local energy kernels

`gliss_energy_gen.f90` defines the two- and three-component energy principle
and perpendicular kinetic energy from the cited orthonormal components. It
derives coefficient contractions, values, argument gradients and JVP/VJP
leaves, retaining strict identities for each accepted simplification. The
generator uses the same directory/report/revision arguments as physical mass:

```sh
gliss_energy_gen src/generated energy-identities.json \
    "$(cat gen/fortsym-revision.txt)"
```

The generated module is also included in the corresponding Enzyme checks.
Independent native energy and derivative tests remain necessary; algebraic
identities alone do not establish discretization or cross-code agreement.
