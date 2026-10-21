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

## Weighted radial FEEC jets

`gliss_weighted_feec_gen.f90` derives the value, first and second radial
derivatives of `x(s)=s**(-alpha)*N(s)`, with `s>0` and finite real `alpha`.
The independent inputs are the coordinate, power and three jets of `N`.
It also derives parameter JVP/VJP leaves. Production normal and tangential
weighting use the generated primal leaves, sharing an emitted include so
the existing loops can inline them and evaluate the weight once per mode.
Strict identities and the native `weighted_feec_compatibility` test cover
custom powers, integer-power kernels and the interior degree-1 obstruction.

These are local derivative kernels. B-spline topology, incidence maps, axis
ties, boundary conditions, quadrature and global assembly remain handwritten.
The [Lean coefficient proof](../benchmarks/feec/lean/README.md) supplies a
separate structural check; neither local code generation nor that bounded
proof establishes a global commuting projection or resolves Solov'ev #35.
The reproduction gate compares both the generated module and primal include.

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

The authoritative coordinate delta-W and perpendicular kinetic definitions live in
[KIN6D](https://github.com/itpplasma/kin6d), in
`gen/kin6d_symbolic_ideal.f90`; its `gen/gen_ideal_mhd_coordinates.f90` emits
`ideal_mhd_energy_generated.f90`. Both consumers compile the byte-identical
checked-in module. GLISS's previous local energy derivation has been removed.
The KIN6D commit is pinned in `ideal-mhd-revision.txt`. Only optional codegen
fetches that source; ordinary builds and wheels need neither KIN6D nor a CAS.
For a local source checkout set `FETCHCONTENT_SOURCE_DIR_GLISS_SHARED_MHD`.

The existing `gliss_energy_gen`, `gliss_check_energy` and regeneration targets
compile and run the upstream driver. They preserve the directory/report/revision
arguments and compare every generated byte. The shared generator proves the
coordinate curl projections, coefficient reconstruction, kinetic norm and a
uniform-field slab identity against the frozen-energy Hessian. These local
identities do not prove general global equivalence of the material and reduced
integrands: equilibrium, integration-by-parts and boundary contracts still apply.
The KS3 Cartesian-jet interface is unchanged.

The generated module is also included in the corresponding Enzyme checks.
Native tests retain independent energy, mass, signed-chart and derivative
oracles. Algebraic identities do not establish discretization convergence.

## Cartesian primitive geometry

`gliss_geometry_gen.f90` derives the metric, signed Jacobian and its three
coordinate derivatives directly from Cartesian dot and cross products. It also
derives magnetic components from the flux slopes and angular metric, and the
oriented second fundamental form from the normalized angular cross product.
Each family has value, JVP and VJP leaves. The generator uses the same
directory/report/revision arguments as physical mass. The production primitive
geometry builder calls these leaves and retains input and regularity guards.
The generator also emits `gliss_geometry_primal.inc`, shared by the generated
module and the primitive geometry builder. This places identical value leaves
in the caller's compilation unit so ordinary optimized builds can inline them.
The reproduction gate checks both emitted files; no link-time optimization is
required to remove these pointwise call and packing costs.

The packed arguments and results are:

| Family | Inputs | Results |
| --- | --- | --- |
| Polynomial | Nine Cartesian three-vectors: `r_s`, `r_theta`, `r_zeta`, `r_ss`, `r_s_theta`, `r_s_zeta`, `r_theta_theta`, `r_theta_zeta`, `r_zeta_zeta` | `g_ss`, `g_s_theta`, `g_s_zeta`, `g_theta_theta`, `g_theta_zeta`, `g_zeta_zeta`, `J`, `J_s`, `J_theta`, `J_zeta`, angular cross product |
| Field | `g_theta_theta`, `g_theta_zeta`, `g_zeta_zeta`, `J`, `Phi_s`, `chi_s`, `1/nfp` | `B^theta`, `B^zeta`, `B_theta`, `B_zeta`, `abs(B)` |
| Second form | Angular cross product and three angular second-derivative vectors; separate scalar orientation | `II_theta_theta`, `II_theta_zeta`, `II_zeta_zeta` |

Position jets use metres and dimensionless coordinates; flux slopes use Wb.
Tangents use the input packing and units, and cotangents use the result packing.
The field-period count is discrete: physical tangents freeze `1/nfp`. The second
form freezes `orientation=sign(J)` within either nonzero signed chart. The field
leaves require a positive angular metric and nonzero field; curvature requires
nonzero angular area. The production builder rejects invalid inputs and returns
zeroed geometry on failure.

The independent native regression uses a quadratic Cartesian map, determinant
finite differences, Cartesian magnetic vectors, three differentiation steps,
scalar-contraction reverse checks, both orientations and physical scaling.
These local derivatives do not supply a force-balanced equilibrium response or
the complete geometry-to-spectrum design chain.

## Cartesian field and spatial derivatives

`generate_cartesian_field.f90` derives Cartesian `B`, its spatial gradient,
curl and divergence from the primitive Cartesian position jets and flux
profiles. The generated field uses `B^theta=-chi_s/(nfp*J)` and
`B^zeta=-Phi_s/J`, with the existing one-period angular convention. FortSym
differentiates this definition, applies the reciprocal-basis chain rule and
proves equivalence with the covariant curl. Its 116 exact identities include
zero divergence; no numerical probe substitutes are accepted. Ordered stages
retain `J_q/J` before multiplication to avoid an unnecessary overflow.

The native `primitive_cartesian_field` module exposes point and spline-grid
evaluators. Packed inputs contain the nine Cartesian position-jet vectors
listed above, followed by `Phi_s`, `chi_s`, `Phi_ss`, `chi_ss` and `1/nfp`.
The sixteen outputs are `B(3)`, column-major `gradient(3,3)`, `curl(3)` and
`divergence`, in tesla and tesla per metre; `gradient(i,j)=dB_i/dx_j`.
Both regular signed charts and magnetic nulls are supported. Invalid inputs
return zeroed output and an error status before division by the Jacobian.

Build `gliss_cartesian_field_gen`, or use the shared `gliss_check_generated`
target and `codegen` CTest label. The directory/report/revision arguments match
the geometry generator. Ordinary builds compile the committed module without
FortSym. Independent tests cover curved charts, a true nested vacuum torus,
a finite-pressure screw pinch and Cartesian finite-difference convergence.
These are spatial derivatives. Parameter JVP/VJP, C/Python bindings and a
force-balanced equilibrium design response remain unimplemented for this API.
