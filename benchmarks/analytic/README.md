# Current analytical benchmark runners

Run `bash benchmarks/analytic/run.sh /absolute/output/path` from a configured
GLISS checkout with `fo`, NumPy, SciPy and Matplotlib. Targets are enabled with
BUILD_TESTING. Outputs are fresh CSVs, PDF/PNG figures, numerical metrics,
source copies, and toolchain/source provenance. Neither archived CSV values
nor fitted references enter the calculation.

The default theta-pinch runner prescribes exact straight-cylinder kernel fields
and calls the **current production** radial FEEC basis, axis factor, split
quadrature, compressible surface assembly and physical mass assembly. The
benchmark's small radial gather loop contains no copied physical-energy
kernel. It tests both parity classes, degrees 1–4 and 8, 16, 32 intervals,
with 8×4 angular points (exact for the constant geometry and m=3,n=1 mode).
All matrices are freshly assembled. Full spectra, scale-normalized eigenpair
residuals and maximum mass-orthogonality defects are retained. This does not
test export reconstruction, public API admission or lowest-mode certification.

In the right-handed chart (s=r²/a², θ/2π, z/L), the exact nonzero fields are
FT′=πa²B, I=LB, sqrt(g)=πa²L, |B|=B and |grad s|²=4s/a². FP′,
flux curvatures, current, pressure gradient, chart shifts and drive vanish.
These follow directly from the cylindrical Euclidean metric and uniform B.
Normal displacement carries the same s^(1/2) axis factor as the current
production m=3 configuration; the H1 trace vanishes at both ends.

For radius a=0.5 m, period L=6π m, B=1 T, p=100 Pa, density 2 kg/m³,
γ=5/3 and (m,n)=(3,1), define vA²=B²/(μ0ρ), cs²=γp/ρ, k=2πn/L,
and q²=k²+(j′₃,l/a)². The two magnetosonic branches solve
ω⁴−(vA²+cs²)q²ω²+vA²cs²k²q²=0. The slow accumulation point is
k²vA²cs²/(vA²+cs²). The Bessel derivative roots enforce zero normal
wall displacement. These formulas follow the uniform ideal-MHD equations
and match the derivation in research `derivations/cylinder_compressional_spectrum.wl`.
We use μ0=4π×10⁻⁷, the explicit convention in `physical_constants.f90`.
Fast branches are the first two eigenvalues above 1.01 times the independent
Alfvén point; the slow diagnostic is the actual minimum, not a nearest match.

The exact-cylinder fast branches converge to the independent values. The
p=1 fast errors show approximately second-order radial convergence and p=2
approximately fourth order; higher-degree errors reach floating-point limits
before reliable rates can be inferred from all three meshes. The CSV records
measured rates separately; no new acceptance tolerance is fitted to these
runs. The narrow slow continuum is ill-conditioned in the full dense spectrum:
at p=4, 32 intervals, its minimum drifts about 2.2×10⁻⁴ relative below the cusp.
Small normalized residuals or mass orthogonality defects alone cannot prove
accuracy of these small eigenvalues in a pencil with a large spectral range.
This is explicitly retained, not counted as optimal FEEC convergence.
Screw-pinch/Newcomb, Suydam, scaling and smooth manufactured derivative/trace
checks remain separate work.

A second, explicitly labeled **invalid legacy-fixture comparison** is retained
as a diagnostic. The optional second argument to `benchmark_theta_pinch`
selects that path. It uses `cylinder_fixture(poloidal_scale=0)` and the full
primitive production assembly. The old fixture supplies cylinder metrics but
position harmonics xhat=(R+r cosθ)cosζ, yhat=(R+r cosθ)sinζ, zhat=r sinθ,
with winding=1. The production primitive path applies the contractual winding
rotation, yielding x=(R+r cosθ)cos(2ζ), y=(R+r cosθ)sin(2ζ). Thus it reconstructs
a double-covered torus with axis period 12π m, not the old supplied cylinder
metric with L=6π m. Its slow minimum is approximately one quarter of the
straight-cylinder cusp, consistent with the doubled length; finite aspect
ratio also remains. This is an invalid physical reference comparison, not
established evidence of a defect in the production energy kernel. It explains
why simply replaying the old fixture cannot qualify the new assembly.

The homogeneous plot tests only `assemble_local_mode` and its 3×3 eigensolve
against the independent plane-wave Alfvén/magnetosonic dispersion at 41 angles.
Its error gate is 10⁻¹² relative to vA². The helical plot checks the analytical
limit function against [Fu, PPPL-3368 (1999)](https://bp-pub.pppl.gov/pub_report/2000/PPPL-3368.pdf),
eqs. 15–16, including both signs around
f=0.4 for elongation two (absolute tolerance 10⁻¹⁴). It is a formula check,
not an independent global solver benchmark. Local success does not resolve
the theta-pinch discrepancies.

Research context read at gvec-stability `822b2bb0296f48cedd246a67cc46ad69266a5cb7`:
README, ROADMAP, validation/pins.env and cylinder_physical acceptance manifest.
That manifest pins historical GLISS `a2cfece6748bd03ea943aba75a3f14eeb9442bbb`;
its old executable is absent from current GLISS. This runner uses the checkout
recorded in its own output provenance instead. Plots use the Okabe–Ito palette
with marker and line-style redundancy; parity curves can coincide.

On the fresh exact-cylinder run, the largest scale-normalized residual is
about 1.8×10⁻¹³ and the largest mass-orthogonality defect about 1.2×10⁻¹⁴.
Residual normalization is ||Kv−ω²Mv||₂ / ((||K||F+|ω²| ||M||F)||v||₂).
These describe the computed finite pencil and do not bound relative error
in the slow cluster. Fast-mode convergence orders are recorded in
`theta_rates.json`, including floating-point plateaus and adverse rates.

## Slow-spectrum conditioning follow-up

`conditioning.py` diagnoses the same assembled pencils independently with
SciPy LAPACK and 50-digit mpmath. Run the matrix dump and script commands in
its module docstring; generated matrices remain outside the repository.
Each JSON case records the raw matrix SHA256 and library versions. The
original frozen benchmark output remains unchanged.

At 50 digits, symmetrized assembled p1,n8 and p4,n8 minima are respectively
9.25732044220955416856 and 9.25732041138864474023 s⁻², both above the cusp.
For p4,n16 the minimum is 9.25732041077774179891 s⁻², also above it. These
high-precision solves operate on rounded assembled matrices; they isolate
spectral roundoff, not assembly accuracy. The first radial analytical slow
branch is 9.25732345592264179252 s⁻²; higher radial branches approach the cusp
from above, so the finite-space minimum is not the first radial branch.

The upper-triangle dense helper used by this benchmark loses relative slow
accuracy. For this component ordering the lower triangle is substantially
better, with or without diagonal scaling. Reversing both rows and columns
swaps the benefit, even after symmetrizing K and M. Thus changing U to L alone
is not a general remedy. The existing allocated helper's equilibration option
also switches triangles; attributing the benefit only to mass scaling is wrong.
This diagnosis does not imply the public certified solver has the same defect.

At p4,n32, mass condition is about 3.91×10⁴ and the whitened spectral norm
is about 2.34×10¹¹ times the cusp. Machine epsilon times that norm is about
4.82×10⁻⁴ s⁻²: global backward accuracy can coexist with poor relative accuracy
near the cusp. The diagnostic also reports ||L⁻¹(Kv−λMv)||₂, M=LLᵀ, for a
mass-normalized vector. This bounds distance to some discrete eigenvalue;
it does not identify a member of the narrow cluster. No production solver
change or fitted acceptance tolerance is included.
