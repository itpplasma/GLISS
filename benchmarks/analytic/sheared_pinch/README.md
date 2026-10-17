# Sheared pinch local marginal control

This benchmark checks the local Newcomb coefficients at a simple rational
surface against production surface stiffness and physical mass. It complements
the [global resonant pinch spectrum](../resonant_pinch/README.md), whose field
line resonance holds throughout the cylinder. The present result qualifies a
local marginal condition, not global finite-mode stability or an imported
equilibrium.

The independent starting point is the raw ideal-MHD energy principle of
[Bernstein et al. (1958)](https://doi.org/10.1098/rspa.1958.0023), equation 3.13,
using cylindrical curl, current and divergence. The reduced radial equation
and its Frobenius discriminant have the Newcomb form discussed by
[Mills, Hole and Dewar (2009)](https://arxiv.org/abs/0902.2612), equations
34–35 and 47–50. Our coefficients are derived from the raw energy with explicit
SI permeability; no numerical coefficient is copied from another solver.
The [strict symbolic derivation](symbolic/derivation.md) states all domains,
the cross-term boundary contribution and the singular tangent reconstruction.

## Exact equilibrium and marginal condition

For physical radius `r`, take

```text
Btheta = A r,
Bz^2 = Zs^2 + C (r^2-rs^2),
p = pe + (2 A^2+C) (a^2-r^2)/(2 mu0),
k = -m A/Zs,  F = m Btheta/r + k Bz.
```

The field satisfies radial force balance exactly and has a simple zero of
`F` at `rs`. The frozen cases use `rs=1 m`, `a=1.1 m`, `A=1 T/m`, `Zs=1 T`,
`pe=17000 Pa`, density `1 kg/m^3`, `m=2`, axial period `2 pi m`, and
`C/A^2` equal to `-9/5`, `4-4 sqrt(2)` or `-1`. Pressure is positive and its
gradient negative away from the axis in all three cases. Axial field squared
is at least `311/500 T^2` throughout the cylinder.

Removing the positive angular and `1/mu0` factors, write the reduced energy
as `integral [fN u'^2 + gN u^2] dr`. Locally
`fN=alpha (r-rs)^2+O((r-rs)^3)` and `gN=beta+O(r-rs)`, with

```text
alpha = A^2 C^2 rs^5 / [Zs^2 (Zs^2+A^2 rs^2)],
beta  = -2 A^2 (2 A^2+C) rs^3 / (Zs^2+A^2 rs^2),
Delta = 1+4 beta/alpha
      = 1-8 Zs^2 (2 A^2+C)/(C^2 rs^2).
```

The local exponents are `(-1 +/- sqrt(Delta))/2`. The exact discriminants
are `41/81`, `0` and `-7`, respectively. Positive discriminant alone does not
certify global stability. Tangential static minimization is performed only
where `F != 0`; its reconstruction can have a pole and does not supply an
admissible eigenfunction with finite physical mass.

## Production comparison and frozen acceptance

The registered native CTest `sheared_pinch_newcomb` constructs independent
normal value/derivative jets and two physical tangential amplitudes. It compares
the full production matrix with raw Bernstein polarization, checks the physical
point mass against `diag(1,0,1,1)`, verifies positive definiteness and conditioning
of the physical tangential block, and takes its Cholesky Schur complement.
The derivative-only jet has zero point value, so full point mass is singular
by construction.

Symmetric radii `rs +/- rs/16`, `rs +/- rs/32` and `rs +/- rs/64` produce
central estimates of `alpha` and `beta`. The test checks second-order central
convergence and uses the two finest levels for Richardson extrapolation.
Three profiles, both parities, both chart orientations, direct/transformed
angular routes and gamma `1`/`5/3` give 48 stencil cases and 288 point matrices.

[protocol.json](protocol.json) preserves the acceptance frozen before the first
native measurement. Independent [80-digit finite-step values](premeasurement-exact-stencils.json)
predicted the truncation before that run. The archived calculation is replayed
by `oracle-stencils.py <output.json>` with mpmath 1.3; this standalone script was
recovered from the original inline calculation after measurement and reproduces
the frozen JSON byte for byte. [native-results.json](native-results.json)
binds the candidate source and frozen library and records the measured errors.
Worst raw matrix and physical mass errors are `1.403e-15` and `2.998e-15`;
the reduced-coefficient absolute error is `4.901e-12` and the largest physical
tangent condition number `3.148e5`. The extrapolated coefficient scaled error
is `2.782e-5` and discriminant absolute error `8.088e-5`, within the unchanged
`1e-4` and `1e-3` limits. Native discriminants differ from independently
evaluated finite-step values by at most `1.953e-10`; truncation dominates.
In particular, the small negative critical-case estimate does not establish
an instability at the exactly marginal parameter.

The fresh isolated Fo audit at `144ec56` rebuilds the committed source, passes
clean lint and optimized compilation with array temporaries treated as errors,
and passes all 186 registered native tests (40.89 s). Its revision and log hash
are recorded in `native-results.json`, closing the earlier candidate archive
qualification.

Build GLISS with its native CMake configuration and run
`ctest --test-dir <build> -R '^sheared_pinch_newcomb$' --output-on-failure`.
The standalone symbolic CMake project under `symbolic/src` requires the
FortSym revision pinned in `gen/fortsym-revision.txt`; ordinary GLISS builds
do not require a CAS. This test does not close global helical thresholds,+rational-surface admissibility, equilibrium import or CAS3D literature gates.
