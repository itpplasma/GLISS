# Regular-quotient displacement spaces

[V] The physical perpendicular-inertia marginality path, including
`solve_axisymmetric`, uses the following space for the physical preset powers.
The low-level builder selects it with `regular_quotient=.true.`. Custom stored
powers, CAS3D coefficient/midpoint replays and the three-component compressible
solver retain their existing spaces. Coefficient vectors from these spaces
are different and must not be interchanged.

## Space and degree

Let s=r², let theta have period one, let T=FT'(s), and let sigma=+1 for a
cosine normal trial and -1 for a sine normal trial. For m>0 define

```text
w = s^(m/2),
x = xi^s = w U,
eta = w E - sigma T x/(4 pi s),
U in continuous piecewise Pp, U(1)=0,
E in broken piecewise P(p-1).
```

[V] Degree p describes the smooth quotient U. Its scalar complex is the
ordinary polynomial FEEC complex, with the wall trace removed and the axis
trace retained. There are p*c normal and p*c companion coefficients per active
Fourier trial on c cells. The generated lift is part of each normal basis
column in both stiffness and physical mass. No independent leading-axis eta
coefficient is eliminated in this representation.

[V] This changes the trial space and the meaning of the normal polynomial
coefficient. The legacy representation uses x=s^(m/2-1) N with a continuous
piecewise degree-p N that vanishes at the axis. In the new representation
N=s U would have degree p+1. The new degree-one result therefore does not
establish consistency of the legacy weighted-P1/weighted-P0 pair. It uses a
specified alternative vector space, without increasing the scalar FEEC degree
or adding undisclosed interior degrees of freedom.

For m=0, x=s U and eta=E. This gives the regular axisymmetric transverse
orders xi_r=O(r) and xi_theta=O(r). The m=0 harmonic has no division by m.

## Exact compatibility and axis regularity

[V] Direct differentiation on s>0 gives

```text
x_s + sigma 2 pi m eta/T = w (U_s + sigma 2 pi m E/T).
```

Consequently, in a uniform axial cylinder, the ordinary FEEC derivative
coefficients select E=-sigma T U_s/(2 pi m) and cancel the full perpendicular
divergence and perturbed field exactly. This works at p=1. The legacy pair
cannot do this for any nonzero boundary-zero P1 normal trial at m=1: its
cellwise derivative has an additional 1/s term, as proved by the existing
weighted-space and Lean obstruction results.

[V] In a circular chart the transverse amplitudes are A=x/(2r) and
B=2 pi r eta/T. Their helicities satisfy

```text
A + sigma B = sigma (2 pi/T) r^(m+1) E,
A - sigma B = r^(m-1) U - sigma (2 pi/T) r^(m+1) E.
```

They have the required r^(m+1) and r^(m-1) orders for a smooth Cartesian
vector. In particular eta's leading m=1 amplitude is tied to x automatically.
The independent companion is of higher radial order. The existing axis rule
s=h*u² integrates the first cell in a smooth coordinate; its point count also
covers the new quotient degree. These statements concern local regularity and
the cylinder derivative complex. They are not a global MHD stability theorem.

## Reproducible checks

[V] `gen/gliss_weighted_feec_gen.f90` checks the derivative, both parities,
helicity identities, axisymmetric derivative and the lift JVP/VJP duality with
native FortSym and SymEngine, together with the existing weighted jets.
It emits the lift and its forward/reverse leaves into the checked-in weighted
FEEC outputs. Ordinary builds need no CAS. The existing pinned Lean development
is unchanged; its theorem scope remains the legacy coefficient obstruction.

[V] `regular_quotient_feec` independently checks the two phase assembly routes,
zero-curl trials for m=1..3 and p=1..4, Cartesian physical mass, coefficient
reconstruction, and the generated lift derivatives. It also solves the m=0
uniform-cylinder fast branch using the production surface assemblers. The
independent exact value is the square of the first positive zero of Bessel J1,
found by bisection of the intrinsic Bessel function. FortSym verifies the
transformed differential equation for U(s)=J1(k sqrt(s))/sqrt(s).
For 2/4/8 cells the coarse-pair eigenvalue rates are 1.74, 3.75, 5.77 and 7.85
for p=1..4, approaching the expected 2p orders before floating-point limits.
No eigenvalue rescaling is used.

[V] The frozen 16-surface q0=1.035 Solov'ev export at 64 uniform FE cells gives
approximately -267.54, -345.53, -345.839 and -345.840 s^-2 for p=1..4.
The q0=1.045 control is positive at every degree. These are eigenvalues with
unit perpendicular inertia. They do not certify the full compressible spectrum.
The registered Solov'ev test checks both signs and their certificates for
all four degrees at 64 cells. The existing Newcomb analytic-limit tests retain
their original tolerances and independent coefficient oracles.

[H] Solov'ev asymptotic orders, regenerated-export convergence and the narrower
marginal bracket require the convergence run. Fixed-export FE refinement must
be distinguished from changing the exported equilibrium. Reduced compression
quadrature still supplies no general eigenvalue upper-bound theorem.
