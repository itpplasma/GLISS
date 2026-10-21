# Weighted radial FEEC derivation

This work connects the physical axis weights to the radial derivative spaces.
It supplies generated local kernels, an independent native compatibility
test and a bounded Lean proof. The [regular-quotient space](regular_quotient.md) now changes the physical
two-component marginality path and resolves the frozen uniform-64 signs for
all degrees. The legacy-space obstruction results below remain applicable to
the legacy representation.

## The compatibility question

The handwritten radial complex has continuous piecewise degree-p normal
coefficients N and broken degree-(p-1) companion coefficients. Its ordinary
derivative map commutes with basis evaluation. Physical amplitudes instead
use the weight w=s^(-alpha), with the separate m=0 preset alpha=0 and
alpha=1-m/2 for positive m. On cell interiors s>0,

```text
x = w N,
x_s = w (N_s - alpha N/s).
```

Multiplying both ordinary spaces by the same weight does not generally put
x_s in the weighted companion space. For an interior degree-1 coefficient
N=a+b*s, the normalized derivative is `(1-alpha)*b-alpha*a/s`; it is not
constant when alpha*a is nonzero. The first axis cell has zero left trace,
so this interior obstruction is separate from the existing axis tie.

In a uniform axial cylinder the physical tangential amplitude
`eta=-T*x_s/(2*pi*m)` cancels divergence and perturbed magnetic field for
nonzero m and no axial dependence. The native test reconstructs a boundary-zero
two-cell tent directly in Cartesian coordinates, checks that cancellation and
the production energy/mass, and exhibits the missing weighted-P0 shape.
This is a concrete approximation-space distinction, not a proof that it
causes the Solov'ev failure. The earlier derivative-generated tangential
candidate still returns positive omega squared on the frozen failing case.
An independent-node eta enrichment also fails: combined with the existing
midpoint compression rule it produces spurious negative selected energy on
the stable control, while the same interpolant has large positive full-node
energy. The next candidate needs a consistent compression projection and
quadrature together with its weighted spaces; see the retained diagnostic
in [ROADMAP.md](../../ROADMAP.md).

## Executable and mathematical evidence

`gen/gliss_weighted_feec_gen.f90` derives the weighted value and first/second
jets, along with JVP/VJP in the coordinate, power and three independent basis
jets. Its strict symbolic checks distinguish the emitted power expression
from the equivalent positive-real exponential definition. Native FortSym
does not decide some symbolic exponent-shift identities; SymEngine proves
the emitted power identities, and both engines check the defining real-log
identities. UNKNOWN is retained as a capability limit, with no numerical
probe substituted for an exact check.

Production weighting includes the emitted primal leaves in the caller, so
the original once-per-mode power evaluation and arithmetic grouping remain.
The generator also emits a public Fortran module for derivative controls.
The reproduction gate compares both outputs. No CAS is needed by ordinary
builds. Run the registered native test with:

```sh
fo test weighted_feec_compatibility
```

The [Lean development](lean/README.md) proves the polynomial derivative bridge,
coefficient kernel/rank, physical-preset injectivity and an explicit interior
tent obstruction. Its first-cell rank theorem is not a whole-mesh rank theorem.
The [written weighted-space audit](../solovev/public/symbolic_derivative/audit.md)
contains the separate interface-gluing and trace argument.

For standard FEEC problems, a subcomplex and a bounded commuting projection
are central stability hypotheses; see
[Arnold, Falk and Winther](https://arxiv.org/abs/0906.4325). This does not by
itself prove stability of the present weighted, constrained MHD operator.
Global finite-energy domains, rational surfaces, boundary constraints,
quadrature consistency and correspondence with the assembled matrices still
require qualification. A new discretization must retain the original
Solov'ev sign/bracket and analytical convergence tests without changing their
acceptance to fit its result.
