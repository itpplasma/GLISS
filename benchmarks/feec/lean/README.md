# Kernel-checked weighted derivative structure

This Lean development proves the finite-dimensional algebra behind the
weighted derivative proposal. It uses Lean 4.16.0 and mathlib commit
`a6276f4c6097675b1cf5ebd49b1146b735f38c02`; the committed Lake manifest pins
all transitive dependencies. There are no added axioms or admitted proofs.
Warnings are errors, and the build prints the axioms used by its main results.

## What is proved

For a **first-cell** polynomial with the left trace removed,

\[
N(s)=\sum_{j=1}^{p}c_js^j,
\qquad A_\alpha N=sN'-\alpha N,
\]

`polynomial_derivative_bridge` proves that the formal polynomial operator
has coefficient vector `((j-alpha)*c_j)`. The coefficient-to-polynomial map
is injective, and `polynomial_kernel_iff` proves that their kernels agree.
The following results are checked by Lean's kernel:

- `kernel_iff`: a zero image can have nonzero coefficients only at `j=alpha`.
- `injective_iff_nonresonant`: injectivity holds exactly when alpha differs
  from all active integer degrees 1 through p.
- `resonant_kernel`, `resonant_nullity`, `resonant_rank`: at an active degree,
  the kernel is precisely the corresponding monomial line, its dimension is
  one, and image dimension is p-1. Nonresonant image dimension is p.
- `physical_injective`, `physical_rank`: alpha <= 1/2 gives an injective map
  and image dimension p. This includes the GLISS physical presets, with the
  separate m=0 preset alpha=0.
- `degree_one_interior_no_old_companion`: on a positive interior cell, an
  affine `N=a+b*s` has normalized derivative `-alpha*a+(1-alpha)*b*s`.
  If `alpha*a` is nonzero, no old weighted-P0 companion can equal it even at
  two distinct points. The formal polynomial operator is checked separately
  in `affine_polynomial_derivative_bridge`.
- `degree_one_tent_obstruction`: the obstruction occurs at alpha=1/2 on the
  descending cell of the two-cell tent `N=s` up to h and
  `N=h*(1-s)/(1-h)` thereafter, for every 0<h<1. This tent is continuous and
  has zero normal trace at both ends. The theorem explicitly checks two
  distinct interior points of the descending branch.

## Connection and limits

On cell interiors **s>0**, real weighted differentiation satisfies

\[
s^{\alpha+1}(s^{-\alpha}N)'=sN'-\alpha N.
\]

The [FortSym derivation](../../solovev/public/symbolic_derivative/README.md)
checks this real-weight normalization strictly. Lean here proves the
polynomial operator and coefficient structure; it does not formalize the
real-power normalization or endpoint limits.

The first-cell rank is p or p-1. It is **not** the whole-mesh rank pc or
pc-1. Extending it to continuous piecewise polynomials requires the weighted
primitive's constants to agree across every positive interface; that gluing
argument remains in the [written audit](../../solovev/public/symbolic_derivative/audit.md).
An interior cell also permits a constant coefficient and therefore has a
different local coefficient model.

These proofs establish a concrete missing derivative in the legacy companion
space. They do not establish that it causes the Solov'ev spectral failure,
nor prove a complete FEEC complex, bounded commuting projection, finite-energy
domain, stability theorem, or correspondence to assembled Fortran matrices.
No generated Fortran kernel is certified by this Lean target alone.

## Reproduce

From this directory, with `elan` available:

```sh
lake exe cache get Mathlib/LinearAlgebra/FiniteDimensional.lean \
  Mathlib/LinearAlgebra/Dimension/Constructions.lean \
  Mathlib/Algebra/Polynomial/Derivative.lean Mathlib/Tactic.lean
lake build
```

Lake fetches the manifest's pinned dependencies. The library limits Lean
compilation to two threads and treats warnings, including admitted-proof
warnings, as errors. For a direct source replay after dependencies are cached:

```sh
lake env lean -j2 -DwarningAsError=true WeightedFeec.lean
```

The local `.lake/` directory is ignored. Shared local caches are an optional
execution optimization and are absent from the committed dependency graph.
