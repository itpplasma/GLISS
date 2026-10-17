# Bounded audit of the derivative-generated tangential space

**First invalid bridge:** the published proposal's injectivity/same-count
claim is true for the physical preset powers, but false for all powers
accepted by the existing APIs. No production fix or Solovev causality is
established by this audit.

## Exact kernel and dimension, without numerical rank tests

Let the mesh partition $[0,1]$ into $c$ cells, degree $p\ge1$, and let
$U=\{N\in C^0[0,1]:N|_K\in P_p,\ N(0)=0\}$. The source knot multiplicity
makes this precisely continuous piecewise degree-$p$ polynomials.
Its dimension is $c(p+1)-(c-1)-1=pc$.
Define $L_\alpha N=(s^{-\alpha}N)'$ on cell interiors. The weight is smooth
and positive for $s>0$, and the weighted primitive is continuous at positive
breaks. Thus $L_\alpha N=0$ makes $s^{-\alpha}N=C$ with the **same** constant
on every cell. Consequently $N=C s^\alpha$.

If $C\ne0$, the first-cell polynomial satisfies
$sN'-\alpha N=0$. Writing $N=\sum_{j=1}^p c_js^j$ yields
$(j-\alpha)c_j=0$ for every $j$. Therefore necessarily
$\alpha\in\{1,\ldots,p\}$. Conversely the global polynomial $s^j$ lies in
$U$ and is a nonzero kernel vector at $\alpha=j$. Hence, exactly,

$$
\dim\ker L_\alpha=\begin{cases}1,&\alpha\in\{1,\ldots,p\},\\0,&\text{otherwise},\end{cases}
\qquad
\dim L_\alpha U=pc-\dim\ker L_\alpha.
$$

This covers all finite real powers in the accepted source input domain.
It does not infer admissible physical energy from the input validator.
The elementary counterexample $\alpha=1,N=s$ works already at degree1.
At the explicitly tested TERPSICHORE power2, $N=s^2$ is a kernel for degrees2–4.

## Source domains and concrete repair

- `compatible_two_component_problem.f90:925` and
  `compatible_three_component_problem.f90:1445` accept finite real powers;
  neither restricts them to the physical presets.
- `terpsichore_topology.f90:70–72,90` assigns arbitrary finite `qn` only to
  shifted m1, with zero elsewhere. `test_terpsichore_topology.f90:22,39–40`
  explicitly accepts qn2; positive and negative quarter powers are also tested.
- `axisymmetric_spectrum.f90:117–120` uses the physical presets
  $\alpha=0$ for m0 and $\alpha=1-m/2$ for positive m. The same formula is
  used by current fixed-boundary/marginality entry points.

For those physical presets $\alpha\le1/2$, $s^{-\alpha}N\to0$ at the axis,
so injectivity also follows immediately from the constant-kernel argument.
At $\alpha=0$, the derivative image equals the original broken $P_{p-1}$
space: integrate any broken polynomial from0 to construct its unique member
of $U$. Thus its dimension is exactly $pc$.

The smallest honest scope repair is to qualify the candidate theorem and
opt-in implementation to the physical presets, preserving legacy/custom
power behavior separately. A general-power implementation must account for
the kernel; keeping all $pc$ derivative generators at integer1..p is linearly
dependent, and dropping one gives $pc-1$ functions. A fixed-count replacement
would require an explicitly specified additional non-derivative function.
Applying $1-m/2$ blindly at m0 would give alpha1 and introduce this kernel;
the actual m0 special case must remain.

## Axis span and tie normalization

On the first cell, write a left-trace generator
$M=\sum_{j=1}^p b_j s^j$. Its derivative is
$\sum_{j=1}^p(j-\alpha)b_j s^{j-1-\alpha}$.
For physical presets every coefficient multiplier is nonzero, so this
first-cell span equals $s^{-\alpha}P_{p-1}$, the legacy local eta span.
Thus the change does not silently introduce a different leading axis power.
For m1/alpha1/2 the leading compression still requires one axis tie.
For m>=2 the source enforces finite-energy admissibility, not every stronger
smooth-Cartesian leading-coefficient identity.

Let the first normal and generating H1 function be the same raw basis
function $B(s)=c_1s+O(s^2)$, with coefficients $n_1,e_1$.
Then $x'\sim(1-\alpha)c_1n_1s^{-\alpha}$ and
$\eta\sim(1-\alpha)c_1e_1s^{-\alpha}$.
For normal cosine/tangential sine phases, the uniform cylinder requires
$\eta=-T x'/(2\pi m)$; opposite parity reverses the sign. Thus at m1

$$
 e_1=-\sigma\frac{T}{2\pi}n_1,\qquad
 \sigma=+1\text{ (normal cosine)},\quad -1\text{ (normal sine)}.
$$

The old weighted-L2 basis has first eta value1 and uses instead
$e_{\rm old}=-\sigma c_1Tn_1/(4\pi)$
(`compatible_axis_regularity.f90:87–101`). Reusing that old coefficient
factor for derivative generators is incorrect unless an explicit basis
rescaling is introduced. Degenerate meshes with no normal generator need
zero leading eta rather than a reference to a nonexistent normal coefficient.

## Wall, zero curl, and implementation boundary

The normal space has both traces removed and embeds in the left-only
eta-generating space by setting the right-edge coefficient0. Retaining the
right-edge generator permits independent eta at the wall. It imposes no
physical normal-wall condition on the generating function. Removing both
traces would reduce dimension and impose $\int_0^1\eta\,ds=0$ for the physical
presets, because the weighted primitive vanishes at both ends.

For a uniform axial cylinder, constant $T=\pi a^2B$, unit-turn theta,
nonzero m and no axial dependence,
$D=x'+2\pi m\eta/T$ for the cosine normal family and $Q=-BD$.
Thus $\eta=-Tx'/(2\pi m)$ gives zero raw curl exactly. This relation cannot
be divided by m0: there $D=x'$ and theta-independent eta does not cancel
normal divergence. At m0/alpha0 the candidate eta span remains the original
broken polynomial space; normal/tangential active ranks must remain separate,
including the zero-harmonic parity exclusions in `trial_space_topology`.

Interior derivative generators can span adjacent cells; for noninteger alpha,
a positive interior cell can see $p+1$ independent local shapes rather than
$p$ weighted polynomial shapes. This is an explicit function-space change
with shared global coefficients, not preservation of the old L2 scatter.
Mass, shear, compression, divergence, evaluation, active indices and axis ties
must use the same new basis. Same global dimension alone proves none of that.

## Strict symbolic check and limit

The scratch native CMake/CTest suite proves13 exact identities with0probes:
the unfactored derivative with native FortSym; its weight-normalized form,
explicit integer kernels, m1 coefficient, physical-preset coefficient and
cylinder cancellation with both native and SymEngine. On s>0 the real weight
is defined as exp(-alpha log s), and arbitrary N,N_s are independent jets.
Both engines prove the normalized derivative identity; the native engine also
proves the unfactored form. Preserved failed logs disclose that direct
symbolic-power normalization is undecided natively, and the SymEngine adapter
reports nonzero for some unfactored exp/log residuals. No upstream repair is
needed for this bounded audit, and no such failure is called a disproof.
Finite-element rank, sign, trace and regularity arguments above are
mathematical arguments, not conclusions of sampled matrices or CAS counts.
