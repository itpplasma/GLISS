# Source-bound algebra of resonant inverse responses

Scope: exact real rational identities and explicit domains for two locally
identified source formulas. This proves neither complete operator equivalence
nor a cause of any spectral discrepancy. There are no numerical spectra,
rootfinding, production edits or upstream repairs in this task.

## Source equations and common-chart mapping

Pinned public TERPSICHORE checkout
`/tmp/gliss-terpsichore-source-audit-20261017`, commit
`04dcf9af3ceea7ef36198508f2beb5fafd9cef20`,
`srcs/tprgl2_block_ap.f:598,602,845–855` initializes a positive delta floor,
sets `djpboo=max(djp,djpboo)`, and in the `LCURRF.ne.9` nonconstant-harmonic
branch computes

$$D=mP-nT,\quad P=\mathrm{FPP},\ T=\mathrm{FTP},\quad
R_T=\frac{D}{D^2+\delta(D+P)^2}.$$

The public QAS deck `QAS3_fx_benchmark/QAS3_fx_cur_bench_n1.data:41–42,52–53`
has LCURRF1 and DELTAJP0.04, above the source floor1e-16; for this fixed
recipe the owner source audit identifies delta0.04. The scratch proof
models that decimal as the exact real rational1/25. It does not assert
bitwise floating-point equality to that rational or audit all call histories
of the saved positive floor variable.

GLISS `export_surface_geometry.f90:726–742,790–804,807–819` uses
$d=D/\sigma$, width $w=W/\sigma$, positive
$\sigma=\max(|P|,|T|)$, and the positive-width inverse
$(1/\sigma)d/(d^2+w^2)$. Algebraically this is

$$R_G=\frac{D}{D^2+W^2},\qquad
 W=\tfrac12|mT|\,\Delta\iota.$$

The spread is a nonnegative variation from equilibrium-data radial intervals
(`primitive_kernel_geometry.f90:326–357`), not a finite-element cell width.
This is the source's finite-resolution regularization, not an exact interval
principal-value integral.

Comparison requires a common flux/Fourier chart. The source-owner mapping,
consistent with `terpsichore_normalization.f90:48–53`, is
$P_T=-N_{fp}P_G$, $T_T=-T_G$, $n_T=N_{fp}n_G$ and equal m. Therefore
$D_T=-N_{fp}D_G$ and $(D_T+P_T)=-N_{fp}(D_G+P_G)$.
For any common scaling $\kappa\ne0$ of D and P, and width scaling
$W_T=|\kappa|W_G$, both responses scale by $1/\kappa$ and the attenuation
$DR$ is invariant. The driver checks these mapping/scaling identities.
The two-pi unit-turn factor and this chart mapping do not by themselves
make unequal dimensionless attenuations equal. Raw P,T,n arrays from the
two codes must not be identified without this mapping.

## Domains, residuals and equality

Let D,P be real, delta>0, W>=0, and put
$A=D^2+\delta(D+P)^2$, $B=D^2+W^2$.
A is positive exactly unless D=P=0. B is positive exactly unless D=W=0.
Thus the common domain is

$$D\ne0\quad\text{or}\quad(P\ne0\ \text{and}\ W>0).$$

On their respective defined domains, the exact residuals are

$$DR_T-1=-\frac{\delta(D+P)^2}{A},\qquad
 DR_G-1=-\frac{W^2}{B}.$$

On the common domain,

$$\boxed{R_T-R_G=\frac{D[W^2-\delta(D+P)^2]}{AB}}.$$

Since both denominators are positive, equality holds exactly when either
D=0 with P!=0 and W>0 (both responses zero), or D!=0 and
$W^2=\delta(D+P)^2$. For fixed P,delta,W the latter nonzero equality
points are $D=-P\pm W/\sqrt\delta$. A fixed positive delta and constant W
cannot give equality on an open interval: their width-matching polynomial
has nonzero quadratic coefficient. In common-chart flux variables the
matching condition is $W^2=\delta[(m+1)P-nT]^2$.
Neither condition establishes that the actual source widths match.

For D!=0, RT equals1/D exactly iff D+P=0; RG equals1/D exactly iff W=0.
With W0 the rational formula is1/D only away from D0. GLISS's actual code
also treats the mean harmonic as zero gauge and, near arithmetic-zero d at
zero width, accepts only compatible forcing before choosing zero gauge;
otherwise it rejects the solve (`export_surface_geometry.f90:728–750`).
The rational identities do not replace those source branches.

## Limits and exceptional points

For each fixed D!=0, rational continuity gives
$R_T\to1/D$ as delta approaches0 from above and $R_G\to1/D$ as W approaches0.
The source's delta floor prevents taking an actual runtime delta all the way
to zero. Letting both regularization parameters vanish makes the difference
vanish away from resonance; varying only one generally does not.

For fixed P!=0, RT is smooth at D0 and
$R_T=D/(\delta P^2)+O(D^2)$, so tends to zero. For fixed W>0, RG is smooth
at D0 and $R_G=D/W^2+O(D^3)$, also tending to zero. The driver checks their
values and derivatives at these defined points.

If P0, instead $R_T=1/[(1+\delta)D]$ for D!=0: positive delta does not remove
this pole. If W0, RG=1/D for D!=0. Neither kernel has a continuous extension
through its exceptional zero-denominator set. Along P=cD,
$R_T=1/\{D[1+\delta(1+c)^2]\}$, showing unbounded approaches to (D,P)=(0,0).
At D0 with fixed P!=0 it is zero; therefore joint vanishing-parameter and
resonance limits must not be interchanged without further hypotheses.
Cancellation in DR does not define an undefined underlying response.

## Reproduction and proof scope

The standalone `src/CMakeLists.txt` finds installed FortSym. Reproduce with

```sh
cmake -S src -B build -G Ninja \
  -DCMAKE_Fortran_COMPILER=/opt/homebrew/bin/gfortran \
  -DCMAKE_PREFIX_PATH=/tmp/fortsym-ordered-outputs-prefix \
  -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 2
ctest --test-dir build --output-on-failure
```

Every submitted zero test uses `check_identity`: no probe fallback. Both
native FortSym and SymEngine check the27 exact residual/mapping identities,
54/54 proved with0probes. Rational limit statements use these checked
identities plus elementary continuity on the explicitly stated domains;
this is not a CAS-computed limit theorem or floating-point enclosure proof.
The receipt binds exact source and tested artifact bytes. Previously frozen
weighted-derivative proof files were not modified.
