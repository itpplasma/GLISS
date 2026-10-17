# Sheared cylindrical pinch: static reduction and local threshold

This consumer-owned derivation starts from the no-gravity Bernstein energy
convention recorded in `test/test_resonant_screw_pinch_energy.f90:22–28`.
It does not prove the general energy principle. The accompanying FortSym
driver checks the displayed algebra strictly with native FortSym and
SymEngine. Identity counts do not establish inequalities, admissibility,
spectral convergence, or an eigenpair.

## Domain and physical amplitudes

Use physical orthonormal cylindrical components and real amplitudes
\(\xi=(f,ig,ih)e^{i(m\theta+kz)}\), with
\(B=(0,t(r),z(r))\). Here the symbol \(z(r)\) denotes the axial magnetic
field, whereas the phase coordinate \(z\) is physical axial distance.
Work on an open radial annulus with \(r>0\), \(\mu_0>0\), \(z(r)>0\),
\(\gamma p>0\), and \(m\ne0\). All derivatives below are with respect to
**physical radius**, not area or normalized flux. For pointwise static
elimination also require \(F\ne0\); the resonant limit is taken afterward.

Set

\[
q_\theta=m/r,\quad K^2=q_\theta^2+k^2,\quad
F=q_\theta t+kz,\quad G=q_\theta z-kt,\quad J=q_\theta z+kt,
\]
\[
C_t=zg-th,\qquad D=f'+f/r-q_\theta g-kh,
\qquad \mu_0p'=-t(t'+t/r)-zz'.
\]

The last equation is radial force balance. Cylindrical cross product and
curl give

\[
\xi\times B=(iC_t,-zf,tf),\quad
Q=\nabla\times(\xi\times B)
=\left(iFf,-kC_t-(tf)',q_\theta C_t-(rzf)'/r\right),
\]
\[
\widehat j=\nabla\times B=(0,-z',t'+t/r).
\]

The physical raw energy density is

\[
\mu_0q=|Q|^2-\widehat j\cdot(Q\times\xi^*)
+\mu_0p'fD+\gamma\mu_0pD^2.
\]

For example, the current term has the exact real form
\(z'f(Q_z-Fh)-(t'+t/r)f(Fg-Q_\theta)\). This fixes the
complex-conjugation and Fourier signs before any minimization.
For a real, nonzero harmonic, \(2\delta W=\pi L\int rq\,dr\);
the complex convention multiplies energy and physical mass by the same two.

## Exact tangent completion away from resonance

Writing \(E=f'+f/r-D\), the change of tangent variables is invertible for
\(F\ne0\):

\[
g=(tE+kC_t)/F,\qquad h=(zE-q_\theta C_t)/F.
\]

Substitution into the raw energy, using force balance, gives the exact
pointwise identity

\[
\mu_0q=K^2(C_t-C_*)^2+\gamma\mu_0pD^2
+a f'^2+2bff'+cf^2,\qquad
C_*=(Gf'+Jf/r)/K^2,
\]
\[
a=F^2/K^2,\qquad b=(k^2z^2-q_\theta^2t^2)/(rK^2),
\]
\[
c=F^2+(z^2-t^2)/r^2-2tt'/r-J^2/(r^2K^2).
\]

Thus the static minimum is attained at \(D=0,C_t=C_*\) on this domain.
The tangent quadratic matrix has trace \(K^2(B^2+\gamma\mu_0p)\) and
determinant \(\gamma\mu_0pK^2F^2\), both positive. If \(\gamma p=0\),
the parallel direction is a semidefinite gauge, so the strict-minimum
statement changes. At \(F=0\), the displayed inverse is unavailable.

Define the magnetic-unit radial coefficients
\(f_N=ra\), \(g_N=rc-(rb)'\). Exactly,

\[
r(af'^2+2bff'+cf^2)=f_Nf'^2+g_Nf^2+(rbf^2)'.
\]

Dropping the last term requires its boundary contribution to vanish:
compact support inside the annulus, or endpoint conditions that give
\([rbf^2]=0\), suffice. Wall normal displacement zero is sufficient at a
regular wall; no tangential wall constraint is inferred. An axis endpoint
would require a separate Cartesian regularity analysis, absent here.
The static radial equation is \(-(f_Nf')'+g_Nf=0\).
Physical energy coefficients are \(f_N/\mu_0,g_N/\mu_0\), since the
algebra above uses \(\mu_0q\).

## Physical-radius chart and the drive

In the right-handed chart with unit-turn poloidal and axial coordinates,
\(\theta_{\rm phys}=2\pi\theta\), \(z_{\rm phys}=L\zeta\), the fields
in the existing 13-field energy-kernel order are

\[
(T,P,T',P',I,J_c,\mathcal J,B,|\nabla r|^2,
\widehat j\!\cdot B,\mu_0p',\widetilde\sigma,\widetilde\beta)
\]
\[
=(2\pi rz,Lt,2\pi(z+rz'),Lt',Lz,2\pi rt,2\pi Lr,
\sqrt{t^2+z^2},1,z(t'+t/r)-tz',\mu_0p',0,0).
\]

Consequently the axial shear changes both **field 3** and **field 10**.
Keeping their constant-axial-field values is inconsistent with this chart.
The transformed tangent amplitudes are \(\eta=iC_t\) and
\(\nu=i(T^2+P^2)h/z\); opposite real sine/cosine conventions move the
common imaginary sign into the tangent coefficients. The generic divergence
leaf uses \(\mu=\nu+\mathcal J(P/T)\eta\), rather than \(\nu\) directly;
this change is implemented by the compressible assembly. The driver checks
the chart's normal, binormal and parallel curl projections and this full
divergence map.
The strict divergence check caught an initial draft transcription that
interchanged the two tangent-derivative groups and omitted the Jacobian on
the eta derivatives. The corrected generic expression is

\[
D=\frac{(\mathcal Jf)'}{\mathcal J}
+\frac{T(\mathcal J\eta)_\theta-P(\mathcal J\eta)_\zeta
+P\mu_\theta+T\mu_\zeta}{\mathcal J(T^2+P^2)}.
\]

The final check substitutes the stated nu shift and proves the physical
divergence with both engines. This was a proof-draft correction; the existing
production assembly contract was unchanged.

Let \(H=\widehat j\times e_r=(0,t'+t/r,z')\). The geometric drive follows
from force balance rather than from its constant-axial-field value:

\[
\mathcal A=\mu_0p''+(e_r\times B)\cdot(\nabla\times H)+|H|^2
=2t(t'+t/r)/r.
\]

The exact pointwise completed-vector identity is
\(\mu_0q=|Q+Hf|^2-\mathcal A f^2+\gamma\mu_0pD^2\).
For \(t=Ar\), this yields **\(\mathcal A=4A^2\)** even when \(z'\ne0\).
Here \(\mathcal J'=2\pi L\), with zero angular derivatives. This chart
mapping is not a claim about a different radial-coordinate fixture.

## Simple-resonance coefficients and the local discriminant

Specialize to

\[
t=Ar,\quad z^2=Z_s^2+C(r^2-r_s^2),\quad
p=p_e+\frac{(2A^2+C)(a^2-r^2)}{2\mu_0},\quad
k=-mA/Z_s.
\]

Require \(0<r_s<a\), \(A,C,m\ne0\), \(Z_s>0\), and positive \(z^2,p\)
on the chosen interval. Choose the positive axial-field branch; then
\(z'=Cr/z\), \(z''=C/z-C^2r^2/z^3\), and
\(\mu_0p'=-(2A^2+C)r\). At \(r_s\),
\(F=0\), \(F'=-mACr_s/Z_s^2\ne0\). With \(x=r-r_s\),

\[
f_N=\alpha x^2+O(x^3),\qquad g_N=\beta+O(x),
\]
\[
\alpha=\frac{A^2C^2r_s^5}{Z_s^2(Z_s^2+A^2r_s^2)}>0,\qquad
\beta=-\frac{2A^2(2A^2+C)r_s^3}{Z_s^2+A^2r_s^2},
\]
\[
\Delta=1+4\beta/\alpha
=1-\frac{8Z_s^2(2A^2+C)}{C^2r_s^2}.
\]

These coefficients are obtained by differentiating the derived \(f_N,g_N\),
not by inserting a convention-dependent published coefficient. Dividing
both by \(\mu_0\) gives their physical-energy versions and leaves
\(\Delta\) unchanged. The leading radial equation has indicial polynomial
\(\nu(\nu+1)=\beta/\alpha\), hence \(\nu_\pm=(-1\pm\sqrt\Delta)/2\).

For compactly supported test functions on one side of the resonance,

\[
\alpha x^2f'^2+\beta f^2
=\alpha(xf'+f/2)^2+(\beta+\alpha/4)f^2
-\tfrac12\alpha(xf^2)'.
\]

The derivative integrates to zero. Thus \(\Delta\ge0\) makes this leading
local form nonnegative. If \(\Delta<0\), set
\(f=x^{-1/2}\psi(\log x)\), with smooth compactly supported \(\psi\).
Its leading energy is
\(\int[\alpha\psi_t^2+(\beta+\alpha/4)\psi^2]dt\).
Long smooth plateaus make it negative; translating the support toward
\(x=0\) makes higher-order smooth-coefficient corrections arbitrarily small.
This is the local threshold argument, not a global stability theorem.
The critical case alone does not decide higher-order or global terms.
Equivalently, the derived threshold is

\[
\frac{r_sZ_s^2}{8}\left(\frac{z'_s}{Z_s}\right)^2
+\mu_0p'_s\ge0,
\]

because \(\beta+\alpha/4\) is its left side times the positive factor
\(2A^2r_s^2/(Z_s^2+A^2r_s^2)\).
For \(A=Z_s=r_s=1\), the exact choices
\(C=-9/5,4-4\sqrt2,-1\) give \(\Delta=41/81,0,-7\), respectively.
These are local discriminant labels, not finite-mode spectra.

## Physical mass and the singular inverse

Physical mass uses \(\rho(f^2+g^2+h^2)\), with \(\rho>0\), not just
\(\rho f^2\). The exact identity

\[
(t^2+z^2)(g^2+h^2)=(zg-th)^2+(tg+zh)^2
\]

shows why static elimination need not produce a finite-mass displacement.
For complex Frobenius branches the corresponding identity uses absolute
squares on both sides; it follows by applying the real identity separately
to the real and imaginary parts.
For a generic nonconstant Frobenius branch \(f\sim x^\nu\),
\(C_*\sim (G_s/K_s^2)\nu x^{\nu-1}\), with \(G_s\ne0\).
Its perpendicular tangent mass is locally integrable only if
\(\operatorname{Re}\nu>1/2\), although radial mass needs only
\(\operatorname{Re}\nu>-1/2\). A zero leading derivative, especially
\(\nu=0\), needs its next term checked separately. The \(1/F\) reconstruction
can also create singular parallel mass. None of the elimination identities
therefore establishes a physical-mass eigenfunction at the pole. They
qualify a local static candidate for subsequent independent behavioral
testing; no generated production kernel or global eigenpair is changed.

## Reproduction

The initial diagnostic also substituted the minimum directly into the fully
expanded raw inverse. Native FortSym reported NONZERO for that redundant
residual while SymEngine proved ZERO. The exact source and log are retained
in the execution receipt. Read-only diagnosis traced this verdict to a
2000-term polynomial-normalization refusal followed by an incorrect NONZERO
fallback; the diagnosis artifacts are also bound by the receipt.
The accepted proof first checks the full raw
completion directly with both engines, then instantiates its proved square
completion at the minimum. The disputed native verdict is not endorsed, and
this evidence does not claim that generic verifier limitation has been fixed.

With the repository-pinned repaired FortSym installed:

```sh
cmake -S benchmarks/analytic/sheared_pinch/symbolic/src \
  -B /tmp/gliss-sheared-pinch-proof -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH=/path/to/fortsym-prefix
cmake --build /tmp/gliss-sheared-pinch-proof --parallel 2
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  ctest --test-dir /tmp/gliss-sheared-pinch-proof --output-on-failure
```
