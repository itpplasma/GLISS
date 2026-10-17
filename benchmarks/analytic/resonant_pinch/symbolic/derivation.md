# Resonant screw-pinch: exact energy and radial eigenvalue equations

This is an independent cylindrical derivation, checked by a reproducible
FortSym driver. It starts from the no-gravity Bernstein energy functional
quoted in `test/test_resonant_screw_pinch_energy.f90:22–28`; it does not
independently prove that general energy principle. All zero-test assertions
are exact identities, not eigenvalue samples or convergence tests.

## Physical fields, real amplitudes, and raw curl

Use physical orthonormal cylindrical components, with
$\xi=(f,ig,ih)e^{i(m\theta+kz)}$, real radial amplitudes, and

$$
B=(0,Ar,Z),\quad p=p_e+\frac{A^2(a^2-r^2)}{\mu_0},\quad
p_r=-\frac{2A^2r}{\mu_0},\quad \nabla\times B=2A e_z.
$$

The prescribed equilibrium obeys $J\times B=\nabla p$ with
$J=(\nabla\times B)/\mu_0$. Direct cylindrical cross product and curl give

$$
\xi\times B=(iC,-Zf,Arf),\qquad C=Zg-Arh,
$$
$$
Q_r=i(mA+kZ)f,\qquad Q_\theta=-kC-A(f+rf'),\qquad
Q_z=-Z(f'+f/r)+mC/r,
$$
$$
D=\nabla\cdot\xi=f'+f/r-mg/r-kh.
$$

At global resonance $k=-mA/Z$, the raw curl is **$Q=-BD$ even before
imposing perpendicularity**. In the complex contraction,

$$
-\frac{(\nabla\times B)\cdot(Q\times\xi^*)}{\mu_0}
=-\frac{2A}{\mu_0}\left[(mA+kZ)fg-Q_\theta f\right]
=p_r fD.
$$

The pressure-gradient energy term is also $p_r fD$. Therefore the exact
real-amplitude raw energy and mass densities are

$$
q=H D^2+2p_r fD,\quad H=B^2/\mu_0+\gamma p,\qquad
n=\rho(f^2+g^2+h^2).
$$

For a real, nonzero harmonic, angular/axial averaging gives the common
factor $C_{\rm ang}=\pi L$:
$2\delta W=C_{\rm ang}\int_0^a rq\,dr$ and
$N=C_{\rm ang}\int_0^a rn\,dr$.
Consequently $\lambda=2\delta W/N$ has the same radial normalization.
The complex amplitude convention instead doubles both common factors.

## Why perpendicular reduction is valid at nonzero eigenvalue

Vary the full functional
$\mathcal L_{\rm full}=r[HD^2+2p_r fD-\lambda\rho(f^2+g^2+h^2)]$.
Writing $P=HD+p_r f$, the two tangential stationarity equations are

$$
\mathcal L_g=-2(mP+\lambda\rho rg)=0,\qquad
\mathcal L_h=2r(mAP/Z-\lambda\rho h)=0.
$$

Their combination is exactly

$$
Ar\mathcal L_g+Z\mathcal L_h
=-2\lambda\rho r(Ar g+Zh).
$$

For $\lambda\ne0$, $\rho>0$, $r>0$, stationarity thus forces
$h=-Ar g/Z$. It follows that

$$
S=\frac{B^2}{Z^2}=1+\frac{A^2r^2}{Z^2},\quad
D=f'+f/r-mSg/r,\quad n=\rho(f^2+Sg^2).
$$

This reduction follows from the unrestricted resonant variation; it is
not merely an imposed trial subspace for nonzero stationary eigenvalues.
It does not apply by division at $\lambda=0$.

## Exact completion and lower bound

For $H>0$,

$$
q=H\left(D+\frac{p_r f}{H}\right)^2-\frac{p_r^2f^2}{H},
$$
$$
q+\frac{p_r^2}{H\rho}n
=H\left(D+\frac{p_r f}{H}\right)^2+\frac{p_r^2Sg^2}{H}\ge0.
$$

Hence every admissible nonzero Rayleigh quotient obeys
$\lambda\ge-\sup_{[0,a]}p_r^2/(\rho H)$.
For this constant-density equilibrium, let

$$
H_0=\frac{Z^2}{\mu_0}+\gamma\left(p_e+\frac{A^2a^2}{\mu_0}\right),
\quad H(r)=H_0+\frac{(1-\gamma)A^2r^2}{\mu_0}.
$$

The driver verifies

$$
\frac{d}{dr}\frac{p_r^2}{\rho H}
=\frac{8A^4rH_0}{\mu_0^2\rho H^2}\ge0.
$$

With $p_e,\gamma,\mu_0,\rho>0$, $Z\ne0$, and $0\le r\le a$,
$p>0$, $H>0$, and $H_0>0$, so

$$
\boxed{\lambda\ge-\frac{4A^4a^2}{\mu_0^2\rho H(a)}},\qquad
H(a)=\frac{Z^2+A^2a^2}{\mu_0}+\gamma p_e.
$$

The inequalities follow from real squares and the stated signs; the
FortSym assertions verify the exact equality decompositions and derivative,
not a separate automated inequality theorem.
The choice $D=-p_rf/H$, realized by
$g=Z^2(rf'+f-rD)/(mB^2)$, gives $q=-p_r^2 f^2/H$ exactly.
It supplies negative-energy trials where $p_rf$ is nonzero; it does not by
itself identify an isolated numerical eigenvalue or its multiplicity.

## Independent Euler–Lagrange derivation and ODE

Vary the reduced real functional

$$
\mathcal L=r[HD^2+2p_r fD-\lambda\rho(f^2+Sg^2)].
$$

The driver differentiates this expression with respect to independent
$f,f',g$, using radial jet differentiation for its total derivative.
With $P=HD+p_r f$, the exact Euler–Lagrange expressions are

$$
\mathcal L_{f'}=2rP,\qquad
\mathcal L_g=-2S(mP+\lambda\rho rg),
$$
$$
\frac d{dr}\mathcal L_{f'}-\mathcal L_f
=2r(P'-p_rD+\lambda\rho f).
$$

Eliminating $g,D$ yields the proposed system without correction:

$$
\boxed{g=-\frac{mP}{\lambda\rho r}},
$$
$$
\boxed{f'=\left[\frac1H-\frac{m^2B^2}{\lambda\rho r^2Z^2}\right]P
-\left(\frac{p_r}{H}+\frac1r\right)f},
$$
$$
\boxed{P'=\frac{p_r}{H}P-\left(\frac{p_r^2}{H}+\lambda\rho\right)f}.
$$

A scalar second-order divergence-form equation avoids division by the
possibly vanishing quantity $p_r^2/H+\lambda\rho$. Define

$$
C=\frac1H-\frac{m^2S}{\lambda\rho r^2},\quad T=C^{-1},\quad
c=\frac{p_r}{H}+\frac1r,\quad b=\frac{p_r^2}{H}+\lambda\rho.
$$

For $\lambda<0$, $C,T>0$. Since $P=T(f'+cf)$, the equivalent scalar ODE is

$$
\boxed{(rTf')'+[(rTc)'-rTc^2+rb]f=0}.
$$

Its coefficients depend on $\lambda$; this formula is not a claim of a
linear Sturm–Liouville problem with eigenvalue-independent coefficients.

## Axis and wall domain

The divided ODE is used on $0<r\le a$, with $m$ a positive integer
(the fixture uses $m=1,2,3$), constant $\rho>0$, $H>0$, $Z\ne0$,
and $\lambda<0$. At the axis, insert $f\sim F r^\alpha$ and
$P\sim P_0r^{\alpha+1}$ into the leading equations. The coefficient equations
are

$$
(\alpha+1)F+\frac{m^2P_0}{\lambda\rho}=0,\qquad
(\alpha+1)P_0+\lambda\rho F=0.
$$

Their indicial polynomial is $(\alpha+1)^2-m^2$. The regular branch has

$$
f\sim Fr^{m-1},\qquad P\sim-\frac{\lambda\rho F}{m}r^m,\qquad
 g\sim Fr^{m-1},\qquad h=O(r^m).
$$

The other branch $\alpha=-m-1$ is excluded by divergent physical mass.
Smooth Cartesian regularity requires even-power radial series after these
leading powers, with $f-g=O(r^{m+1})$, $f+g=O(r^{m-1})$, and $h=O(r^m)$.
In particular, $m=1$ requires the finite axis tie $f(0)=g(0)$.
The driver checks the indicial polynomial and regular leading coefficient
identities. It does not prove convergence of a Frobenius series.

For a rigid conducting cylindrical wall the condition is $f(a)=0$.
The variational boundary term is $2[rP\delta f]_0^a$; no additional
$P(a)=0$ condition follows at a fixed wall. The trial fixture's stronger
$f=f'=g=h=0$ at the wall is a trial construction, not extra eigenmode data.
A moving/free interface requires its interface/vacuum energy and is outside
this volume-only derivation.

## Reproduction and limits

```sh
cmake -S src -B build -G Ninja \
  -DCMAKE_Fortran_COMPILER=/opt/homebrew/bin/gfortran \
  -DCMAKE_PREFIX_PATH=/tmp/fortsym-ordered-outputs-prefix \
  -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 2
ctest --test-dir build --output-on-failure
```

The test uses FortSym `2f0571e39f36a0606b3ae150c8b3b3df92ce6750`, native and
SymEngine strict identity checks, GNU Fortran 16.1.0 with `-O3 -Wall -Wextra
-Werror=array-temporaries -fcheck=all`, and OMP/BLAS thread counts one.
The frozen receipt records the final strict counts, exact source/library
hashes, and native test results. No eigenvalue rootfinding or numerical
convergence evidence is an input to this proof.
