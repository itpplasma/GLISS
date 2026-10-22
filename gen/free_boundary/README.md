# Exact edge and cylindrical wall identities

[V] `prove_free_boundary.f90` checks 20 identities with FortSym's native exact
zero test. `check_identity` rejects unknown and numerical-probe evidence. This
is a derivation and proof program, with no production numerical kernel to emit.
Build it with an installed, compiler-compatible FortSym:

```sh
cmake -S gen/free_boundary -B build-free-boundary-proof -G Ninja \
  -DCMAKE_PREFIX_PATH=/path/to/fortsym-prefix
cmake --build build-free-boundary-proof -j8
ctest --test-dir build-free-boundary-proof --output-on-failure
```

[V] For the public Solov'ev family, write u=sqrt(psi_n),
R²=R0²+2 a R0 u cos(theta), Z=e a R0 u sin(theta)/R. Implicit
differentiation gives d(R,Z)/d(psi_n,theta)=e a² R0²/(2 R²).
The proof parametrizes the unit circle rationally; the angle derivatives in
the Jacobian are with respect to theta, not the rational parameter.
Toroidal flux is F integral(dR dZ/R), hence
Phi_psi=F e a² R0² integral(R^-3 dtheta)/2=2 pi psi0 q,
where psi0=e F a²/(2 q0 R0). With chi=-2 pi psi0 psi_n and
s=Phi/Phi_edge, the chain rule gives Phi_s=Phi_edge,
chi_s=-Phi_edge/q and chi_s/Phi_s=-1/q. The independent Decimal80 series in
`qualify_export.py` evaluates the contour integral; the symbolic proof checks
the contour Jacobian and chain rule, not the infinite-series truncation error.

[V] In a current-free cylindrical annulus a<r<b, the k=0 harmonic potential
f(r) cos(m theta), f=r^m+b^(2m) r^-m, solves Laplace's equation and f'(b)=0.
For prescribed radial field g cos(m theta) at r=a, the positive exterior
Neumann-to-Dirichlet response is

    -f(a)/f'(a) = (a/m) (b^(2m)+a^(2m))/(b^(2m)-a^(2m)).

[V] Green's identity gives the absolute vacuum energy
Wv=pi L a² g²/(2 mu0 m) times the wall factor above. The proof checks the
radial equation, wall condition, response and its strictly negative derivative
with b for m=1,2,3. The existing native scalar-potential test separately checks
absolute toroidal harmonic vacuum energies and refinement; a cylinder is not
being substituted for the finite-aspect-ratio Solov'ev torus.

[V] A uniform-current screw pinch has Btheta=A r, constant Bz,
mu0 p=A²(a²-r²), and a rigid helical translation
xi=xi0(cos(k z),sin(k z),0). It is incompressible, with
Q=curl(xi cross B)=(A-k Bz) xi0 (sin(k z),-cos(k z),0) and curl(B)=2 A ez.
The raw energy integrand Q²-curl(B) dot (Q cross xi) is therefore
xi0²(k² Bz²-A²). Its plasma energy is
Wp=pi L a² xi0²(k² Bz²-A²)/(2 mu0).

[V] In the long-wavelength limit |k a|,|k b| -> 0 at fixed q=k Bz/A,
the annular response supplies g=(A-k Bz) xi0. For w=b²/a² the total energy,
divided by pi L a² xi0² A²/(2 mu0), is

    q²-1+(1-q)²(w+1)/(w-1) = 2(1-q)(1-q w)/(w-1).

[V] Thus 0<q<1 is stabilized for 1<w<1/q and marginal at b/a=1/sqrt(q).
The no-wall expression is 2 q(q-1), the Kruskal-Shafranov limit. These are
leading long-wavelength limits; finite-k cylindrical vacuum fields require
modified Bessel functions. The proof checks these energy identities exactly.
[H] Convergence of the complete imported GLISS operator to this rigid-kink
limit remains a separate test; the symbolic threshold alone does not establish
that convergence or qualify a toroidal wall crossing.
