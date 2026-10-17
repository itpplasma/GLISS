# Bounded QAS3 current-regularization source audit

Date: 2026-10-17. Read-only audit; no code builds or numerical physics runs.
GLISS source anchor: 3ad9b1b5614ff9e4ff2fdd50de1f2712b6c27c4c.
TERPSICHORE anchor: 04dcf9af3ceea7ef36198508f2beb5fafd9cef20.
Source bindings and line anchors are frozen in receipt.json.

## Finding

The checked QAS3 recipe compares different equilibrium-current regularizations.
This is an operator mismatch, not merely a kinetic normalization difference.
It is a plausible explanation to isolate, not a demonstrated cause of the
historical count difference, and establishes no defect in either implementation.
The strongest tractable seam is the radial force-balance MDE; compressional
elimination, complete axis-domain equivalence and the whole operator remain unaudited.

## Common continuum equation and exact Fourier inverse

Let e_i=partial_i r, J=e_s dot(e_theta cross e_zeta),
B=(P e_theta+T e_zeta)/J, and B_i=(beta,G(s),I(s)) in a radian chart.
With jhat=curl B, jhat^theta=(beta_zeta-I')/J and
jhat^zeta=(G'-beta_theta)/J. Therefore

  (jhat cross B)_s=(P beta_theta+T beta_zeta-T I'-P G')/J.

Radial force balance jhat cross B=mu0 p' grad s consequently gives

  (P partial_theta+T partial_zeta) beta=mu0 p' J+T I'+P G'.

For phase m theta-n zeta and nonconstant forcing F_c cos+F_s sin,
D=mP-nT; the exact periodic solution has beta_s=F_c/D,
beta_c=-F_s/D when D!=0. GLISS unit-turn angles introduce the explicit
2*pi divisor at export_surface_geometry.f90:738/740/751/753.
At D=0 a classical periodic solution requires zero forcing in that harmonic;
choosing beta=0 for incompatible forcing does not solve the ideal equation.
Mean-force balance and the zero-mode gauge are separate contracts.

## Actual pinned TERPSICHORE response

BOPHYS tprgl2_block_ap.f:598/602 sets delta=max(DELTAJP,1e-16).
At 848-855 it evaluates, for LCURRF!=9 and every nonconstant equilibrium mode,

  R_T(D)=D/[D^2+delta(D+P)^2], P=FPP, T=FTP,
  FSIGBS=(PP*FBJAC+FGPARP)*R_T,
  FPARJP=(m CI-n CJ)*PP*FBJAC/(CJ*P-CI*T)*R_T.

The VMEC text converter sets hot-pressure coefficients PBPC/PPPC to zero and
the mean SIGC/TAUC to one (vmecv92terps.f:163-174). This isolates the
isotropic radial-force-balance term; this is not evidence for arbitrary hot plasmas.
The public QAS3 deck has LCURRF=1 (41-42), DELTAJP=0.04 (52-53), MODELK=0 (55).
GLISS run_case.sh:31-32 uses that template; make_deck.py changes tables,
IVAC/shift/resolution only and preserves DELTAJP and LCURRF.
These are current recipe defaults, not recovered hashes of historical executed decks.

The response enters potential energy: FSIGBS/FPARJP reconstruct SIGBS/PARJP
(tprgl2:901/906); MODELK=0 sets A1/A4/A5 from them (tprgl5:664-674);
A1/A4/A5 enter C1/C4/C5 (750-775,824-832), hence C9 (853-859) and
LHSMAT bulk couplings (340-395). It is not confined to the mass matrix.
There is also separate first-cell SIGBS extrapolation (tprgl2:930-932);
it is not accounted for by the regularization comparison alone.

## Actual current GLISS response

For S=max(abs(P),abs(T)), d=D/S, its positive-width branch is
r_G=d/(d^2+w^2), with physical R_G=r_G/S=D/(D^2+W^2), W=S*w.
At export_surface_geometry:726-740/790-803,

  W=(1/2)*abs(m T)*Delta_iota.

primitive_kernel_geometry:326-357 computes Delta_iota from the flux-spline
ratio at the containing equilibrium-data interval's endpoints and the assembly
point, divided by Nfp. This is not the FE radial-cell width: refining FE cells
at fixed equilibrium data does not remove this input-data regularization scale.
At zero width, GLISS uses the exact inverse if cancellation is resolved; an
unresolved nonconstant harmonic is accepted only with compatible forcing,
then assigned zero gauge (export_surface_geometry:741-753).
At positive width incompatible resonant forcing is attenuated rather than rejected.

The actual primitive path calls this MDE (primitive_kernel_geometry:81/88;
export_surface_geometry:319-324). Its beta derivatives enter jhat dot B
(field10) and the drive (344-347,364-381), and so stiffness.
For has_chart_metric=true, field13 beta is reconstructed geometrically at
355-356: do not claim every beta field comes directly from the regularized solve.
The generated local energy uses field10 and drive (gen/gliss_energy_gen:91-112).

Raw P,T,n arrays are not directly comparable across charts. Existing internal
flip-pol normalization maps T_T=-T_G, P_T=-Nfp*P_G
(terpsichore_normalization:48-53); equilibrium harmonics n_T=Nfp*n_G
(tprgl0:571), hence D_T=-Nfp*D_G. Compare common-chart attenuation, with
W_T^2=Nfp^2 W_G^2. Signs and 2*pi factors cannot remove its different denominator.

## Exact rational consequences and domain

Put A=D^2+delta(D+P)^2, B=D^2+W^2. For A*B!=0,

  D R_T-1=-delta(D+P)^2/A,
  D R_G-1=-W^2/B,
  R_T-R_G=D[W^2-delta(D+P)^2]/(A B).

Thus at D!=0 the responses agree exactly iff W^2=delta(D+P)^2.
At fixed delta>0 the TERPSICHORE inverse generally does not become 1/D
under radial refinement; the GLISS limit W->0 is 1/D only off resonance.
The formal delta->0 limit is likewise off-resonance, but TERPSICHORE enforces
a 1e-16 floor. At m=0 GLISS W=0, while TERPSICHORE attenuation still differs
generically. For delta>0, A=0 iff D=P=0, so the source formula is undefined
there; active QAS3 nonzero poloidal flux avoids that corner. B=0 iff D=W=0;
GLISS handles that corner by the separate compatibility branch above.
Neither regularized inverse is an exact interval principal-value average.
No convergence through a forced rational surface is established here.
Strict FortSym receipts are being prepared independently for these identities;
this source report itself is an algebraic derivation, not an executable proof receipt.

## Smallest useful matched control

First qualify one common equilibrium on the actual trial support (force balance,
metric/Jacobian, Fourier and radial reconstruction, mode conventions, finite
axis energy and fixed-wall traces), and bind source, deck and export bytes.
Freeze a small retained harmonic set, common radial points and forcing, then
compare source responses and force-balance residuals after the existing chart map.
Hold geometry/profiles, Fourier table, trial space, axis constraints and quadrature
fixed; only substitute these two explicit MDE responses into a bounded trial-energy
control. Measure energy change before any inertia/eigenvalue attribution.
An off-resonant test may require both delta(D+P)^2/D^2 and W^2/D^2 below a
frozen error budget; near rational surfaces require a declared regularized domain
or compatible forcing, not impossible equality to the singular ideal inverse.
A matched-response diagnostic does not certify either complete independent code.
