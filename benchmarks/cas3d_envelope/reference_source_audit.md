# Reference source-contract audit receipt
Date: 2026-10-17. Read-only bounded audit; no builds, solves, numerical measurements, source edits or worker commits.
Scope: accessible VMEC++ discrete geometry/field/current/residual contracts and BOOZ 0.0.9 full-to-half geometry interpolation; not CAS3D implementation certification.
VMEC checkout root: /Users/ert/code/vmecpp
VMEC current source/wrapper commit: a40d22b0181c7378d44a5a1115a32dc136b35bc6
Installed VMEC metadata: 0.5.5.dev120+g1b36dfca2.d20260627.
Metadata-associated source commit: 1b36dfca2da37bb9b44952916efdf41af2f8da8d; exact native binary-to-source build association remains UNPROVEN.
Installed native path: /Users/ert/code/vmecpp/.venv/lib/python3.12/site-packages/vmecpp/cpp/_vmecpp.cpython-312-darwin.so
Installed native SHA256: ffdfd300ec932be420e3feed4cbe399ced4fd82b78eb9c3b818a61af341452df (previous frozen provenance receipt; no rebuild).
Provenance receipt: /tmp/gliss-w7x-source-radial289-extended-20261017/settings-frozen.json SHA256 4b3939c834a1e9b584ce04b1dbc419a8a51ca076b619b8283a285ed2ea004b81.
Audited source files below: SHA256(current checkout bytes); SHA256(metadata-associated commit bytes); equality is source evidence, not binary attestation.
- src/vmecpp/cpp/vmecpp/vmec/ideal_mhd_model/jacobian_kernel.h: 1d74c923171911532b95d66cbcb9e4431dd95ba2a9212a1cf391a21cd1542932; 1d74c923171911532b95d66cbcb9e4431dd95ba2a9212a1cf391a21cd1542932; identical=True.
- src/vmecpp/cpp/vmecpp/vmec/ideal_mhd_model/metric_kernel.h: 898f3753300802757aab56429a112b4bb7637ccfb0d4657c38628945683d83ab; 898f3753300802757aab56429a112b4bb7637ccfb0d4657c38628945683d83ab; identical=True.
- src/vmecpp/cpp/vmecpp/vmec/ideal_mhd_model/bcontra_kernel.h: f47b9b96d064793ce78e08c781f7148dc967082c804cfae1e6a74296c01a69a8; f47b9b96d064793ce78e08c781f7148dc967082c804cfae1e6a74296c01a69a8; identical=True.
- src/vmecpp/cpp/vmecpp/vmec/ideal_mhd_model/mhdforce_kernel.h: 8fd5cabdc368b2410975683f1c260726fdcae321164268c435d49650e8841b04; 8fd5cabdc368b2410975683f1c260726fdcae321164268c435d49650e8841b04; identical=True.
- src/vmecpp/cpp/vmecpp/vmec/ideal_mhd_model/ideal_mhd_model.cc: 2f4ae7d8e3f5bb9fec0ccc23e45c2b0030e9a838556ad14b1b419b84e749ea0f; 4b223c481186378573c53b6351dd966a34b5beea83764a0334437d5e1924d806; identical=False.
- src/vmecpp/cpp/vmecpp/vmec/fourier_forces/fourier_forces.cc: 30916befc59524190b3f4331f6813d52467dec4eac2921c8d5c9844601c6069b; 30916befc59524190b3f4331f6813d52467dec4eac2921c8d5c9844601c6069b; identical=True.
Revision comparison: git diff 1b36dfc a40d22b reports no changes in the five kernel/FourierForces files; ideal_mhd_model.cc has only six added lines in constraintForceMultiplier (1933 onward), adding an optional freeze branch. Audited geometry/current/residual arithmetic is unchanged.
BOOZ source root: /Users/ert/.cache/uv/sdists-v9/pypi/booz-xform/0.0.9/US9JLQ3OXA1SsSIK6_wyM/src/src/_booz_xform
- read_wout.cpp: SHA256 d444c59deee7dade7928afb52ad57bf0b8e77c485be8f592d55d07870b407843.
- init_from_vmec.cpp: SHA256 1023316194f30127f37e7975c953b919d046aa87c2b5477cc621e62f839f50ac.
BOOZ cached build extension SHA256: 4107084cc9a11eab78e2619815e8b2ca164527b585fbe18d94e02a996ef90e17; prior reference-assets audit matched the frozen native extension to this cached build.
VMEC geometry contract: jacobian_kernel.h:50-60 averages parity-separated R/Ru/Zu, differences adjacent full surfaces for Rs/Zs, then tau=tau1+(1/4)tau2; metric_kernel.h:43-75 sets gsqrt=R12*tau and separately averages quadratic metric products.
Exact algebraic distinction: for even-only affine tangent a=(Ru,Zu), g_uu,discrete-|a_mid|^2=|a_outer-a_inner|^2/4. Thus discrete half-grid geometry/metric need not equal the metric/determinant of a separately interpolated smooth geometry.
VMEC field contract: bcontra_kernel.h:19-41 averages normalized angular lambda derivatives and divides by discrete gsqrt; ideal_mhd_model.cc:1291-1324 adds phi-prime before averaging, and :1393-1398 adds chi-prime/gsqrt (lambda_u plus, lambda_v minus convention).
VMEC current/residual contract: ideal_mhd_model.cc:1518-1547 uses half-to-full finite differences of angularly averaged covariant B and pressure; equiF is flux-surface averaged radial force balance. FourierForces::residuals at fourier_forces.cc:105-158 sums squared finite-mode force coefficients; model:322-324,1636-1638 normalizes them. These are not pointwise continuous Cartesian force norms.
BOOZ geometry contract: init_from_vmec.cpp:211-219 assumes uniform h=1/(ns-1), s_half=(k+1/2)h; :228-257 averages even-m f and odd-m q=f/sqrt(s), restoring sqrt(s_half).
BOOZ m=1 first-half contract: init_from_vmec.cpp:260-280 writes column0 as sqrt(h/2)*(1.5*q(h)-0.5*q(2h)), equivalent to q_axis=2*q(h)-q(2h); correct for linear q, error=-3*h^2*q_second(h/2)/8+O(h^3). No normal-grid m=1 indexing defect found.
BOOZ higher-m limit: for Cartesian-regular f=s^(m/2)*F(s), first-half leading-value ratios are 2^(m/2-1) for even m>=2 and 2^((m-3)/2) for odd m>=3. Relative error can be order one although absolute amplitudes vanish; high-order radial derivative accuracy is not established.
BOOZ NONBLOCKING UNRELATED STATIC FINDING: init_from_vmec.cpp:34 accepts ns=2, but the m=1 branch :270-280 reads index2. Requires a separate upstream reproducer. Not observed on these larger grids and not an explanation of the W7X discrepancy.
Control limit: the isolated positive-node-only m=1 reconstruction control worsened continuous force by 8-20% at all four studied shells (solver-validation receipt); no fix attribution or production change is justified by this control.
Acceptance limit: discrete solver convergence and reference agreement do not certify continuous force/current consistency after export/reconstruction. Require matched geometry, flux/J field and derivatives plus independent continuous force checks; exact installed VMEC build-source identity remains open.
No claim: neither this conceptual audit nor BOOZ representation controls certify unseen CAS3D implementation, a global spectrum, or a self-consistent corrected equilibrium.
