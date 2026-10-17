# Uniform-grid axis and BOOZ interpolation identities

This is an exact interpolation result under the assumptions below. It does
not establish the cause of a production force residual or prescribe a physics
repair. No repository source, upstream source, equilibrium or exporter was
changed.

## Assumptions and notation

Let the full radial nodes be `s_0=0`, `s_1=h`, `s_2=2h`, with **real `h>0`**.
Write `q_i=q(s_i)` and `M_i` for the interpolating cubic spline's second
derivative at node `i`. All coefficients are real. The later half-grid result
holds for any real `j`; subsequent physical half-grid cells have integer
`j>=1` and center `(j+1/2)h`.

The conditional spline result assumes exactly these two equations:

\[
K=M_0-2M_1+M_2=0,\qquad
R=M_0+4M_1+M_2-\frac{6(q_2-2q_1+q_0)}{h^2}=0.
\]

The first is the not-a-knot condition on the first two equal intervals. The
second is the cubic spline continuity row after division by `h`. The positive
step assumption justifies both division by `h` and division by `h^2`.

## Eliminating the two rows

Put `D=q_2-2q_1+q_0`. Direct subtraction gives

\[
R-K=6M_1-6D/h^2,\qquad
M_1=(R-K)/6+D/h^2.
\]

Thus `K=R=0` implies `M_1=D/h^2`. If the axis value is manufactured by
`qhat_0=2q_1-q_2`, then `D=0` exactly and therefore **`M_1=0`**.

This is the curvature at the first positive node `s=h`. It does **not** force
the curvature at the axis `M_0` to vanish: the remaining local relation is
only `M_0=-M_2`. No uniqueness assumption about the full spline system is
needed for the implication for `M_1`.

## Exact quadratic model

Take `q(s)=a+b s+c s^2`. Its true samples and curvature are

\[
q_0=a,\quad q_1=a+bh+ch^2,\quad
q_2=a+2bh+4ch^2,\quad q''(s)=2c.
\]

The extrapolated axis value is `qhat_0=a-2ch^2`, so

\[
\widehat q_0-q(0)=-2ch^2.
\]

The true samples have `D=2ch^2`, and the spline row gives `M_1=2c`.
The extrapolated axis sample instead gives `D=0`, hence `M_1=0` and

\[
M_1-q''(h)=-2c.
\]

These are exact equalities for this model, not truncation estimates. When
`c=0`, the model is linear and all the interpolation biases below vanish.

## BOOZ first and later half-grid values

The inspected BOOZ_XFORM 0.0.9 source constructs odd-mode geometry using
`q=C/sqrt(s)`. For `m=1`, its first half-grid quotient is

\[
\widehat q_{1/2}=\tfrac32q_1-\tfrac12q_2
=\tfrac12(\widehat q_0+q_1).
\]

For the quadratic model,

\[
\widehat q_{1/2}-q(h/2)=-\tfrac34ch^2.
\]

On subsequent cells the quotient is the ordinary endpoint average:

\[
\widehat q_{j+1/2}=\tfrac12(q(jh)+q((j+1)h)),\qquad
\widehat q_{j+1/2}-q((j+1/2)h)=+\tfrac14ch^2.
\]

Consequently the first-cell quotient bias differs from the ordinary-cell
bias by `-ch^2`. These formulas concern `q`; the corresponding `m=1` geometry
coefficient bias also carries the factor `sqrt(s_half)`. They do not by
themselves quantify derivative, magnetic-field, current or force errors.

## Source correspondence and limits

The read-only source is
`/Users/ert/.cache/uv/sdists-v9/pypi/booz-xform/0.0.9/US9JLQ3OXA1SsSIK6_wyM/src/src/_booz_xform/init_from_vmec.cpp`.

- Lines 197–208 describe the parity quotient and axis extrapolation.
- Lines 210–219 construct the uniform full and half grids.
- Lines 246–258 average adjacent odd-mode quotients.
- Lines 260–280 implement the special first half-grid `m=1` value.

Its SHA256 is
`1023316194f30127f37e7975c953b919d046aa87c2b5477cc621e62f839f50ac`.
The existing read-only installed-reader audit is recorded separately in
`/tmp/gliss-booz-interpolation-audit-20261017/receipt.json`.

The spline implication is conditional on the stated rows and input closure;
this proof does not claim that BOOZ itself solves that spline system. The
explicit BOOZ code correspondence establishes its quotient interpolation
formula. Accuracy of the producer equilibrium and physical causality require
their own evidence. A nonzero interpolation bias alone does not establish a
production bug or that removing the closure improves force balance.

## Strict verification and reproduction

`src/axis_interpolation_proof.f90` records `h>0` in a native FortSym assumption
context and checks its implied reality and nonzero facts. It checks 17 exact
identities through strict `check_identity` using the native engine and the
SymEngine adapter: **34/34 proved, 0 probes, 0 failures**. UNKNOWN or a failed
identity terminates the driver. The exact rational fractions use integer
numerators and denominators; no floating-point samples enter the proof.

The native registered CTest passed with GNU Fortran 16.1.0, Release `-O3`,
`-Wall -Wextra -Werror=array-temporaries -fcheck=all`, and at most two build
jobs. The dependency is FortSym revision
`2f0571e39f36a0606b3ae150c8b3b3df92ce6750` installed at
`/tmp/fortsym-ordered-outputs-prefix`. File hashes and commands are in
`receipt.json`; the strict suite counters are in `build/strict-identities.json`.

```sh
cmake -S src -B build -G Ninja \
  -DCMAKE_Fortran_COMPILER=/opt/homebrew/bin/gfortran \
  -DCMAKE_PREFIX_PATH=/tmp/fortsym-ordered-outputs-prefix \
  -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 2
ctest --test-dir build --output-on-failure
```
