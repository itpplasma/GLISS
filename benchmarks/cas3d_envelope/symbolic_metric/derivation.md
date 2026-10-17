# Exact tangent metric averaging identity

For arbitrary real two-component vectors $a_0=(x_0,z_0)$ and
$a_1=(x_1,z_1)$, expansion of squares gives

$$
\frac{|a_0|^2+|a_1|^2}{2}
-\left|\frac{a_0+a_1}{2}\right|^2
=\frac{(x_1-x_0)^2+(z_1-z_0)^2}{4}
=\frac{|a_1-a_0|^2}{4}.
$$

The right side is nonnegative for real components, and is zero exactly when
$a_0=a_1$. FortSym's native engine with explicit real assumptions and the
SymEngine exact zero-test adapter each prove the displayed equality as a
zero residual: **2/2 strict checks proved, 0 probes, 0 failures**. The sign
and equality characterization follow from squares of real numbers; those
inequalities were not submitted as separate FortSym identity checks.

## Available VMEC++ source contract

Read-only source checkout `/Users/ert/code/vmecpp`, HEAD
`a40d22b0181c7378d44a5a1115a32dc136b35bc6`; exact file hashes are in the
separate receipt. The available kernels are under
`src/vmecpp/cpp/vmecpp/vmec/ideal_mhd_model/`.

- `metric_kernel.h:45–50` computes `guu` by separate quadratic averages,
  including parity-dependent terms. In the even-only slice
  `ruo_i=ruo_o=zuo_i=zuo_o=0`, set
  $a_0=(\mathrm{rue_i},\mathrm{zue_i})$ and
  $a_1=(\mathrm{rue_o},\mathrm{zue_o})$; the result is the first term
  in the displayed identity.
- `jacobian_kernel.h:51–52` computes `ru12` and `zu12` by averaging
  tangent components, including odd-parity contributions. In the same
  even-only slice, their squared norm is the second term in the identity.
- `jacobian_kernel.h:56` subsequently contracts those interpolated
  tangent components with radial derivatives in `tau1`.

Thus this exact algebra explains a possible discrete distinction between
averaging quadratic metric quantities and forming the metric from averaged
tangents. It does not assert that the complete parity-dependent VMEC++
metric equals this even-only specialization, that either discretization is
incorrect, or that this distinction causes a force residual. No scientific
sampling, production changes, or exhaustive VMEC proof are involved.

## Reproduction

```sh
cmake -S src -B build -G Ninja \
  -DCMAKE_Fortran_COMPILER=/opt/homebrew/bin/gfortran \
  -DCMAKE_PREFIX_PATH=/tmp/fortsym-ordered-outputs-prefix \
  -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 2
ctest --test-dir build --output-on-failure
```

The native build uses GNU Fortran 16.1.0, FortSym
`2f0571e39f36a0606b3ae150c8b3b3df92ce6750`, `-O3 -Wall -Wextra
-Werror=array-temporaries -fcheck=all`; the registered test passed in 0.27 s.
This proof and receipt are independent of the frozen 17-identity axis proof.
