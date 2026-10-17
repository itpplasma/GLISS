# Weighted derivative-space proof

This is a bounded proof for the candidate tangential space described in the
[readiness audit](../derivative_space_audit.md). The [full argument](audit.md)
proves its exact kernel and dimension, distinguishes accepted custom powers
from physical defaults, and derives the axis normalization and wall-generator
requirements. It does not establish a Solov'ev repair.

The current native CMake/CTest suite checks 14 identities strictly, with zero
numerical probes: the original unfactored derivative identity and six
normalized or algebraic identities with both native FortSym and SymEngine.
The additional SymEngine check exercises the repaired exponent normalization;
its [validation receipt](repair-validation.json) records the tested source,
repair patch and library separately from the historical evidence.
The function-space rank and trace conclusions follow from the written
mathematical argument, rather than numerical ranks or the identity count.

With a FortSym build containing the recorded repair installed, reproduce with:

```sh
cmake -S benchmarks/solovev/public/symbolic_derivative/src \
  -B /tmp/gliss-weighted-derivative-proof -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/path/to/fortsym-prefix
cmake --build /tmp/gliss-weighted-derivative-proof --parallel 2
ctest --test-dir /tmp/gliss-weighted-derivative-proof --output-on-failure
```

The original [receipt](receipt.json), [13-check result](strict-identities.json)
and [audit](audit.md) remain the frozen evidence from FortSym `2f0571e`. Their
source hashes refer to that historical source, retained in Git history; the
new validation records the four-line source addition and the 14-check result.
Build and diagnostic log paths in the original receipt refer to its isolated
execution directory and are not required to reproduce the current suite.
The audit's unsuccessful unfactored normalizations describe that earlier
execution, including the inconsistent SymEngine NONZERO verdict. They are
not treated as mathematical disproofs or rewritten as successful checks.
