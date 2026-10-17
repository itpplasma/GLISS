# Weighted derivative-space proof

This is a bounded proof for the candidate tangential space described in the
[readiness audit](../derivative_space_audit.md). The [full argument](audit.md)
proves its exact kernel and dimension, distinguishes accepted custom powers
from physical defaults, and derives the axis normalization and wall-generator
requirements. It does not establish a Solov'ev repair.

The native CMake/CTest suite checks 13 identities strictly, with zero numerical
probes: one unfactored derivative identity with native FortSym, and six
normalized or algebraic identities with both native FortSym and SymEngine.
The function-space rank and trace conclusions follow from the written
mathematical argument, rather than numerical ranks or the identity count.

With the repository-pinned FortSym installed, reproduce with:

```sh
cmake -S benchmarks/solovev/public/symbolic_derivative/src \
  -B /tmp/gliss-weighted-derivative-proof -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/path/to/fortsym-prefix
cmake --build /tmp/gliss-weighted-derivative-proof --parallel 2
ctest --test-dir /tmp/gliss-weighted-derivative-proof --output-on-failure
```

The original [receipt](receipt.json) binds sources, library and results from
the isolated execution. Its build and diagnostic log paths refer to that
execution directory; they are not required to reproduce the committed suite.
The audit records unsuccessful unfactored normalizations separately from the
proved identities, including an inconsistent SymEngine NONZERO verdict. No
failed normalization is treated as a disproof or as evidence for a repair.
