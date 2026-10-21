# Sparse stability

Configure with `-DGLISS_ENABLE_FORTSPARSE=ON` to use pinned
[FortSparse](https://github.com/lazy-fortran/fortsparse) native LDLT for the
shifted self-adjoint pencils. This option retains global pivoting across radial
interfaces without allocating a dense global matrix. The default build keeps
the existing block/dense path for reference. Native external backends are
not required. Offline builds may provide the pinned FortNum, FortSparse and
FortGen runtime trees through CMake's `FETCHCONTENT_SOURCE_DIR_*` options.

For SPD mass M, Sylvester's law gives the number of eigenvalues below sigma
from the negative pivots of K-sigma M. An uncertain or singular pivot is not
accepted as a stability sign; existing bracketing handles unresolved probes.
FortSparse's reliability threshold is a floating-point heuristic, not an
outward-rounded proof of the exact continuous operator's sign. Mesh, mode and
equilibrium convergence remain separate requirements.

The adapter uses `fortnum_precision` stage policies. This first increment
supports `cpu-fp64` only; candidate evaluation, factors, accepted residuals,
estimators and certificates all remain FP64. No mixed-precision result is
claimed. Lower-precision policies are explicitly rejected by the adapter's
admission query until corresponding kernels and error controls are supported.

C0 radial assembly uses two colors of nonadjacent cells after the first cell
fixes axis orientation. Each color writes disjoint scalar entries. Block matrix
applications also assign disjoint output rows to OpenMP threads. Set
`OPENBLAS_NUM_THREADS=1` when measuring `OMP_NUM_THREADS` strong scaling to
avoid nested BLAS oversubscription. The LDLT factorization itself is currently
serial; assembly speedup is not a claim of parallel factorization.

`-DGLISS_BUILD_BENCHMARKS=ON` builds `gliss_sparse_performance`. It reports
actual team size, assembly/solve wall times, eigenvalue, certificate and inertia
for a fixed Solov'ev export, cell count and degree. Its final argument selects
sparse (1) or dense (0) matrix storage. The sparse build changes factorization;
the default build supplies the independent dense/block comparison. Compare
physical eigenvalues within the reported certificates and solver stopping
accuracy, rather than requiring last-digit equality.

The opt-in native backend uses bounded BLAS-3 LDLT panels with OpenMP assembly
and scatter. Panels needing an interchange recover through unrestricted sparse
Bunch–Kaufman pivoting, preserving the inertia and near-zero rejection contract.
The block backend remains the default: the native panels improve the old sparse
path but have not beaten the block solver on the measured large controls.
