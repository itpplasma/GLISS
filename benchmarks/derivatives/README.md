# Fixed-geometry pressure reverse assembly

`pressure_vjp.py` compares the pressure trace VJP with a gradient constructed
from one public pressure trace JVP per sample. Both actions use the same
admitted fixed-boundary equilibrium, discretization and mass-normalized
cluster basis. The public C ABI is called through the Python-owned problem;
the eigensolve and spectral-gap admission happen before the timed actions.

```sh
GLISS_LIB="$PWD/build/libgliss_c.so" PYTHONPATH=python \
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  python benchmarks/derivatives/pressure_vjp.py \
  --output /tmp/pressure-vjp-threads1.json
```

Use `libgliss_c.dylib` on macOS. The default committed input is
`test/data/solovev_q1.045.nc`, with mode `(1,1)`, degree 1, two radial cells,
64×16 angular points and the highest isolated eigenvalue. Pass a fixture
as the positional argument, repeat `--mode M N` to change the mode table,
and use `--cluster-size 2` for an admitted two-dimensional cluster.
`--degree`, `--radial-cells`, `--angular`, `--parity` and `--gap` select the
same public solver controls. Free-boundary cases are outside this benchmark.

Each method has an untimed warm-up; three measured repetitions alternate
execution order. `--warmup` and `--repeat` change these counts. Every
repetition checks the reverse gradient against all cardinal JVPs and the
public Python sensitivity with relative tolerance `2e-9` and absolute
tolerance `1e-11`. A failed admission or comparison stops the run. There is
no timing pass threshold.

The JSON records individual times and their distribution, separate setup
stages, agreement errors, gradient, input/library/runner hashes, source
revision and dirty status, linked-library paths/hashes, host, Python/NumPy,
and the OpenMP/BLAS thread environment. Set thread counts explicitly in each
fresh process when comparing scaling. CMake cache and affected compilation
commands are recorded when available; `--build-dir` locates them for an
installed library. These files describe the contemporary build and are not
embedded binary attestations. The library hash must stay fixed during a run.

The timing ratio compares two algorithms for this exact imported operator.
It does not certify geometry or mesh convergence, differentiate force
balance, or establish superiority over other solvers or hardware.
