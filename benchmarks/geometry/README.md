# Cartesian geometry value-kernel timing

The generated primal include matches the original caller's speed on this
single-host regular-point workload. The recorded medians are 0.192735 s for
the generated caller and 0.198593 s for the original: a time ratio of 0.97050.
Treat the small difference as comparable performance, not a general speedup.

The baseline is `src/cartesian_primitive_geometry.f90` from
`f2da7f4301e9207db38d3fd5fc3c1488b344fb56`. The recipe obtains it with
`git show` and changes only its module name so both implementations can link.
No historical source copy is maintained here.

Run from a quiet host with GNU Fortran and Python 3:

```sh
python3 benchmarks/geometry/run.py --compiler gfortran --output /tmp/geometry-results.json
```

The script compiles frozen copies of both callers, the current generated
module/include and [benchmark.f90](benchmark.f90) using the same `-O3
-ffree-line-length-none -Igenerated` flags. It uses separate compilation,
without LTO, fast-math or OpenMP compilation, and runs with `OMP_NUM_THREADS=1`.
The temporary build directory, exact commands, source hashes and raw timings
are retained and recorded in JSON. It does not run the project's test suite.
The baseline commit must be present in the local Git history.

The harness visits 4096 changing, regular, nonorthogonal Cartesian jets,
including arbitrary symmetric second derivatives, changing flux slopes and
field periods, and 2048 samples of each Jacobian sign. Before timing, all
22 output components must be finite and agree with the original within
`1e-13 * max(1, abs(original), abs(current))`. The maximum recorded normalized
difference was `1.4432899320127035e-15`. This is an implementation comparison;
the independent physical and derivative oracles belong to the native tests.

After warmup, eleven pairs alternate execution order. Each block makes
8,388,608 calls and consumes a checksum of Jacobian, field magnitude, metric
and second form. All recorded checksums agree. The recipe rejects blocks
shorter than 100 ms; a faster host may require increasing the fixed harness
pass count and recording the changed source hash.

[results.json](results.json) records the actual product sources measured on
an Apple M3 Ultra with Homebrew GNU Fortran 16.1.0, including all eleven raw
pairs and SHA256 provenance. Its harness bytes are identical to the committed
benchmark; the recorded scratch filename was `geometry_timing.f90`.

The shared include preserves the generated arithmetic while allowing the
caller's private value routines to inline. During diagnosis, separate module
calls took 1.87 times the original time, and a global LTO control took 1.42
times. GNU's inline report identified three unavailable procedure bodies;
the private generated include removed those calls and restored the original
400-byte stack frame. These controls explain the include choice and are not
the recorded product result.
