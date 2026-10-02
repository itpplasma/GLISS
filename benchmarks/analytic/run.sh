#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
out=${1:?Usage: benchmarks/analytic/run.sh ABSOLUTE_OUTPUT_DIRECTORY}
mkdir -p "$out"
export FO_JOBS=2 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
fo exec benchmark_theta_pinch "$out/theta_spectrum.csv"
fo exec benchmark_theta_pinch "$out/legacy_spectrum.csv" legacy
fo exec benchmark_local_limits "$out/homogeneous.csv" "$out/helical.csv"
python3 benchmarks/analytic/plot.py "$out"
python3 benchmarks/analytic/plot.py "$out" --legacy
{
    git rev-parse HEAD
    git diff --binary HEAD | sha256sum
    rg '^(CMAKE_(Fortran|BUILD_TYPE)|BLAS_|LAPACK_|BLA_)' build/CMakeCache.txt
    rg '^set\(CMAKE_Fortran_(COMPILER|COMPILER_ID|COMPILER_VERSION)' \
        build/CMakeFiles/*/CMakeFortranCompiler.cmake
    ldd build/benchmark_theta_pinch
    sha256sum build/benchmark_theta_pinch build/benchmark_local_limits
    python3 -c 'import numpy, scipy, matplotlib; print(numpy.__version__, scipy.__version__, matplotlib.__version__)'
    sha256sum benchmarks/analytic/*
} > "$out/provenance.txt"
cp benchmarks/analytic/*.py benchmarks/analytic/*.f90 benchmarks/analytic/run.sh "$out/"
cp benchmarks/analytic/README.md "$out/README.md"
git diff --binary HEAD > "$out/source.patch"
