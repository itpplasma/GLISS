#!/usr/bin/env bash
# Run inside a Slurm allocation: bash ci/cluster-validate.sh CHECKOUT RUN_ROOT
set -euo pipefail
: "${SLURM_JOB_ID:?Run this script inside a Slurm job}"
gliss_source=$(realpath "${1:?Git checkout required}")
gliss_run=$(realpath -m "${2:?Run root required}")
# Benchmark provenance tests require the real checkout, including its history.
git -C "$gliss_source" rev-parse --verify HEAD >/dev/null
if [ -n "$(git -C "$gliss_source" status --porcelain)" ]; then
    echo 'Cluster validation requires a clean Git checkout' >&2
    exit 1
fi
mkdir -p "$gliss_run"
exec >"$gliss_run/validation.log" 2>&1
unset DISPLAY LD_LIBRARY_PATH
source /home/ert/opt/kin6d-env.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg
export PYTHONPATH="$gliss_source/python"
gliss_build="$gliss_run/build-$SLURM_JOB_ID"
git -C "$gliss_source" rev-parse HEAD >"$gliss_run/revision.txt"
"$FC" --version | head -1 >"$gliss_run/compiler.txt"
cmake -S "$gliss_source" -B "$gliss_build" -G Ninja \
    -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_FLAGS=-Werror \
    -DGLISS_REQUIRE_OPENBLAS=ON "${@:3}" >"$gliss_run/configure.log" 2>&1
cmake --build "$gliss_build" -j "${SLURM_CPUS_PER_TASK:-8}" \
    >"$gliss_run/build.log" 2>&1
export GLISS_LIB="$gliss_build/libgliss_c.so"
sha256sum "$GLISS_LIB" >"$gliss_run/binary.sha256"
ldd "$GLISS_LIB" >"$gliss_run/binary.ldd"
ctest --test-dir "$gliss_build" -j "${SLURM_CPUS_PER_TASK:-8}" \
    --output-on-failure >"$gliss_run/native.log" 2>&1
python3 -m venv "$gliss_run/venv"
"$gliss_run/venv/bin/pip" install 'numpy>=1.23,<2' pytest scipy \
    'booz-xform>=0.0.9,<0.1' >"$gliss_run/python-install.log" 2>&1
cd "$gliss_source"
"$gliss_run/venv/bin/python" -m pytest python/tests -q \
    >"$gliss_run/python.log" 2>&1
"$gliss_run/venv/bin/pip" freeze >"$gliss_run/python-packages.txt"
