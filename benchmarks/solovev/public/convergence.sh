#!/usr/bin/env bash
# Radial convergence of the Solov'ev n=1 eigenvalue in FEEC degree and ns.
# usage: benchmarks/solovev/public/convergence.sh OUTPUT_DIRECTORY [Q0]
# Reuses OUTPUT/q$Q0 from run.sh (it runs run.sh when missing), exports the
# same GVEC state at NS_LIST surfaces (default "16 32 64 128") with M = 8, and
# prints the lowest eigenvalue of the axisymmetric family (poloidal_max 6)
# for degrees 1 to 4. Environment: GVEC_VENV, GLISS_LIB, THREADS (4).
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../../.." && pwd)
out=$(mkdir -p "${1:?usage: convergence.sh OUTPUT_DIRECTORY [Q0]}" && cd "$1" && pwd)
q0=${2:-1.035}
venv=${GVEC_VENV:-$out/venv}
threads=${THREADS:-4}
run=$out/q$q0
[ -f "$run/export.nc" ] || "$here/run.sh" "$out" "$q0"
for ns in ${NS_LIST:-16 32 64 128}; do
    [ -f "$run/export_ns$ns.nc" ] || (cd "$run" && OMP_NUM_THREADS=$threads \
        "$venv/bin/pygvec" to-cas3d --ns "$ns" --MN_out 8 0 --winding -1 \
        --stellsym -o "export_ns$ns.nc" > "export_ns$ns.log" 2>&1)
done
export GLISS_LIB=${GLISS_LIB:-$repo/build/libgliss_c.so}
OMP_NUM_THREADS=$threads PYTHONPATH=$repo/python python3 - "$run" ${NS_LIST:-16 32 64 128} <<'PY'
import sys
import gliss
run, sizes = sys.argv[1], [int(value) for value in sys.argv[2:]]
print("ns,degree,negative_count,lowest_eigenvalue")
for ns in sizes:
    with gliss.Equilibrium(f"{run}/export_ns{ns}.nc") as equilibrium:
        for degree in (1, 2, 3, 4):
            result = gliss.solve_axisymmetric(
                equilibrium, poloidal_max=6, degree=degree
            )
            print(f"{ns},{degree},{result.negative_count},"
                  f"{result.lowest_eigenvalue:.10e}", flush=True)
PY
