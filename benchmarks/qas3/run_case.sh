#!/bin/bash
# Run one QAS3 case through VMEC, TERPSICHORE and GLISS.
#
# usage: run_case.sh NAME PRESSURE_FACTOR CURRENT_FACTOR
#
# Environment: WORK (from build.sh), GLISS_LIB (libgliss_c.so), MMS (8),
# NSMAX (5), AL0 (-1e-5, TERPSICHORE shift), INDEPENDENT=1 to also run the
# independent GLISS path (convert_vmec + CAS3D marginality; 15-35 min).
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
WORK=${WORK:-$HERE/work}
NAME=$1; PRESSURE=$2; CURRENT=$3
MMS=${MMS:-8}; NSMAX=${NSMAX:-5}; AL0=${AL0:--1e-5}
BENCH=$WORK/terpsichore/QAS3_fx_benchmark
CASE=$WORK/cases/$NAME
mkdir -p "$CASE"
cd "$CASE"

python3 "$HERE/make_input.py" "$BENCH/input.QAS3_fx_cur75_bench" \
    "input.$NAME" "$PRESSURE" "$CURRENT"
[ -f "wout_$NAME.nc" ] || "$WORK/STELLOPT/VMEC2000/Release/xvmec2000" \
    "input.$NAME" > vmec.log 2>&1

# TERPSICHORE: public VMEC text converter, fixed boundary, MODELK=0.
cp "wout_$NAME.txt" fort.8
"$WORK/terpsichore/srcs/vmecv92terps.x" > vmecv92terps.log 2>&1
python3 "$HERE/make_deck.py" "$BENCH/QAS3_fx_cur_bench_n1.data" deck.data \
    "$MMS" "-$NSMAX" "$NSMAX" 0 "$AL0" > /dev/null
"$WORK/terpsichore/srcs/tpr_ap.x" < deck.data > terpsichore.log 2>&1
grep -E "EIGENVALUE FROM|NUMBER OF NEGATIVE EIGEN|ITERATIONS DONE|NON CONVERGED" \
    fort.16 | tee terpsichore.txt
# Shifted inverse iteration separates eigenvalues only when the shift is close
# to the lowest one; an iteration count equal to NITMAX with nonconverged
# components means the reported pair is a mixture, not an eigenpair.
if grep -q "NON CONVERGED Y COMPONENTS = *[1-9]" fort.16; then
    echo "WARNING: TERPSICHORE did not converge; move AL0 closer to the lowest eigenvalue"
fi

python3 "$HERE/compare.py" replay fort.23 | tee replay.txt
python3 "$HERE/compare.py" mercier "wout_$NAME.nc" "$NAME" | tee mercier.txt
if [ "${INDEPENDENT:-0}" = 1 ]; then
    python3 "$HERE/compare.py" inertia "${NAME}_gliss.nc" deck.data.modes \
        | tee inertia.txt
fi
rm -f fort.8 fort.17 fort.22 fort.23 fort.73
