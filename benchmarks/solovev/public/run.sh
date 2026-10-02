#!/usr/bin/env bash
# Regenerate the GPEC Solov'ev n=1 comparison from public sources only.
# usage: benchmarks/solovev/public/run.sh OUTPUT_DIRECTORY [Q0 ...]
# Environment: GVEC_VENV (default OUTPUT/venv), NS (64), EXPORT_M (24),
# MMAX (8), DEGREE (2), GLISS_BIN (build/fo/bin), THREADS (4).
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../../.." && pwd)
out=$(mkdir -p "${1:?usage: run.sh OUTPUT_DIRECTORY [Q0 ...]}" && cd "$1" && pwd)
shift
q0s=("$@")
[ ${#q0s[@]} -gt 0 ] || q0s=(1.035 1.039062 1.039843 1.045)
venv=${GVEC_VENV:-$out/venv}
ns=${NS:-64}; export_m=${EXPORT_M:-24}; mmax=${MMAX:-8}; degree=${DEGREE:-2}
bin=${GLISS_BIN:-$repo/build/fo/bin}; threads=${THREADS:-4}
# Preserve the export and trial truncations while admitting the cubic MDE
# forcing and the displacement products. Round upward to a power of two.
case $export_m:$mmax in
    *[!0-9:]*|:*|*:) echo "EXPORT_M and MMAX must be integers" >&2; exit 2 ;;
esac
required_theta=$((6 * 10#$export_m + 1))
trial_theta=$((2 * 10#$mmax + 10#$export_m + 1))
if [ "$trial_theta" -gt "$required_theta" ]; then required_theta=$trial_theta; fi
angular_theta=64
while [ "$angular_theta" -lt "$required_theta" ]; do
    angular_theta=$((2 * angular_theta))
done

if [ ! -x "$venv/bin/pygvec" ]; then
    python3 -m venv "$venv"
    "$venv/bin/pip" install -q numpy scipy netCDF4 xarray
    "$venv/bin/pip" install -q "gvec==1.5.0"
    package=$("$venv/bin/python" -c 'import gvec, os; print(os.path.dirname(gvec.__file__))')
    patch -d "$package" -p3 < "$here/gvec-1.5.0-cas3d.patch"
fi
csv=$out/gliss_solovev.csv
echo "q0,negative_count,lowest_eigenvalue,certificate,force_balance_residual,angular_theta,angular_zeta" > "$csv"
for q0 in "${q0s[@]}"; do
    run=$out/q$q0
    "$venv/bin/python" "$here/make_gvec_input.py" "$q0" "$run"
    (cd "$run" && OMP_NUM_THREADS=$threads "$venv/bin/pygvec" run parameter.toml > gvec.log 2>&1)
    # which_hmap=1 has y = -R sin(zeta), so the Boozer frame needs winding -1.
    (cd "$run" && OMP_NUM_THREADS=$threads "$venv/bin/pygvec" to-cas3d --ns "$ns" \
        --MN_out "$export_m" 0 --winding -1 -o export.nc > export.log 2>&1)
    OMP_NUM_THREADS=$threads OPENBLAS_NUM_THREADS=1 "$bin/gliss_axisymmetric" \
        "$run/export.nc" 1 "$mmax" --degree "$degree" \
        --angular "$angular_theta" 8 > "$run/gliss.log"
    tail -1 "$run/gliss.log" | awk -F, -v q="$q0" \
        '{print q "," $11 "," $8 "," $9 "," $12 "," $14 "," $15}' >> "$csv"
done
{
    git -C "$repo" rev-parse HEAD
    git -C "$repo" diff HEAD | sha256sum
    "$venv/bin/pip" show gvec | grep -i '^version'
    sha256sum "$here"/*.py "$here"/*.patch "$here"/run.sh
} > "$out/provenance.txt"
cat "$csv"
