#!/usr/bin/env bash
# Critical conformal-wall distance of the n=1 free-boundary kink of the GPEC
# Solov'ev family at q0=1.5, from public sources: GPEC DCON (VACUUM ishape=6)
# and GLISS on a GVEC export of the same analytic equilibrium.
# usage: benchmarks/solovev/free_boundary/run.sh OUTPUT_DIRECTORY
# Environment: GVEC_VENV (OUTPUT/venv), GLISS_LIB (bundled library),
# THREADS (4), MMAX (8), EDGE ("128 128"). Needs what
# benchmarks/solovev/dcon/run.sh and benchmarks/solovev/public/run.sh need.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../../.." && pwd)
out=$(mkdir -p "${1:?usage: run.sh OUTPUT_DIRECTORY}" && cd "$1" && pwd)
q0=1.5
threads=${THREADS:-4}; mmax=${MMAX:-8}
walls_dcon=(0.10 0.12 0.14 0.15 0.16 0.17 0.18 0.20)
edge=${EDGE:-128 128}

# DCON: build GPEC once through the fixed-boundary script, then rerun the
# regression example with the vacuum on and a conformal wall at distance
# a times the plasma half-width (a = 20 means no wall).
"$repo/benchmarks/solovev/dcon/run.sh" "$out/dcon" "$q0" > /dev/null
gpec=$out/dcon/GPEC
example=$gpec/docs/examples/solovev_ideal_example
echo "wall_half_widths,total_energy" > "$out/dcon_walls.csv"
for a in "${walls_dcon[@]}" 20; do
    run=$out/dcon_free/a$a
    mkdir -p "$run"
    cp "$example"/{equil,dcon,sol,vac}.in "$run"/
    sed -i 's/^\(\s*\)vac_flag=[tf]/\1vac_flag=t/; s/^\(\s*\)bin_euler=t/\1bin_euler=f/' "$run/dcon.in"
    sed -i "s/^\(\s*\)q0 = [0-9.]*/\1q0 = $q0/" "$run/sol.in"
    sed -i "s/^\(\s*\)a = [0-9.]*/\1a = $a/" "$run/vac.in"
    grep -q "vac_flag=t" "$run/dcon.in"
    grep -q "q0 = $q0" "$run/sol.in"
    grep -q "a = $a" "$run/vac.in"
    (cd "$run" && OMP_NUM_THREADS=1 "$gpec/bin/dcon" > stdout.log 2>&1)
    energy=$(sed -n 's/.*Energies:.*real = *\([^,]*\),.*/\1/p' "$run/stdout.log")
    echo "$a,$energy" >> "$out/dcon_walls.csv"
done

# GLISS: the GVEC export of the public Solov'ev pipeline (48 surfaces, M=16).
venv=${GVEC_VENV:-$out/venv}
if [ ! -x "$venv/bin/pygvec" ]; then
    python3 -m venv "$venv"
    "$venv/bin/pip" install -q numpy scipy netCDF4 xarray "gvec==1.5.0"
    package=$("$venv/bin/python" -c 'import gvec, os; print(os.path.dirname(gvec.__file__))')
    patch -d "$package" -p3 < "$repo/benchmarks/solovev/public/gvec-1.5.0-cas3d.patch"
fi
equilibrium=$out/gvec/q$q0
"$venv/bin/python" "$repo/benchmarks/solovev/public/make_gvec_input.py" "$q0" "$equilibrium"
(cd "$equilibrium" && OMP_NUM_THREADS=$threads "$venv/bin/pygvec" run parameter.toml > gvec.log 2>&1)
(cd "$equilibrium" && OMP_NUM_THREADS=$threads "$venv/bin/pygvec" to-cas3d --ns 48 \
    --MN_out 16 0 --winding -1 -o export.nc > export.log 2>&1)
# Bisect the critical wall between a stable 0.12 and an unstable 0.20.
OMP_NUM_THREADS=$threads python3 "$here/scan.py" "$equilibrium/export.nc" \
    --bisect 0.12 0.20 --poloidal-max "$mmax" --edge $edge \
    > "$out/gliss_walls.csv"
{
    git -C "$repo" rev-parse HEAD
    cat "$out/dcon/provenance.txt"
    "$venv/bin/pip" show gvec | grep -i '^version'
    sha256sum "$here"/*.py "$here"/run.sh
} > "$out/provenance.txt"
cat "$out/dcon_walls.csv" "$out/gliss_walls.csv"
