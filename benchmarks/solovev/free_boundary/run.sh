#!/usr/bin/env bash
# Critical conformal-wall distance of the n=1 free-boundary kink of the GPEC
# Solov'ev family at q0=1.5, from public sources: GPEC DCON (VACUUM ishape=6)
# and GLISS on a GVEC export of the same analytic equilibrium.
# usage: benchmarks/solovev/free_boundary/run.sh OUTPUT_DIRECTORY
# Environment: GLISS_EXPORT/GLISS_PARAMETER (qualified fixture), GLISS_LIB,
# THREADS (4), MMAX (8), EDGE ("128 128"). Needs what
# benchmarks/solovev/dcon/run.sh needs, plus NumPy, SciPy and netCDF4.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../../.." && pwd)
out=$(mkdir -p "${1:?usage: run.sh OUTPUT_DIRECTORY}" && cd "$1" && pwd)
q0=1.5
threads=${THREADS:-4}; mmax=${MMAX:-8}
walls_dcon=(0.10 0.12 0.14 0.15 0.16 0.17 0.18 0.20)
edge=${EDGE:-128 128}

# GLISS: use the qualified native true-edge export of the same public family.
# Overrides must supply a matching generator parameter file. Admission is a
# hard gate; extrapolated half-grid flux slopes cannot qualify this benchmark.
equilibrium=${GLISS_EXPORT:-$repo/test/data/solovev_q1.5_edge.nc}
parameter=${GLISS_PARAMETER:-$repo/test/data/solovev_q1.5_edge.parameter.toml}
python3 "$here/qualify_export.py" "$parameter" "$equilibrium" \
    --require-explicit-boundary > "$out/qualification.json"

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

# Bisect the critical wall between a stable 0.12 and an unstable 0.20.
OMP_NUM_THREADS=$threads python3 "$here/scan.py" "$equilibrium" \
    --bisect 0.12 0.20 --poloidal-max "$mmax" --angular 256 8 --edge $edge \
    > "$out/gliss_walls.csv"
{
    git -C "$repo" rev-parse HEAD
    cat "$out/dcon/provenance.txt"
    sha256sum "$equilibrium" "$parameter"
    sha256sum "$here"/*.py "$here"/run.sh
} > "$out/provenance.txt"
cat "$out/dcon_walls.csv" "$out/gliss_walls.csv"
