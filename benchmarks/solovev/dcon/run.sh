#!/usr/bin/env bash
# Rebuild GPEC DCON from its public repository and count fixed-boundary n=1
# Newcomb zero crossings for the GPEC analytic Solov'ev family.
# usage: benchmarks/solovev/dcon/run.sh OUTPUT_DIRECTORY [Q0 ...]
# Needs gfortran, OpenBLAS/LAPACK and the NetCDF C and Fortran libraries.
set -euo pipefail
commit=f5595c0689c624834d720068ddc8b5a7e4028248
out=$(mkdir -p "${1:?usage: run.sh OUTPUT_DIRECTORY [Q0 ...]}" && cd "$1" && pwd)
shift
q0s=("$@")
[ ${#q0s[@]} -gt 0 ] || q0s=(1.030 1.035 1.039062 1.039843 1.045)
multiarch=/usr/lib/$(gcc -print-multiarch 2>/dev/null || echo x86_64-linux-gnu)
export FC=gfortran CC=gcc
export NETCDFHOME=$multiarch NETCDF_DIR=$multiarch NETCDF_FORTRAN_HOME=$multiarch
export NETCDF_C_HOME=$multiarch NETCDFINC=/usr/include
export LAPACKHOME=$multiarch OPENBLASHOME=$multiarch
build() {  # build DIRECTORY
    local tree=$1
    [ -x "$tree/bin/dcon" ] && return
    if [ ! -d "$tree/.git" ]; then
        git init -q "$tree"
        git -C "$tree" remote add origin https://github.com/PrincetonUniversity/GPEC
    fi
    git -C "$tree" fetch -q --depth 1 origin "$commit"
    git -C "$tree" checkout -q FETCH_HEAD
    make -C "$tree/install" dcon > "$tree.build.log" 2>&1
}
build "$out/GPEC"
gpec=$out/GPEC
example=$gpec/docs/examples/solovev_ideal_example
csv=$out/dcon_solovev.csv
# DCON classifies a sign change of its critical eigenvalue as a Newcomb zero
# (not a pole) from one linearly interpolated midpoint.  With the example's
# tolerances (tol_nr=1e-6, tol_r=1e-7) the output steps are coarse enough that
# genuine crossings are rejected at isolated q0 (1.030, 1.035); tolerances of
# 1e-8 give a monotone count.  Both are reported.
echo "q0,newcomb_zero_count,tight_tolerance_newcomb_zero_count" > "$csv"
count() {  # count BINARY RUN_DIRECTORY
    (cd "$2" && OMP_NUM_THREADS=1 "$1" > stdout.log 2>&1)
    awk -F, '$1 == "newcomb_zero_count" {printf "%d", $3}' "$2/fixed_boundary.csv"
}
for q0 in "${q0s[@]}"; do
    row=$q0
    for variant in default tight; do
        run=$out/$variant/q$q0
        mkdir -p "$run"
        cp "$example"/{equil,dcon,sol,vac}.in "$run"/
        # Fixed boundary (no vacuum), axis scan from q=0.5, fixed-boundary output.
        sed -i 's/^\(\s*\)vac_flag=[tf]/\1vac_flag=f/; s/^\(\s*\)qlow=[0-9.]*/\1qlow=0.5/' "$run/dcon.in"
        sed -i 's/^\(\s*\)netcdf_out=t/\1netcdf_out=t\n\1out_fixed=t/' "$run/dcon.in"
        sed -i "s/^\(\s*\)q0 = [0-9.]*/\1q0 = $q0/" "$run/sol.in"
        grep -q "out_fixed=t" "$run/dcon.in"
        grep -q "q0 = $q0" "$run/sol.in"
        if [ $variant = tight ]; then
            sed -i 's/^\(\s*\)tol_nr=[0-9.eE+-]*/\1tol_nr=1e-8/; s/^\(\s*\)tol_r=[0-9.eE+-]*/\1tol_r=1e-8/' "$run/dcon.in"
            grep -q "tol_nr=1e-8" "$run/dcon.in"
        fi
        row=$row,$(count "$gpec/bin/dcon" "$run")
        rm -f "$run"/*.bin
    done
    echo "$row" >> "$csv"
done
{
    echo "GPEC $commit"
    gfortran --version | head -1
    sha256sum "$0"
} > "$out/provenance.txt"
cat "$csv"
