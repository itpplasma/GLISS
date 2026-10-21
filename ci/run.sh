#!/bin/bash
# Execute the repository's original native CI checks in isolated Slurm scratch.
set -euo pipefail
: "${SLURM_JOB_ID:?ci/run.sh must run in Slurm}"
tier=$(cat "$(dirname "$(dirname "$(realpath "$0")")")/../tier.txt")
src=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
root=$(dirname "$src")
cd "$src"
. "$root/task.env"
. /home/ert/opt/kin6d-env.sh
export PATH=/home/ert/.local/bin:$PATH
export GIT_CONFIG_GLOBAL=$root/gitconfig
: > "$GIT_CONFIG_GLOBAL"
git config --global protocol.file.allow always
shopt -s nullglob
for mirror in "$root"/mirrors/*/*.git; do
    slug=${mirror#"$root/mirrors/"}; slug=${slug%.git}
    git config --global --add "url.file://$mirror.insteadOf" "https://github.com/$slug.git"
    git config --global --add "url.file://$mirror.insteadOf" "https://github.com/$slug"
    git config --global --add "url.file://$mirror.insteadOf" "git@github.com:$slug.git"
done
export GITHUB_WORKSPACE=$src
export GITHUB_ENV=$root/github-env GITHUB_PATH=$root/github-path GITHUB_OUTPUT=$root/github-output
: > "$GITHUB_ENV"; : > "$GITHUB_PATH"; : > "$GITHUB_OUTPUT"
export CMAKE_BUILD_PARALLEL_LEVEL=${SCI_CPUS:-4} CTEST_PARALLEL_LEVEL=${SCI_CPUS:-4}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
mkdir -p ci-out
{ hostname; /usr/bin/gfortran --version | head -1; cat "$root/revision.txt"; } | tee ci-out/toolchain.txt
export CMAKE_TOOLCHAIN_FILE=$src/ci/cluster-toolchain.cmake
export RUNNER_TEMP=$root/tmp
mkdir -p "$RUNNER_TEMP"
load_job_environment() {
    while IFS='=' read -r key value; do
        [[ $key =~ ^[A-Za-z_][A-Za-z0-9_]*$ ]] && export "$key=$value"
    done < "$GITHUB_ENV"
    while IFS= read -r item; do [ -z "$item" ] || export PATH="$item:$PATH"; done < "$GITHUB_PATH"
}
case "$tier" in
  ci-gfortran)
    export OMP_NUM_THREADS="1"
    export OPENBLAS_NUM_THREADS="1"
    export MKL_NUM_THREADS="1"
    export BLIS_NUM_THREADS="1"
    echo Configure
    (
cmake -S . -B build -G Ninja

    )
    load_job_environment
    echo Build
    (
cmake --build build

    )
    load_job_environment
    echo Test
    (
ctest --test-dir build --output-on-failure

    )
    load_job_environment
    ;;
  ci-sparse)
    export OMP_NUM_THREADS="1"
    export OPENBLAS_NUM_THREADS="1"
    export MKL_NUM_THREADS="1"
    export BLIS_NUM_THREADS="1"
    echo 'Configure native sparse consumer'
    (
cmake -S . -B build -G Ninja -DGLISS_ENABLE_FORTSPARSE=ON \
  -DGLISS_REQUIRE_OPENBLAS=ON


    )
    load_job_environment
    echo 'Build and test sparse path'
    (
cmake --build build --parallel 2
ctest --test-dir build --output-on-failure


    )
    load_job_environment
    ;;
  ci-gfortran-runtime)
    export OMP_NUM_THREADS="1"
    export OPENBLAS_NUM_THREADS="1"
    export MKL_NUM_THREADS="1"
    export BLIS_NUM_THREADS="1"
    echo 'Configure with runtime checks'
    (
cmake -S . -B build -G Ninja \
  -DCMAKE_Fortran_COMPILER=gfortran \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_Fortran_FLAGS_DEBUG="-O1 -g -fcheck=all -fbacktrace" \
  -DGLISS_REQUIRE_OPENBLAS=ON


    )
    load_job_environment
    echo Build
    (
cmake --build build --parallel 2

    )
    load_job_environment
    echo 'Test with runtime checks'
    (
ctest --test-dir build --output-on-failure

    )
    load_job_environment
    ;;
  ci-lint)
    export OMP_NUM_THREADS="1"
    export OPENBLAS_NUM_THREADS="1"
    export MKL_NUM_THREADS="1"
    export BLIS_NUM_THREADS="1"
    echo 'Configure with warnings as errors'
    (
cmake -S . -B build -G Ninja -DCMAKE_Fortran_FLAGS=-Werror

    )
    load_job_environment
    echo Build
    (
cmake --build build

    )
    load_job_environment
    ;;
  ci-python)
    python3 -m venv "$root/venv"
    . "$root/venv/bin/activate"
    export OMP_NUM_THREADS="1"
    export OPENBLAS_NUM_THREADS="1"
    export MKL_NUM_THREADS="1"
    export BLIS_NUM_THREADS="1"
    echo 'Install dependencies'
    (
python -m pip install --upgrade pip
pip install build "numpy<2" pytest scipy booz-xform==0.0.9

    )
    load_job_environment
    echo 'Build shared library'
    (
cmake -S . -B build -G Ninja
cmake --build build --target gliss_c gliss_wheel_fixture
build/gliss_wheel_fixture build/installed_wheel_fixture.nc


    )
    load_job_environment
    echo 'Run Python binding tests'
    (
GLISS_LIB=$PWD/build/libgliss_c.so python -m pytest python/tests -v

    )
    load_job_environment
    echo 'Check the public analytic equilibrium generator'
    (
python benchmarks/solovev/public/test_solovev.py

    )
    load_job_environment
    echo 'Build distributions'
    (
python -m build

    )
    load_job_environment
    echo 'Test installed wheel'
    (
python -m venv /tmp/gliss-wheel
/tmp/gliss-wheel/bin/pip install "$(echo dist/*.whl)[vmec]" pytest
/tmp/gliss-wheel/bin/pip check
/tmp/gliss-wheel/bin/python ci/check_installed.py \
  --tests python/tests --equilibrium build/installed_wheel_fixture.nc


    )
    load_job_environment
    echo 'Test installed sdist'
    (
python -m venv /tmp/gliss-sdist
/tmp/gliss-sdist/bin/pip install dist/gliss-*.tar.gz
/tmp/gliss-sdist/bin/pip check
cd /tmp
/tmp/gliss-sdist/bin/python - <<'PY'
from pathlib import Path
import gliss

assert gliss.version() == gliss.__version__ == "0.0.3"
assert (Path(gliss.get_include()) / "gliss.h").is_file()
PY


    )
    load_job_environment
    ;;
  *) echo "unknown CI tier: $tier" >&2; exit 2 ;;
esac
