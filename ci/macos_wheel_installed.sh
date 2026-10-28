#!/usr/bin/env bash
# Install the repaired wheel from artifacts/ into a fresh environment of the
# python on PATH and test it outside the checkout; see ci/macos-wheels.md.
set -euo pipefail
architecture=${1:?usage: ci/macos_wheel_installed.sh arm64|x86_64 [python-version]}
workspace=$(pwd)
scratch=$(mktemp -d "${RUNNER_TEMP:-/tmp}/gliss-wheel.XXXXXX")
if [ -n "${2:-}" ]; then
  # Outside a workflow step (nightly bisection), select the runner's Python here.
  GITHUB_PATH=$scratch/path GITHUB_ENV=$scratch/env GITHUB_OUTPUT= \
    /opt/itp-ci-macos/setup-python.sh "$2"
  PATH="$(cat "$scratch/path"):$PATH"
fi
unset GLISS_LIB PYTHONPATH
for loader_variable in ${!DYLD_@}; do unset "$loader_variable"; done
python -m venv "$scratch/venv"
python_exe="$scratch/venv/bin/python"
"$python_exe" -m pip install artifacts/*.whl pytest scipy
"$python_exe" -m pip check
cd "$scratch"
"$python_exe" "$workspace/ci/check_macos_wheel.py" \
  --architecture "$architecture" \
  --c-consumer "$workspace/cmake/test_gliss_header.c"
"$python_exe" "$workspace/ci/check_installed.py" \
  --tests "$workspace/python/tests" \
  --equilibrium "$workspace/artifacts/installed_wheel_fixture.nc"
