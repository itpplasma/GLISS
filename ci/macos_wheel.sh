#!/usr/bin/env bash
# Build, test and repair the macOS wheel for one architecture inside the locked
# micromamba toolchain (ci/macos-<arch>.lock); see ci/macos-wheels.md.
set -euo pipefail
architecture=${1:?usage: ci/macos_wheel.sh arm64|x86_64}
test "$(uname -m)" = "$architecture"
export MACOSX_DEPLOYMENT_TARGET=15.0
SDKROOT="$(xcrun --sdk macosx --show-sdk-path)"
export SDKROOT
export CMAKE_PREFIX_PATH="$CONDA_PREFIX"
export PKG_CONFIG_PATH="$CONDA_PREFIX/lib/pkgconfig"
python -m pip install build==1.3.0 delocate==0.13.0 twine==6.2.0 scipy
mkdir -p artifacts
micromamba list --explicit > artifacts/toolchain.lock
python -m pip freeze > artifacts/python-build.txt
cp pyproject.toml artifacts/build-metadata.toml
"$FC" --version > artifacts/fortran.txt
sw_vers > artifacts/macos.txt
xcodebuild -version > artifacts/xcode.txt
cmake -S . -B build -G Ninja -DGLISS_REQUIRE_OPENBLAS=ON
cmake --build build --parallel 3
ctest --test-dir build --output-on-failure
build/gliss_wheel_fixture artifacts/installed_wheel_fixture.nc
python -m build --wheel
delocate-wheel --require-archs "$architecture" \
  --sanitize-rpaths -w artifacts dist/*.whl
delocate-listdeps --all artifacts/*.whl > artifacts/dependencies.txt
python -m twine check artifacts/*.whl
