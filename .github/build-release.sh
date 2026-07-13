#!/usr/bin/env bash
set -euo pipefail

root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$root"

tag=$(git describe --exact-match --match 'v[0-9]*' HEAD)
version=${tag#v}
epoch=$(git show -s --format=%ct HEAD)
image="gliss-manylinux-release:$version"
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

if [[ -n $(git status --porcelain) ]]; then
    echo "release build requires a clean worktree" >&2
    exit 1
fi
if ! rg -q "^version = \"$version\"$" pyproject.toml; then
    echo "tag $tag does not match pyproject.toml" >&2
    exit 1
fi

mkdir -p "$tmp/sdist-a" "$tmp/sdist-b" dist
rm -f dist/*
for output in "$tmp/sdist-a" "$tmp/sdist-b"; do
    SOURCE_DATE_EPOCH=$epoch uv tool run --from build pyproject-build \
        --sdist --outdir "$output"
done
cmp "$tmp/sdist-a/gliss-$version.tar.gz" \
    "$tmp/sdist-b/gliss-$version.tar.gz"
cp "$tmp/sdist-a/gliss-$version.tar.gz" dist/

docker build -f .github/manylinux.Dockerfile -t "$image" .github

build_wheel() {
    local output=$1
    mkdir -p "$output/raw" "$output/repaired"
    docker run --rm -e SOURCE_DATE_EPOCH="$epoch" \
        -v "$root/dist:/input:ro" -v "$output/raw:/output" "$image" \
        bash -lc "mkdir /src && tar -xzf /input/gliss-$version.tar.gz \
        --strip-components=1 -C /src && \
        /opt/python/cp39-cp39/bin/python -m pip wheel /src --no-deps \
        --wheel-dir /output"
    docker run --rm -e SOURCE_DATE_EPOCH="$epoch" \
        -v "$output/raw:/input:ro" -v "$output/repaired:/output" "$image" \
        bash -lc "/opt/python/cp39-cp39/bin/python -m auditwheel repair \
        --plat manylinux_2_28_x86_64 --wheel-dir /output \
        /input/gliss-$version-py3-none-linux_x86_64.whl"
}

build_wheel "$tmp/wheel-a"
build_wheel "$tmp/wheel-b"
wheel_a=$(find "$tmp/wheel-a/repaired" -name '*.whl' -print -quit)
wheel_b=$(find "$tmp/wheel-b/repaired" -name '*.whl' -print -quit)
cmp "$wheel_a" "$wheel_b"
cp "$wheel_a" dist/

uv tool run --from twine twine check dist/*
sha256sum dist/*
