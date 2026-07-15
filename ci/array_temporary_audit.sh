#!/usr/bin/env bash

set -euo pipefail

repository=$(git rev-parse --show-toplevel)
if ! git -C "$repository" diff --quiet --ignore-submodules -- \
        || ! git -C "$repository" diff --cached --quiet --ignore-submodules --; then
    printf '%s\n' \
        'array-temporary audit requires a clean tracked tree; commit the code first' >&2
    exit 2
fi

temporary_base=${GLISS_AUDIT_TMPDIR:-${TMPDIR:-/tmp}}
mkdir -p "$temporary_base"
audit_worktree=$(mktemp -d "$temporary_base/gliss-array-temporary-audit.XXXXXX")
rmdir "$audit_worktree"

cleanup() {
    git -C "$repository" worktree remove --force "$audit_worktree" \
        >/dev/null 2>&1 || true
    rm -rf -- "$audit_worktree"
}
trap cleanup EXIT HUP INT TERM

git -C "$repository" worktree add --detach "$audit_worktree" HEAD
export FO_CACHE_DIR="$audit_worktree/.fo-cache"

cd "$audit_worktree"
fo clean --cache
fo build --flag '-O3 -Warray-temporaries -Werror=array-temporaries'
fo test
