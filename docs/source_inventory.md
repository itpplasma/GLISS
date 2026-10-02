# Source inventory

The machine-readable {download}`source inventory <../source_inventory.jsonl>`
accounts for the native library, Python interfaces, C declarations, applications,
benchmark runners, verification sources, and build/CI configuration. Each JSON
line identifies one source file, its SHA-256, and the exact source snapshot base.
The inventory is a starting point for assigning equation and interface reviews.
It does not establish that those reviews have been completed.

The `mathematical_contract`, `units`, `assumptions`, `independent_oracle`,
`derivative_coverage`, and `reviewer` fields explicitly retain unmapped or
unreviewed whole-component obligations. `targeted_contract_evidence` records
only bounded fixes and oracle descriptions at cited commits. A review of a
single resonance derivative or inertia probe cannot qualify the whole module.
The recorded regressions are evidence sources, rather than new executions by
the inventory generator. Subsequent edits require a fresh scoped review.

## Reference map

`dependencies` and `caller_references` are reciprocal static references. The
generator resolves Fortran `use` statements, Python imports, and C includes to included
source files. It also records C ABI token references and literal file paths
with their source line numbers. Those latter references can occur in comments
or strings; they identify candidates for review rather than prove execution.
Python package re-exports, dynamic imports/calls, runtime compiler selection,
and generated command paths still require manual tracing. External and
unresolved imports remain visible in each record.
Fortran interface declarations and named external C bindings also remain
visible as external candidates when no included binding definition is found;
their library linkage needs review. Native named C binding definitions are
recorded separately, including benchmark bindings whose names lack `gliss_`.

The native-library and header records list each detected named
`gliss_*` entry point in the source. This does not establish compiled symbol
visibility. Their individual contracts remain explicitly unmapped
unless a bounded evidence record applies. Binary equilibrium fixtures and
external solver/data repositories are outside this source-file inventory.
Their identities and benchmark acceptance evidence remain separate obligations
in [PLAN.md](https://github.com/itpplasma/GLISS/blob/main/PLAN.md) and
[PROVENANCE.md](https://github.com/itpplasma/GLISS/blob/main/PROVENANCE.md).

## Updating the inventory

Run these commands from the repository root after changes to included sources:

```sh
python3 ci/source_inventory.py
python3 ci/source_inventory.py --check
```

Writing pins the current Git commit and hashes the current included sources;
an uncommitted source patch therefore remains distinguishable by its hashes.
`--base COMMIT` can retain a deliberately frozen review base. `--check` checks
only metadata consistency with the current files and static reference map. It
does not test scientific behavior, oracle independence, derivative correctness,
or source-review completion. The generator writes and flushes JSON lines in
bounded batches so an interrupted write leaves recoverable partial work.

The evidence catalog in `ci/source_inventory_evidence.json` must preserve the
review scope, cited fix commit, independent oracle, derivative scope, reviewer,
and remaining gaps. New evidence must describe an actually performed bounded
review; test counts alone do not establish equation coverage.
