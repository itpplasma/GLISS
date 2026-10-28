# Contributing

## Merging

- Pull requests, including every layer of a stacked PR, run the fast CI tier.
  Label a PR `full-ci` to run the full tier on it before queueing.
- Review a stack layer by layer. When it is ready, collapse it into one PR into
  `main` by retargeting the top layer to `main`; the lower PRs then show as
  merged.
- Click "Merge when ready" (or run `gh pr merge --merge --auto`). The PR enters
  the merge queue, the full CI tier runs on the exact merge result, and on
  success it lands on `main` as a merge commit. `main` itself is not re-tested
  on push.
- A branch need not be up to date with `main`. Merge `main` into a PR only when
  it conflicts.

Before a release, run the array-temporary audit described in
[docs/building.md](docs/building.md).
