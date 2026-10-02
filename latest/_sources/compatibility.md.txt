# TERPSICHORE compatibility

`solve_terpsichore_fixed_boundary` and `solve_terpsichore_pseudoplasma`
replay the matrices stored in TERPSICHORE FORT.23 (and FORT.24) files and
certify the lowest eigenpair with GLISS's inertia bracketing. They
reproduce TERPSICHORE's eigenvalue to 1e-8 relative on the QAS3 family
(`benchmarks/qas3`). A stable file returns `negative_count == 0` with the
lowest nonnegative eigenpair.

These replays validate GLISS against TERPSICHORE's own discretization and
normalization. They are not GLISS's physical fixed-boundary operator and not
a production plasma-vacuum solve.
