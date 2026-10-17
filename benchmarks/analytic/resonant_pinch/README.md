# Unstable resonant screw-pinch eigenvalue

The registered `resonant_screw_pinch_spectrum` test checks the most-negative
isolated eigenvalue of the exact equilibrium
`B_theta=A*r`, `B_z=Z`, `p=p_edge+A^2*(a^2-r^2)/mu0`, at `k=-m*A/Z`.
Inputs are a=0.7 m, A=0.8 T/m, Z=1.2 T, p_edge=17000 Pa,
rho=1.7e-7 kg/m^3, m=n=2, L=2*pi*Z/A, and gamma=1 or 5/3.
The axis branch is Cartesian regular and the rigid wall imposes f(a)=0.

The [exact derivation](symbolic/derivation.md) starts from the Bernstein
energy principle and cylindrical vector calculus. Its reproducible FortSym
driver proves 45 identities with native FortSym and SymEngine: 90 exact
checks, zero probes. At nonzero eigenvalue the unrestricted variation
forces perpendicular displacement. The independent radial ODE supplies
references in physical units of s^-2:

| gamma | Ground eigenvalue | Dimensionless ell |
| --- | ---: | ---: |
| 1 | -3.34986873378943e11 | -0.02435118850778054 |
| 5/3 | -3.11741695965623e11 | -0.02266142767811291 |

Here `ell=lambda*mu0*rho*a^2/Z^2`. Continuous Pruefer phase brackets the
ground root using the completed-square lower bound and an independently
negative trial Rayleigh upper bound. The first overtone is separately
identified. DOP853 axis/accuracy studies and 80-digit classical RK4 with
graded axis steps and three refinements agree to 1.52e-15 relatively.
See the frozen [oracle receipt](oracle/receipt.json) and protocols.
Its `shoot.py` and `shooting.json` hashes bind the preserved
`shoot-premeasurement.py` and `shooting-premeasurement.json` bytes. A later
[output-only correction](oracle/oracle-metadata-correction.json) fixes the
reported lower-endpoint phase after a loop reused its variable. All roots,
axis studies and other output values are exactly unchanged; the corrected
scripts and results retain their own hashes. Acceptance settings are unchanged.

The native fixture gathers the global radial pencil using the current FEEC
basis, axis factors, full sqrt(s) axis quadrature, split interior quadrature,
production surface stiffness and physical mass. Current eigenpair refinement,
mass-residual bounds and inertia brackets certify its ground mode. All 36
cases pass: degrees 1–3, cells 16/32/64, both parities and both gamma values.
Worst finest errors are 1.403e-4, 5.626e-9 and 1.476e-11, within budgets fixed
before production measurement of 1e-3, 1e-5 and 1e-7 respectively.
The [native results](native-results.json) retain every case and source hashes.

Run the registered test from a configured checkout:

```sh
cmake --build build --target test_resonant_screw_pinch_spectrum --parallel 2
ctest --test-dir build -R '^resonant_screw_pinch_spectrum$' --output-on-failure
```

To reproduce the numerical oracle, copy `oracle/shoot.py` and
`oracle/precision.py` into a fresh scratch directory. Run `shoot.py`, then
`precision.py` with Python containing NumPy, SciPy and mpmath; both scripts
write results beside themselves. The recorded environment and hashes are
in the receipt. Symbolic CMake/CTest reproduction is described separately
in the derivation and runs in pinned generated-kernel CI.

This qualifies a global radial eigenvalue calculation starting at prescribed
exact cylindrical surface fields. Native equilibrium import and public
problem orchestration are separate gates. Infinitely many unstable radial
modes approach zero alongside the parallel neutral branch; this test makes
no total negative-count or gap-at-zero claim. Marginal thresholds, other
mode families and CAS3D numerical agreement remain open.
