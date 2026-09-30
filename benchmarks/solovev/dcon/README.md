# DCON reference for the Solov'ev n=1 benchmark

`run.sh OUTPUT [Q0 ...]` fetches GPEC at the pinned commit
`f5595c0689c624834d720068ddc8b5a7e4028248` from
<https://github.com/PrincetonUniversity/GPEC>, builds `dcon` with gfortran,
OpenBLAS and NetCDF, and counts fixed-boundary Newcomb zero crossings for the
analytic Solov'ev equilibrium (`docs/examples/solovev_ideal_example`) with
three switches changed: `vac_flag=f`, `qlow=0.5`, `out_fixed=t`. These are the
settings of GPEC's own `regression/test/fixed_boundary_diagnostics.py`.
Building takes a few minutes; each q0 takes about one second.

| q0 | 1.020 | 1.025 | 1.030 | 1.032 | 1.035 | 1.037 | 1.039062 | 1.0394 | 1.0397 | 1.039843 | 1.045 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| example tolerances | 1 | 1 | **0** | 1 | **0** | 1 | 1 | 1 | 0 | 0 | 0 |
| `tol_nr=tol_r=1e-8` | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 0 | 0 | 0 |

The marginal point lies in (1.0394, 1.0397), consistent with GLISS (1.039062
unstable, 1.039843 stable) and the converged DCON sensitivity study
(1.03956-1.03959).

## Upstream defect

With the example tolerances (`tol_nr=1e-6`, `tol_r=1e-7`), DCON misses the
Newcomb crossing at isolated q0 values such as 1.030 and 1.035, although both
neighbours count one. `ode_output_monitor` (`dcon/ode_output.f`) takes the sign
change of `crit`, the dominant eigenvalue of the inverse response matrix, as a
zero crossing only if a linearly interpolated midpoint lies between the
endpoint values *and* has less than half their smaller magnitude. Over coarse
output steps `crit` is strongly nonlinear. At q0=1.035 a genuine root
(-0.050 -> +0.173 between psi 0.670 and 0.683) interpolates to +0.163 and is
rejected. Dropping the magnitude test instead admits poles: at q0=1.037 the
sign change 2.615 -> -2.515 interpolates between its endpoints. The
classification needs a step-size independent test, for example bisection of
the interpolated interval until |crit| clearly tends to zero or to infinity.
Tighter tolerances hide the defect here, but they do not remove it. GLISS
tracks the report in its issue tracker; it has not been filed upstream yet.
