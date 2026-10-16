#!/usr/bin/env python3
"""Execute one bounded serial GVEC attempt; never promote its final state."""

# Parse CLI before optional scientific imports; set thread limits before loading BLAS.
# ruff: noqa: E402
import argparse
import hashlib
import json
import os
from pathlib import Path
import time

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("parameter", type=Path)
parser.add_argument("output", type=Path)
parser.add_argument("--iterations", type=int, default=20000)
parser.add_argument("--restart-parameter", type=Path)
parser.add_argument("--restart-state", type=Path)
parser.add_argument(
    "--target-frame", type=Path, help="axisNB source; defaults to boundary_filename"
)
args = parser.parse_args()
if args.iterations < 1 or args.output.exists():
    parser.error("use a positive budget and a fresh output directory")
if bool(args.restart_parameter) != bool(args.restart_state):
    parser.error("provide both restart paths together")
import numpy as np
import gvec
from boundary import check_boundary

params = gvec.util.read_parameters(args.parameter.resolve())
for name in ["hmap_ncfile", "boundary_filename"]:
    if name in params:
        path = Path(params[name])
        params[name] = str(
            path if path.is_absolute() else args.parameter.resolve().parent / path
        )
frame = (
    args.target_frame.resolve()
    if args.target_frame
    else Path(params["boundary_filename"])
)
goal_path = Path(__file__).with_name("acceptance-goals.json")
goals = json.loads(goal_path.read_text())
params.update(
    totalIter=args.iterations,
    maxIter=args.iterations,
    MinimizerType=10,
    minimize_tol=goals["minimizer_max_force"],
    outputIter=1000,
    logIter=100,
)
if "I_tor" in params:
    params["maxIter"] = min(1000, args.iterations)
restart = (
    gvec.State(args.restart_parameter.resolve(), args.restart_state.resolve())
    if args.restart_state
    else None
)
if restart is not None:
    edge = check_boundary(restart, frame, goals["boundary_position_max_error_m"])
    if not edge["boundary_goal_passed"]:
        parser.error(f"restart retains a different boundary: {edge}; initialize fresh")
args.output.resolve().mkdir(parents=True)
started = time.monotonic()
run = gvec.run(
    params,
    restartstate=restart,
    runpath=args.output.resolve(),
    keep_intermediates="stages",
    quiet=False,
    loglevel="INFO",
)
current = getattr(run, "rms_iota", None)
fit = getattr(run, "rms_iota_fit", None)
current = float(current) if current is not None else None
fit = float(fit) if fit is not None else None
current_ok = (
    None
    if "I_tor" not in params
    else bool(
        current is not None
        and fit is not None
        and np.isfinite(current)
        and np.isfinite(fit)
        and current <= goals["picard_rms_delta_iota"]
        and fit <= goals["picard_rms_iota_fit"]
    )
)
force_ok = bool(
    np.isfinite(run.max_force) and run.max_force <= goals["minimizer_max_force"]
)
edge = check_boundary(run.state, frame, goals["boundary_position_max_error_m"])
passed = force_ok and current_ok is not False and edge["boundary_goal_passed"]
report = {
    "status": "solver-goals-passed" if passed else "nonconverged-or-boundary-mismatch",
    "iterations_used": int(run.GVEC_iter_used),
    "elapsed_seconds": time.monotonic() - started,
    "max_force": float(run.max_force),
    "rms_delta_iota": current,
    "rms_iota_fit": fit,
    "minimizer_converged": force_ok,
    "current_picard_converged": current_ok,
    "parameterfile": str(run.state.parameterfile),
    "statefile": str(run.state.statefile),
    "gvec_package_version": gvec.__version__,
    "omp_threads": 1,
    "acceptance_goals_sha256": hashlib.sha256(goal_path.read_bytes()).hexdigest(),
    "boundary": edge,
    "physical_equilibrium_qualified": False,
}
(args.output / "minimizer-report.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report, indent=2))
raise SystemExit(0 if passed else 2)
