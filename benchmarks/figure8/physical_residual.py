"""Off-knot physical residual oracle independent of minimizer force norms."""

from pathlib import Path
import hashlib
import argparse
import json
import os

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"

import numpy as np
import xarray as xr
from numpy.polynomial.legendre import leggauss
import gvec

parser = argparse.ArgumentParser(description=__doc__)
for name in ["parameter", "state", "output"]:
    parser.add_argument(name, type=Path)
parser.add_argument("--toroidal-multiplier", type=int, default=1)
parser.add_argument(
    "--force-free",
    action="store_true",
    help="fixed-iota p=0 case; do not require zero current",
)
args = parser.parse_args()
if args.toroidal_multiplier < 1:
    parser.error("toroidal multiplier must be positive")
goal_path = Path(__file__).with_name("acceptance-goals.json")
goals = json.loads(goal_path.read_text())
state = gvec.State(args.parameter.resolve(), args.state.resolve())


def norm(vector):
    return np.sqrt((vector * vector).sum("xyz"))


def cartesian_curl(ev):
    return sum(
        xr.cross(ev[gradient], ev[derivative], dim="xyz")
        for gradient, derivative in [
            ("grad_rho", "dB_dr"),
            ("grad_theta", "dB_dt"),
            ("grad_zeta", "dB_dz"),
        ]
    )


levels = [(6, 33, 49), (8, 49, 65), (10, 65, 81)]
physical = []
for nq, nt, nz in levels:
    if args.toroidal_multiplier > 1:
        nz = nz * args.toroidal_multiplier + 1
    x, w = leggauss(nq)
    edges = state.get_radial_gridpoints()
    theta = 2 * np.pi * (np.arange(nt) + 0.371) / nt
    zeta = 2 * np.pi * (np.arange(nz) + 0.293) / (nz * state.nfp)
    sums = dict(
        volume=0.0, force_l1=0.0, force_l2=0.0, scale=0.0, current_l2=0.0, field_l2=0.0
    )
    max_force = 0.0
    min_jac = np.inf
    max_curl_disagreement = 0.0
    max_pressure_slope = 0.0
    for a, b in zip(edges[:-1], edges[1:]):
        rho = (a + b) / 2 + (b - a) * x / 2
        wr = (b - a) * w / 2
        ev = state.evaluate(
            "Jac",
            "B",
            "J",
            "dp_dr",
            "grad_rho",
            "grad_theta",
            "grad_zeta",
            "dB_dr",
            "dB_dt",
            "dB_dz",
            "mod_B",
            "grad_mod_B",
            rho=rho,
            theta=theta,
            zeta=zeta,
            extrapolate=False,
        )
        weights = xr.DataArray(wr, dims="rad") * ev.Jac / (nt * nz)
        independent_current = cartesian_curl(ev) / ev.mu0
        independent_force = (
            xr.cross(independent_current, ev.B, dim="xyz") - ev.dp_dr * ev.grad_rho
        )
        force_norm = norm(independent_force)
        scale = ev.mod_B * norm(ev.grad_mod_B) / ev.mu0
        if not all(np.isfinite(q).all() for q in [weights, force_norm, scale]):
            raise ValueError("Nonfinite independent residual or normalization")
        sums["volume"] += float(weights.sum())
        sums["force_l1"] += float((weights * force_norm).sum())
        sums["force_l2"] += float((weights * force_norm**2).sum())
        sums["scale"] += float((weights * scale).sum())
        sums["current_l2"] += float((weights * norm(independent_current) ** 2).sum())
        sums["field_l2"] += float((weights * ev.mod_B**2).sum())
        max_force = max(max_force, float(force_norm.max()))
        min_jac = min(min_jac, float(ev.Jac.min()))
        max_curl_disagreement = max(
            max_curl_disagreement,
            float(norm(independent_current - ev.J).max()),
        )
        max_pressure_slope = max(max_pressure_slope, float(abs(ev.dp_dr).max()))
    normalization = sums["scale"] / sums["volume"]
    if not np.isfinite(normalization) or normalization <= 0 or min_jac <= 0:
        raise ValueError("Invalid physical normalization or off-axis Jacobian")
    result = {
        "quadrature_per_cell_theta_zeta": [nq, nt, nz],
        "volume_m3": sums["volume"] * (2 * np.pi) ** 2,
        "normalization_N_per_m3": normalization,
        "force_relative_L1": sums["force_l1"] / sums["scale"],
        "force_relative_RMS": np.sqrt(sums["force_l2"] / sums["volume"])
        / normalization,
        "force_relative_max": max_force / normalization,
        "jacobian_min": min_jac,
        "cartesian_vs_metric_current_max_A_per_m2": max_curl_disagreement,
        "pressure_slope_max": max_pressure_slope,
        "vacuum_current_relative_RMS": (
            (sums["volume"] * (2 * np.pi) ** 2) ** (1 / 3)
            * float(ev.mu0)
            * np.sqrt(sums["current_l2"] / sums["field_l2"])
        ),
    }
    physical.append(result)
    print(json.dumps(result), flush=True)

coordinates = {
    "rho": np.array([0.137, 0.379, 0.631, 0.893]),
    "theta": np.array([0.127, 0.763, 1.923]),
    "zeta": np.array([0.139, 0.613, 1.317]) / state.nfp,
}
ev = state.evaluate(
    "B",
    "J",
    "dp_dr",
    "grad_rho",
    "grad_theta",
    "grad_zeta",
    "dB_dr",
    "dB_dt",
    "dB_dz",
    **coordinates,
    extrapolate=False,
)
finite_difference = []
for step in [3e-4, 1e-4, 3e-5, 1e-5]:
    curl_fd = 0
    for coordinate, gradient in [
        ("rho", "grad_rho"),
        ("theta", "grad_theta"),
        ("zeta", "grad_zeta"),
    ]:
        plus, minus = coordinates.copy(), coordinates.copy()
        plus[coordinate] = coordinates[coordinate] + step
        minus[coordinate] = coordinates[coordinate] - step
        bp = state.evaluate("B", **plus, extrapolate=False).B
        bm = state.evaluate("B", **minus, extrapolate=False).B
        # Coordinates label shifted points; align by position, not xarray labels.
        slope = xr.DataArray(
            (bp.values - bm.values) / (2 * step),
            dims=ev.B.dims,
            coords=ev.B.coords,
        )
        curl_fd += xr.cross(ev[gradient], slope, dim="xyz")
    fd_force = xr.cross(curl_fd / ev.mu0, ev.B, dim="xyz") - ev.dp_dr * ev.grad_rho
    mismatch_force = xr.cross((curl_fd - cartesian_curl(ev)) / ev.mu0, ev.B, dim="xyz")
    result = {
        "central_difference_step": step,
        "force_relative_point_RMS": float(np.sqrt((norm(fd_force) ** 2).mean()))
        / normalization,
        "analytic_vs_difference_force_relative_RMS": float(
            np.sqrt((norm(mismatch_force) ** 2).mean())
        )
        / normalization,
    }
    finite_difference.append(result)
    print(json.dumps(result), flush=True)

diagnostics = {
    "status": "independent residual evidence; equilibrium qualification remains incomplete",
    "gvec_package_version": gvec.__version__,
    "input_sha256": {
        name: hashlib.sha256(path.read_bytes()).hexdigest()
        for name, path in [("parameter", args.parameter), ("state", args.state)]
    },
    "acceptance_goals_sha256": hashlib.sha256(goal_path.read_bytes()).hexdigest(),
    "omp_threads": int(os.environ["OMP_NUM_THREADS"]),
    "blas_thread_environment": int(os.environ["OPENBLAS_NUM_THREADS"]),
    "physical_quadrature_refinement": physical,
    "finite_difference_step_sweep": finite_difference,
    "physical_force_RMS_goal_passed": bool(
        physical[-1]["force_relative_RMS"] <= goals["physical_force_relative_rms"]
    ),
    "vacuum_current_RMS_goal_passed": None
    if args.force_free
    else (
        bool(
            physical[-1]["vacuum_current_relative_RMS"]
            <= goals["vacuum_current_relative_rms"]
        )
    ),
    "equilibrium_resolution_convergence_qualified": False,
    "embedding_and_nestedness_qualified": False,
    "export_qualified": False,
}
args.output.write_text(json.dumps(diagnostics, indent=2) + "\n")

raise SystemExit(
    0
    if diagnostics["physical_force_RMS_goal_passed"]
    and diagnostics["vacuum_current_RMS_goal_passed"] is not False
    else 2
)
