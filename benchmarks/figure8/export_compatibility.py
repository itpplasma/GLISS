#!/usr/bin/env python3
"""Export an unqualified state with the pinned consumer fix and check its frame."""

# Parse CLI before scientific imports; bound native threads before loading BLAS.
# ruff: noqa: E402
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
parser = argparse.ArgumentParser(description=__doc__)
for name in ["exporter", "parameter", "state", "output"]:
    parser.add_argument(name, type=Path)
args = parser.parse_args()
if args.output.exists():
    parser.error("output must be fresh")
evidence = json.loads(
    Path(__file__).with_name("compatibility-evidence.json").read_text()
)
exporter_hash = hashlib.sha256(args.exporter.read_bytes()).hexdigest()
if exporter_hash != evidence["producer"]["exporter_sha256"]:
    parser.error("exporter does not match the pinned reviewed consumer patch")
import gvec
import numpy as np

spec = importlib.util.spec_from_file_location("private_cas3d", args.exporter)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
state = gvec.State(args.parameter.resolve(), args.state.resolve())
ev = gvec.EvaluationsBoozer(
    [0.37], 25, 73, state=state, radial_derivative=False, MNfactor=3
)
state.compute(ev, "N_FP", "pos")
gvec.surface.get_xyz_hat(ev, 1)


def reconstruction_error():
    x = ev.xhat * np.cos(ev.zeta_B) - ev.yhat * np.sin(ev.zeta_B)
    y = ev.xhat * np.sin(ev.zeta_B) + ev.yhat * np.cos(ev.zeta_B)
    return float(
        np.sqrt(
            (x - ev.pos.sel(xyz="x")) ** 2
            + (y - ev.pos.sel(xyz="y")) ** 2
            + (ev.zhat - ev.pos.sel(xyz="z")) ** 2
        ).max()
    )


unpatched = reconstruction_error()
module._set_boozer_position_frame(ev, 1)
patched = reconstruction_error()
if not np.isfinite(patched) or patched > 1e-12:
    raise ValueError("Boozer position reconstruction failed before export")
module.gvec_to_cas3d(
    state, args.output.resolve(), **evidence["producer"]["export_options"]
)
report = {
    "status": "UNQUALIFIED physical equilibrium; geometry compatibility diagnostic only",
    "producer_commit": evidence["producer"]["commit"],
    "exporter_sha256": exporter_hash,
    "state_sha256": hashlib.sha256(args.state.read_bytes()).hexdigest(),
    "export_sha256": hashlib.sha256(args.output.read_bytes()).hexdigest(),
    "export_options": evidence["producer"]["export_options"],
    "native_angle_rotated_hat_reconstruction_error_m": unpatched,
    "patched_boozer_angle_rotated_hat_reconstruction_error_m": patched,
    "physical_equilibrium_qualified": False,
}
args.output.with_suffix(".json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report, indent=2))
