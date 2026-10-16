#!/usr/bin/env python3
"""Prepare the pinned N2-12 geometry, or the explicitly weak QUASR control."""

# Parse CLI before optional scientific imports; set thread limits before loading BLAS.
# ruff: noqa: E402
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("gvec_source", type=Path)
parser.add_argument("output", type=Path)
parser.add_argument("--quasr-json", type=Path, help="optional weak-shape control")
parser.add_argument("--resolution", type=int, nargs=3, metavar=("M", "N", "CELLS"))
parser.add_argument("--quadrature", type=int, nargs=2, metavar=("THETA", "ZETA"))
args = parser.parse_args()
protocol = json.loads(Path(__file__).with_name("protocol.json").read_text())
source, output = args.gvec_source.resolve(), args.output.resolve()
quasr_json = args.quasr_json.resolve() if args.quasr_json else None
revision = subprocess.check_output(
    ["git", "-C", str(source), "rev-parse", "HEAD"], text=True
).strip()
if revision != protocol["gvec_commit"]:
    parser.error("GVEC source checkout does not match protocol.json")
if output.exists():
    parser.error("output directory must be fresh")
if args.resolution and min(args.resolution) < 1:
    parser.error("resolution values must be positive")
if args.quadrature and min(args.quadrature) < 1:
    parser.error("quadrature values must be positive")
if args.resolution and not args.quadrature:
    parser.error("an explicit quadrature is required with a custom resolution")


def verify(path, expected):
    if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
        raise ValueError(f"Pinned input hash mismatch: {path}")


import gvec

output.mkdir(parents=True)
if quasr_json:
    candidate = protocol["weak_control"]
    verify(quasr_json, candidate["sha256"])
    from gvec.scripts.quasr import load_quasr

    with gvec.util.chdir(output):
        params, _ = load_quasr(
            quasr_json,
            filetype="json",
            tol=1e-7,
            nt=81,
            nz=81,
            clean=0,
            stellsym=False,
            cutoff=-1,
            name="quasr-0189686",
            quiet=False,
        )
    frame = output / "quasr-0189686-Gframe.nc"
else:
    candidate = protocol["primary"]
    parameter, frame_input = (
        source / candidate["parameter"],
        source / candidate["frame"],
    )
    verify(parameter, candidate["parameter_sha256"])
    verify(frame_input, candidate["frame_sha256"])
    shutil.copy2(parameter, output / "parameter.upstream.ini")
    frame = output / frame_input.name
    shutil.copy2(frame_input, frame)
    params = gvec.util.read_parameters(parameter)
    # Restore the paper boundary's N=10 support; CI's N=5 truncates the boundary.
    args.resolution = args.resolution or [2, 10, 2]
    args.quadrature = args.quadrature or [9, 49]
if args.resolution:
    if not args.quadrature:
        parser.error("an explicit quadrature is required with a custom resolution")
    m, n, cells = args.resolution
    for name in ["X1_mn_max", "X2_mn_max", "LA_mn_max"]:
        params[name] = [m, n]
    params["sgrid"]["nElems"] = cells
if args.quadrature:
    params.update(fac_nyq=-1, mn_nyq=args.quadrature)
params.update(
    hmap_ncfile=str(frame),
    boundary_filename=str(frame),
    which_visu=0,
    visu1D=0,
    visu2D=0,
    visu3D=0,
)
gvec.util.write_parameters(params, output / "parameters.toml")
manifest = {
    "status": "prepared-unqualified",
    "candidate": candidate,
    "source_checkout_commit": revision,
    "gvec_package_version": gvec.__version__,
    "installed_source_identity": "Build from the pinned checkout using README; version alone is not attestation.",
    "omp_threads": 1,
    "resolution": params["X1_mn_max"],
    "radial_grid": params["sgrid"],
    "sha256": {
        p.name: hashlib.sha256(p.read_bytes()).hexdigest()
        for p in [frame, output / "parameters.toml"]
    },
}
(output / "manifest.json").write_text(
    json.dumps(gvec.util.serialize(manifest), indent=2) + "\n"
)
print(output / "parameters.toml")
