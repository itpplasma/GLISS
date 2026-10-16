#!/usr/bin/env python3
"""Plot unmodified GVEC input frames on shared Cartesian scales; no solve."""

# Parse CLI before optional scientific imports; set thread limits before loading BLAS.
# ruff: noqa: E402
import argparse
import hashlib
import json
import os
from pathlib import Path

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "output", type=Path, help="PNG or PDF; adjacent JSON records inputs"
)
parser.add_argument("frames", type=Path, nargs="+")
parser.add_argument("--labels", nargs="+")
args = parser.parse_args()
if args.labels and len(args.labels) != len(args.frames):
    parser.error("provide one label per frame")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from gvec.gframe import read_Gframe_ncfile, to_axis, to_surface

frames, manifest = [], []
for index, path in enumerate(args.frames):
    data = read_Gframe_ncfile(path)
    axis = to_axis(data, nzeta=201)["axis"]["xyz"].T
    boundary = to_surface(data, nzeta=201, ntheta=97)["xyz"]
    label = args.labels[index] if args.labels else path.name
    frames.append((label, np.vstack([axis, axis[0]]), boundary))
    manifest.append(
        {
            "label": label,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "nfp": int(data["nfp"]),
            "full_turn_axis_samples": len(axis),
            "theta_samples": 97,
        }
    )
points = np.concatenate([surface.reshape(-1, 3) for _, _, surface in frames])
lo, hi = points.min(axis=0), points.max(axis=0)
padding = 0.05 * (hi - lo)
lo, hi = lo - padding, hi + padding
fig = plt.figure(figsize=(16, 3.7 * len(frames)), layout="constrained")
for row, (label, axis, surface) in enumerate(frames):
    for col, (i, j) in enumerate([(0, 1), (0, 2), (1, 2)]):
        ax = fig.add_subplot(len(frames), 4, 4 * row + col + 1)
        for k in range(0, surface.shape[1], 12):
            curve = np.vstack([surface[:, k], surface[0, k]])
            ax.plot(curve[:, i], curve[:, j], color="#56B4E9", lw=0.6, alpha=0.5)
        for k in range(0, surface.shape[0], 20):
            curve = np.vstack([surface[k], surface[k, 0]])
            ax.plot(curve[:, i], curve[:, j], color="#56B4E9", lw=0.5, alpha=0.4)
        ax.plot(axis[:, i], axis[:, j], color="black", lw=1.6)
        ax.set(
            xlim=(lo[i], hi[i]),
            ylim=(lo[j], hi[j]),
            xlabel="xyz"[i] + " [m]",
            ylabel="xyz"[j] + " [m]",
            title=label + "\n" + "xyz"[i] + "xyz"[j] + " projection",
        )
        ax.set_aspect("equal", adjustable="box")
        ax.grid(alpha=0.2)
    ax = fig.add_subplot(len(frames), 4, 4 * row + 4, projection="3d")
    ax.plot_wireframe(
        *surface.transpose(2, 0, 1),
        rstride=20,
        cstride=12,
        color="#56B4E9",
        lw=0.5,
        alpha=0.5,
    )
    ax.plot(*axis.T, color="black", lw=1.6)
    ax.set(
        xlim=(lo[0], hi[0]),
        ylim=(lo[1], hi[1]),
        zlim=(lo[2], hi[2]),
        xlabel="x [m]",
        ylabel="y [m]",
        zlabel="z [m]",
        title="Full 3D input geometry",
    )
    ax.set_box_aspect(hi - lo)
    ax.view_init(elev=25, azim=-55)
fig.suptitle(
    "GVEC input frames on shared Cartesian scales\nBlue: input boundary; black: frame centreline. Geometry alone does not qualify equilibrium."
)
fig.savefig(args.output, dpi=150, bbox_inches="tight")
args.output.with_suffix(".json").write_text(
    json.dumps(
        {
            "method": "GVEC frame interpolation; no centering, rescaling, rotation or filtering",
            "shared_cartesian_limits_m": [lo.tolist(), hi.tolist()],
            "frames": manifest,
        },
        indent=2,
    )
    + "\n"
)
print(args.output)
