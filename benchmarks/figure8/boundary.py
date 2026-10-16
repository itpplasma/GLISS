#!/usr/bin/env python3
"""Check an axisNB edge against independent interpolation of its source frame."""

import argparse
import hashlib
import json
import os
from pathlib import Path

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
import numpy as np
import xarray as xr


def interpolate(values, points, origin, period, axis):
    """Evaluate the DFT of the original samples, retaining their grid offset."""
    n = values.shape[axis]
    modes = np.fft.fftfreq(n) * n * 2 * np.pi / period
    coefficients = np.moveaxis(np.fft.fft(values, axis=axis) / n, axis, 0)
    return np.real(
        np.einsum(
            "tn,n...->t...", np.exp(1j * np.outer(points - origin, modes)), coefficients
        )
    )


def target_boundary(frame, theta, zeta):
    nfp = int(xr.load_dataset(frame)["NFP"])
    axis = xr.load_dataset(frame, group="axis")
    boundary = xr.load_dataset(frame, group="boundary")
    za = axis["zeta(:)"].values
    nz = int(axis["nzeta"])
    za = np.concatenate([za + 2 * np.pi * j / nfp for j in range(nfp)])
    basis = {}
    # hmap_axisNB.F90:195-291 interpolates field-periodic hat components,
    # then rotates back; a Cartesian full-turn DFT aliases different modes.
    for name in ["xyz(::)", "Nxyz(::)", "Bxyz(::)"]:
        vectors = axis[name].values.T
        hat = vectors.copy()
        hat[:, 0] = np.cos(za) * vectors[:, 0] + np.sin(za) * vectors[:, 1]
        hat[:, 1] = -np.sin(za) * vectors[:, 0] + np.cos(za) * vectors[:, 1]
        hat = interpolate(hat[:nz], zeta, za[0], 2 * np.pi / nfp, 0)
        vectors = hat.copy()
        vectors[:, 0] = np.cos(zeta) * hat[:, 0] - np.sin(zeta) * hat[:, 1]
        vectors[:, 1] = np.sin(zeta) * hat[:, 0] + np.cos(zeta) * hat[:, 1]
        basis[name] = vectors

    def surface(name, sign):
        return interpolate(
            interpolate(
                boundary[name].values,
                sign * theta,
                float(boundary["theta(:)"][0]),
                2 * np.pi,
                0,
            ),
            sign * zeta,
            float(boundary["zeta(:)"][0]),
            2 * np.pi / nfp,
            1,
        )

    x, y = surface("X(::)", 1), surface("Y(::)", 1)
    removed = {"X": 0.0, "Y": 0.0}
    if int(boundary["lasym"]) == 0:
        # boundaryFromFile.F90:207-240: declared lasym=0 means Xcos/Ysin.
        x_negative, y_negative = surface("X(::)", -1), surface("Y(::)", -1)
        removed = {
            "X": float(abs((x - x_negative) / 2).max()),
            "Y": float(abs((y + y_negative) / 2).max()),
        }
        x, y = (x + x_negative) / 2, (y - y_negative) / 2
    return (
        basis["xyz(::)"][:, None, :]
        + x[:, :, None] * basis["Nxyz(::)"][:, None, :]
        + y[:, :, None] * basis["Bxyz(::)"][:, None, :]
    ), removed


def check_boundary(state, frame, tolerance):
    theta = 2 * np.pi * (np.arange(65) + 0.371) / 65
    zeta = 2 * np.pi * (np.arange(163) + 0.293) / (163 * state.nfp)
    target, removed = target_boundary(frame, theta, zeta)
    actual = state.evaluate("pos", rho=[1.0], theta=theta, zeta=zeta).pos
    actual = actual.isel(rad=0).transpose("tor", "pol", "xyz").values
    error = float(np.linalg.norm(actual - target, axis=-1).max())
    return {
        "source_frame_sha256": hashlib.sha256(Path(frame).read_bytes()).hexdigest(),
        "samples_theta_zeta": [65, 163],
        "source_declared_symmetry_projection_removed_m": removed,
        "boundary_position_max_error_m": error,
        "tolerance_m": tolerance,
        "boundary_goal_passed": bool(np.isfinite(error) and error <= tolerance),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["frame", "parameter", "state", "output"]:
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    import gvec

    goals = json.loads(Path(__file__).with_name("acceptance-goals.json").read_text())
    state = gvec.State(args.parameter.resolve(), args.state.resolve())
    report = check_boundary(state, args.frame, goals["boundary_position_max_error_m"])
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    raise SystemExit(0 if report["boundary_goal_passed"] else 2)
