#!/usr/bin/env python3
"""Measure errors of the public R0=1, a=.33, e=1.6, F=1 Solov'ev family.

This reports generator input, stored half-grid data and an independent
SciPy estimate of GLISS's edge extrapolation. It does not admit an export
or certify physical matching of the plasma to the vacuum.
"""
import argparse
from decimal import Decimal, localcontext
import hashlib
import json
import math
from pathlib import Path
import tomllib

import netCDF4
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("parameter", type=Path)
parser.add_argument("export", type=Path)
parser.add_argument("--q0", type=float, default=1.5)
parser.add_argument("--theta-points", type=int, default=32768)
args = parser.parse_args()
if args.theta_points < 256:
    parser.error("use at least 256 angular points")
if not math.isfinite(args.q0) or args.q0 <= 0:
    parser.error("q0 must be positive and finite")

parameter_bytes = args.parameter.read_bytes()
export_bytes = args.export.read_bytes()
parameter = tomllib.loads(parameter_bytes.decode())
data = netCDF4.Dataset("qualification-memory", memory=export_bytes)
if np.asarray(data["n"][:]).tolist() != [0]:
    raise ValueError("this axisymmetric qualifier requires only n=0")
if int(np.asarray(data["N_FP"][:]).item()) != 1:
    raise ValueError("this public Solov'ev family requires N_FP=1")

major, minor, elongation, field = 1.0, 0.33, 1.6, 1.0
psi0 = elongation * field * minor**2 / (2 * args.q0 * major)
p0 = (2 * psi0**2 * (elongation**2 + 1)
      / (minor * major * elongation)**2 / (4e-7 * np.pi))


def decimal_reference(psi_n, q0):
    """Independent 80-digit positive hypergeometric series."""
    with localcontext() as context:
        context.prec = 80
        d = Decimal
        pi = d("3.141592653589793238462643383279502884197169399375105820"
               "974944592307816406286208998628")
        a, e, r0, f = d("0.33"), d("1.6"), d(1), d(1)
        u = d(str(psi_n))
        z = 4 * a**2 * u / r0**2

        def hyper(c):
            total = term = d(1)
            for n in range(1, 350):
                term *= ((d(n) - d("0.25")) * (d(n) + d("0.25"))
                         * z / (d(n) * (d(n) - 1 + c)))
                total += term
            return total

        return d(str(q0)) * hyper(d(1)), pi * f * e * a**2 * u / r0 * hyper(d(2))


edge_q_decimal, edge_flux_decimal = decimal_reference(1, args.q0)
edge_flux = float(edge_flux_decimal)
reference_t = 2 * np.pi * np.arange(2048) / 2048


def reference(psi_n):
    if psi_n == 0:
        return args.q0, 0.0
    r = np.sqrt(major**2 + 2 * minor * major * np.sqrt(psi_n)
                * np.cos(reference_t))
    inverse_cube = r**-3
    return (args.q0 * major**3 * np.mean(inverse_cube),
            2 * np.pi * field * elongation * minor**2 * major**2 * psi_n
            * np.mean(np.sin(reference_t)**2 * inverse_cube))


def psi_at_s(s):
    if s <= 0:
        return 0.0
    if s >= 1:
        return 1.0
    return brentq(lambda u: reference(u)[1] / edge_flux - s, 0, 1,
                  xtol=2e-15, rtol=1e-14)


def implicit(r, z):
    return (r**2 * z**2 / elongation**2
            + (r**2 - major**2)**2 / 4) / (minor * major)**2


def maximum(values):
    return float(np.max(np.abs(values)))


theta = np.arange(args.theta_points) / args.theta_points
t = 2 * np.pi * theta
modes = np.asarray(data["m"][:])
cosine, sine = np.cos(modes[:, None] * t), np.sin(modes[:, None] * t)
s = np.asarray(data["rho"][:])**2
if len(s) < 4 or not np.all(np.isfinite(s)) or not np.all(np.diff(s) > 0):
    raise ValueError("at least four ordered finite half-grid surfaces are required")
if s[0] <= 0 or s[-1] >= 1:
    raise ValueError("this diagnostic requires a strictly interior half-grid")
parity_power = np.where(modes == 0, 0, np.where(modes % 2 == 1, 0.5, 1))
factor = s[:, None]**parity_power


def reconstruct(name, edge=False):
    c = np.asarray(data[name + "_mnc"][:, :, 0])
    z = np.asarray(data[name + "_mns"][:, :, 0])
    if edge:
        c = CubicSpline(s, c / factor, axis=0, bc_type="not-a-knot")(1)
        z = CubicSpline(s, z / factor, axis=0, bc_type="not-a-knot")(1)
    return c @ cosine + z @ sine


x, y, z = [reconstruct(name) for name in ("xhat", "yhat", "zhat")]
r = np.hypot(x, y)
psi_chi = -np.asarray(data["chi"][:]) / (2 * np.pi * psi0)
psi_s = np.array([psi_at_s(value) for value in s])
geometry_psi = implicit(r, z)
residual = geometry_psi - psi_chi[:, None]
surface_index, angle_index = np.unravel_index(np.argmax(abs(residual)), residual.shape)

input_r = sum(value * np.cos(int(key.strip("()").split(",")[0]) * t)
              for key, value in parameter["X1_b_cos"].items())
input_z = sum(value * np.sin(int(key.strip("()").split(",")[0]) * t)
              for key, value in parameter["X2_b_sin"].items())
exact_r = np.sqrt(major**2 + 2 * minor * major * np.cos(t))
exact_z = elongation * minor * major * np.sin(t) / exact_r
input_s = np.asarray(parameter["iota"]["rho2"])
input_psi = np.array([psi_at_s(value) for value in input_s])
input_reference_iota = np.array([1 / reference(value)[0] for value in input_psi])
input_pressure_s = np.asarray(parameter["pres"]["rho2"])
input_pressure_psi = np.array([psi_at_s(value) for value in input_pressure_s])

xe, ye, ze = [reconstruct(name, edge=True) for name in ("xhat", "yhat", "zhat")]
re = np.hypot(xe, ye)
edge_geometry_psi = implicit(re, ze)
chi_spline = CubicSpline(s, np.asarray(data["chi"][:]), bc_type="not-a-knot")
phi_spline = CubicSpline(s, np.asarray(data["Phi"][:]), bc_type="not-a-knot")
iota_spline = CubicSpline(s, -np.asarray(data["iota"][:]), bc_type="not-a-knot")
edge_psi_chi = -chi_spline(1) / (2 * np.pi * psi0)
edge_iota_flux = -chi_spline(1, 1) / phi_spline(1, 1)
edge_shape_t = np.arctan2(re * ze / (elongation * minor * major),
                        (re**2 - major**2) / (2 * minor * major))
edge_exact_r = np.sqrt(major**2 + 2 * minor * major * np.cos(edge_shape_t))
edge_exact_z = elongation * minor * major * np.sin(edge_shape_t) / edge_exact_r

generator = {
    "R_boundary_max_abs_error_m": maximum(input_r - exact_r),
    "Z_boundary_max_abs_error_m": maximum(input_z - exact_z),
    "implicit_boundary_max_abs_error": maximum(implicit(input_r, input_z) - 1),
    "iota_nodes_max_abs_error": maximum(np.asarray(parameter["iota"]["vals"])
                                        - input_reference_iota),
    "pressure_nodes_max_abs_error_Pa": maximum(np.asarray(parameter["pres"]["vals"])
                                               - p0 * (1 - input_pressure_psi)),
    "phi_edge_abs_error_Wb": float(abs(Decimal(str(parameter["PHIEDGE"]))
                                      - edge_flux_decimal)),
}
halfgrid = {
    "s_min": float(s[0]), "s_max": float(s[-1]),
    "phi_over_s_max_abs_error_Wb": maximum(np.asarray(data["Phi"][:]) / s - edge_flux),
    "psi_chi_vs_reference_s_max_abs_error": maximum(psi_chi - psi_s),
    "implicit_vs_chi_max_abs_error": maximum(residual),
    "implicit_vs_reference_s_max_abs_error": maximum(geometry_psi - psi_s[:, None]),
    "iota_vs_reference_s_max_abs_error": maximum(-np.asarray(data["iota"][:])
                        - np.array([1 / reference(value)[0] for value in psi_s])),
    "pressure_vs_reference_s_max_abs_error_Pa": maximum(np.asarray(data["p"][:])
                                                      - p0 * (1 - psi_s)),
    "axis_near_max_abs_implicit_error": maximum(residual[0]),
    "central_max_abs_implicit_error": maximum(residual[len(s) // 2]),
    "maximum_location": {
        "surface_one_based": int(surface_index + 1),
        "s": float(s[surface_index]), "theta_periods": float(theta[angle_index]),
        "R_m": float(r[surface_index, angle_index]),
        "Z_m": float(z[surface_index, angle_index]),
        "analytic_psi_from_chi": float(psi_chi[surface_index]),
        "implicit_psi_from_geometry": float(geometry_psi[surface_index, angle_index]),
        "signed_error": float(residual[surface_index, angle_index]),
    },
}
edge = {
    "method": "source-equivalent parity-quotient not-a-knot extrapolation to s=1",
    "implicit_vs_boundary_1_max_abs_error": maximum(edge_geometry_psi - 1),
    "implicit_vs_chi_max_abs_error": maximum(edge_geometry_psi - edge_psi_chi),
    "psi_from_chi": float(edge_psi_chi),
    "phi_slope_Wb": float(phi_spline(1, 1)),
    "phi_slope_abs_error_Wb": float(abs(Decimal(str(float(phi_spline(1, 1))))
                                       - edge_flux_decimal)),
    "exported_iota_spline_positive": float(iota_spline(1)),
    "exported_iota_spline_abs_error": float(abs(iota_spline(1)
                                               - float(1 / edge_q_decimal))),
    "flux_ratio_iota_positive": float(edge_iota_flux),
    "flux_ratio_iota_abs_error": float(abs(edge_iota_flux - float(1 / edge_q_decimal))),
    "flux_ratio_minus_exported_iota": float(edge_iota_flux - iota_spline(1)),
    "shape_parameter_matched_R_max_abs_error_m": maximum(re - edge_exact_r),
    "shape_parameter_matched_Z_max_abs_error_m": maximum(ze - edge_exact_z),
    "shape_parameter_matched_distance_max_m": maximum(np.hypot(re - edge_exact_r,
                                                               ze - edge_exact_z)),
}
summary = {
    "export_sha256": hashlib.sha256(export_bytes).hexdigest(),
    "parameter_sha256": hashlib.sha256(parameter_bytes).hexdigest(),
    "qualifier_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "q0": args.q0, "theta_points": args.theta_points,
    "explicit_boundary_schema_version": getattr(data, "boundary_schema_version", None),
    "reference": {
        "phi_edge_Wb_decimal80": str(edge_flux_decimal),
        "q_edge_decimal80": str(edge_q_decimal), "psi0": psi0, "p0_Pa": p0,
        "half_width_m": float((np.sqrt(major**2 + 2 * minor * major)
                               - np.sqrt(major**2 - 2 * minor * major)) / 2),
    },
    "generator_input": generator,
    "export_halfgrid": halfgrid,
    "reconstructed_edge_estimate": edge,
    "domain_notes": [
        "No acceptance decision or matched exterior-field qualification is made.",
        "Export contains half-grid surfaces; edge metrics use native-equivalent extrapolation.",
        "Cylindrical R=hypot(xhat,yhat); positive analytic iota is minus flipped export iota.",
        "Parameter-matched shape distances are not nearest-point distances.",
        "Exported iota and reconstructed chi'/Phi' are separate quantities.",
    ],
}
data.close()
print(json.dumps(summary, indent=2))
