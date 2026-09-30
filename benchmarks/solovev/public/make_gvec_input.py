#!/usr/bin/env python3
"""Write a GVEC fixed-boundary parameter file for the GPEC Solov'ev q0."""
import argparse
from pathlib import Path

import numpy as np
from scipy.interpolate import make_interp_spline

import solovev


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("q0", type=float)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--modes", type=int, default=20)
    parser.add_argument("--elements", type=int, default=32)
    parser.add_argument("--degree", type=int, default=5)
    parser.add_argument("--tolerance", type=float, default=1e-10)
    args = parser.parse_args()
    s, iota, pressure, psi_n, edge_flux = solovev.profiles(args.q0)
    nodes = np.linspace(0.0, 1.0, 101)
    iota_nodes = make_interp_spline(s, iota, k=3)(nodes)
    pressure_nodes = make_interp_spline(s, pressure, k=3)(nodes)
    r_cos, z_sin = solovev.boundary_harmonics(args.modes)

    def array(values):
        return "[" + ", ".join(repr(float(v)) for v in values) + "]"

    m = args.modes
    lines = [
        'ProjectName = "SOL"', "whichInitEquilibrium = 0", "init_LA = true",
        f"PHIEDGE = {float(edge_flux)!r}", "which_hmap = 1", "nfp = 1",
        f"X1X2_deg = {args.degree}", f"LA_deg = {args.degree}",
        f"degGP = {args.degree + 3}",
        f"X1_mn_max = [{m}, 0]", f"X2_mn_max = [{m}, 0]", f"LA_mn_max = [{m}, 0]",
        'X1_sin_cos = "_cos_"', 'X2_sin_cos = "_sin_"', 'LA_sin_cos = "_sin_"',
        "MinimizerType = 10", "PrecondType = 1", "start_dt = 0.5",
        "maxIter = 40000", f"minimize_tol = {args.tolerance}", "",
        "[sgrid]", "grid_type = 0", f"nElems = {args.elements}", "",
        "[iota]", 'type = "interpolation"', f"rho2 = {array(nodes)}",
        f"vals = {array(iota_nodes)}", "",
        "[pres]", 'type = "interpolation"', f"rho2 = {array(nodes)}",
        f"vals = {array(pressure_nodes)}", "",
        "[X1_b_cos]", *[f'"({i}, 0)" = {float(r_cos[i])!r}' for i in range(m + 1)],
        "", "[X2_b_sin]",
        *[f'"({i}, 0)" = {float(z_sin[i])!r}' for i in range(1, m + 1)],
        "", "[X1_a_cos]", f'"(0, 0)" = {solovev.R0!r}', "",
    ]
    args.directory.mkdir(parents=True, exist_ok=True)
    (args.directory / "parameter.toml").write_text("\n".join(lines))
    np.savetxt(args.directory / "profiles.txt",
               np.c_[s, iota, pressure, psi_n], header="s iota p_Pa psi_n")


if __name__ == "__main__":
    main()
