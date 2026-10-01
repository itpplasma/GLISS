"""GLISS side of the QAS3 benchmark.

usage:
  compare.py replay FORT23             certified replay of TERPSICHORE's matrices
  compare.py mercier WOUT NAME         GLISS Mercier terms against VMEC's DMerc*
  compare.py inertia EXPORT MODEFILE   independent GLISS FEEC inertia count
  compare.py refine EXPORT MODEFILE DEGREE CELLS...
                                       lowest FEEC eigenvalue on each count of
                                       radial cells, independent of the
                                       equilibrium surfaces

GLISS_LIB selects the native library. Mercier and inertia use the export
written by convert_vmec (M = N = 8 by default; GLISS_MN overrides).
"""

import io
import os
import subprocess
import sys
import time
import warnings

import numpy as np

import gliss
import gliss.vmec


def replay(path):
    result = gliss.solve_terpsichore_fixed_boundary(path)
    print(
        f"GLISS replay: negative_count={result.negative_count} "
        f"eigenvalue={result.eigenvalue:.8e} "
        f"reference={result.reference_eigenvalue:.8e} "
        f"overlap={result.mode_overlap:.8f} certificate={result.certificate:.2e}"
    )


def mercier(wout, name):
    from scipy.io import netcdf_file

    order = int(os.environ.get("GLISS_MN", "8"))
    export = f"{name}_gliss.nc"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        gliss.vmec.convert_vmec(
            wout, export, poloidal_max=order, toroidal_max=order,
            overwrite=True, force_balance_policy="warn",
        )
    s, d_gliss = gliss.mercier_profile(export)
    variables = netcdf_file(wout, "r", mmap=False).variables
    s_vmec = np.linspace(0.0, 1.0, int(variables["ns"].data))
    d_vmec = np.array(variables["DMerc"].data)
    inside = (s_vmec > max(0.1, s[0])) & (s_vmec < min(0.9, s[-1]))
    interpolated = np.interp(s_vmec[inside], s, d_gliss)
    reference = d_vmec[inside]
    error = np.abs(interpolated - reference) / np.maximum(np.abs(reference), 1e-300)
    agree = np.mean(np.sign(interpolated) == np.sign(reference))
    print(
        f"Mercier 0.1<s<0.9: median |G-V|/|V|={np.median(error):.3e} "
        f"p90={np.quantile(error, 0.9):.3e} sign agreement={100 * agree:.1f}% "
        f"NaN surfaces={int(np.isnan(d_gliss).sum())}"
    )


def inertia(export, modefile):
    modes = [tuple(map(int, line.split())) for line in open(modefile) if line.strip()]
    with gliss.Equilibrium(export) as equilibrium:
        result = gliss.cas3d_marginality_inertia(
            equilibrium, modes=modes, parity_class=1, degree=1,
            angular_theta=96, angular_zeta=64,
        )
    print(f"GLISS independent: negative_count={result.negative_count}")


def refine(export, modefile, degree, *cells):
    modes = [tuple(map(int, line.split())) for line in open(modefile) if line.strip()]
    with gliss.Equilibrium(export) as equilibrium:
        for count in cells:
            started = time.perf_counter()
            result = gliss.solve_cas3d_marginality(
                equilibrium, modes=modes, parity_class=1, degree=int(degree),
                angular_theta=96, angular_zeta=64, radial_cells=int(count),
            )
            print(
                f"degree {degree} cells {result.radial_surfaces}: "
                f"negative_count={result.negative_count} "
                f"lowest={result.lowest_eigenvalue:.6e} "
                f"seconds={time.perf_counter() - started:.0f}",
                flush=True,
            )


if __name__ == "__main__":
    command, *arguments = sys.argv[1:]
    commands = {
        "replay": replay, "mercier": mercier, "inertia": inertia, "refine": refine,
    }
    commands[command](*arguments)
