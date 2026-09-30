"""Memory and quotient scaling of the sparse CAS3D2MN coefficient solve (#11).

usage: python run.py EXPORT OUTPUT.json [--entries 10 20 35 70]
       [--base M N] [--angular 64 64] [--surfaces STRIDE]

Solves the CAS3D2MN labeled phase envelope with the coefficient
normalization for envelope tables of increasing size, each in a fresh
process, and writes one manifest entry per table: labeled and physical mode
counts, quotient rank, labeled nullity, peak block width, inertia, lowest
eigenvalue, certificate, peak resident memory and wall time. The envelope
table is the first ENTRIES pairs (m, n) of the ordered list (0,0), then
increasing m + |n|, m, n with m >= 0, |n| <= 4 and n > 0 for m = 0 (envelope
n counts field periods; |n| <= 4 keeps the table resolved by 64 toroidal
points on a five-period equilibrium). The published CAS3D W7-X L139 table is not public; 70
entries (139 labels) reproduce its size, not its sidebands.
"""

import argparse
import json
import resource
import subprocess
import sys
import time


def envelope_table(entries):
    modes = [(0, 0)]
    candidates = [
        (m, n)
        for m in range(0, 16)
        for n in range(-4, 5)
        if (m, n) != (0, 0) and not (m == 0 and n < 0)
    ]
    candidates.sort(key=lambda mode: (mode[0] + abs(mode[1]), mode[0], mode[1]))
    modes.extend(candidates[: entries - 1])
    return modes


def solve(export, entries, base, angular):
    import gliss

    table = envelope_table(entries)
    start = time.perf_counter()
    with gliss.Equilibrium(export) as equilibrium:
        result = gliss.solve_cas3d_phase_envelope(
            equilibrium,
            base_mode=tuple(base),
            envelope_modes=table,
            degree=1,
            angular_theta=angular[0],
            angular_zeta=angular[1],
            normalization="cas3d2mn_coefficient",
            coefficient_angular_resolution=tuple(angular),
            reference_length=1.0,
        )
    elapsed = time.perf_counter() - start
    return {
        "envelope_entries": entries,
        "labeled_sidebands": result.labeled_sideband_count,
        "radial_surfaces": result.radial_surfaces,
        "quotient_rank": result.quotient_rank,
        "labeled_nullity": result.labeled_nullity,
        "peak_block_width": result.peak_block_width,
        "negative_count": result.negative_count,
        "lowest_eigenvalue": result.lowest_eigenvalue,
        "certificate": result.certificate,
        "peak_resident_mib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        / 1024.0,
        "wall_seconds": elapsed,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("export")
    parser.add_argument("output")
    parser.add_argument("--entries", type=int, nargs="+", default=[10, 20, 35, 70])
    parser.add_argument("--base", type=int, nargs=2, default=[1, 1])
    parser.add_argument("--angular", type=int, nargs=2, default=[64, 64])
    parser.add_argument("--single", action="store_true", help=argparse.SUPPRESS)
    arguments = parser.parse_args()
    if arguments.single:
        (entries,) = arguments.entries
        print(json.dumps(solve(arguments.export, entries, arguments.base,
                               arguments.angular)))
        return
    manifest = {"export": arguments.export, "base_mode": arguments.base,
                "angular_resolution": arguments.angular, "runs": []}
    for entries in arguments.entries:
        # A fresh process per table so that peak memory is per run.
        output = subprocess.run(
            [sys.executable, __file__, arguments.export, "-", "--single",
             "--entries", str(entries), "--base", *map(str, arguments.base),
             "--angular", *map(str, arguments.angular)],
            capture_output=True, text=True,
        )
        if output.returncode != 0:
            sys.exit(f"{entries} entries failed:\n{output.stderr}")
        run = json.loads(output.stdout.strip().splitlines()[-1])
        manifest["runs"].append(run)
        print(json.dumps(run), flush=True)
        with open(arguments.output, "w", encoding="utf-8") as handle:
            json.dump(manifest, handle, indent=2)


if __name__ == "__main__":
    main()
