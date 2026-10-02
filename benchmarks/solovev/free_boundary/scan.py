#!/usr/bin/env python3
"""Free-boundary n=1 stability of a Solov'ev export against a conformal wall.

Every wall distance is given in units of the plasma half-width
``(R_max - R_min) / 2``, the convention of the VACUUM code's conformal
``ishape=6`` wall used by DCON. One CSV row per wall distance reports the
negative eigenvalue count and lowest omega^2 of parity class 1. With
``--bisect LOW HIGH`` the walls are chosen by bisection of the critical
distance between a stable LOW and an unstable HIGH to ``--tolerance``, and a
final row ``critical,LOW,HIGH`` reports the bracket.
"""
import argparse
import csv
import sys
from pathlib import Path

import gliss


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("export", type=Path)
    parser.add_argument("walls", type=float, nargs="*")
    parser.add_argument("--half-width", type=float, default=0.35245,
                        help="plasma half-width in metres")
    parser.add_argument("--poloidal-max", type=int, default=8)
    parser.add_argument("--degree", type=int, default=2)
    parser.add_argument("--angular", type=int, nargs=2, default=(128, 8))
    parser.add_argument("--edge", type=int, nargs=2, default=(128, 128))
    parser.add_argument("--bisect", type=float, nargs=2, metavar=("LOW", "HIGH"))
    parser.add_argument("--tolerance", type=float, default=1e-3)
    args = parser.parse_args()
    modes = [(0, 1)] + [
        (m, n) for m in range(1, args.poloidal_max + 1) for n in (-1, 1)
    ]
    writer = csv.writer(sys.stdout)
    writer.writerow(["wall_half_widths", "wall_distance_m", "negative_count",
                     "lowest_eigenvalue", "certificate"])
    with gliss.Equilibrium(args.export) as equilibrium:

        def unstable(wall):
            distance = wall * args.half_width
            vacuum = gliss.VacuumModel(tuple(args.edge), distance)
            with gliss.StabilityProblem(
                # A floor far below the kink eigenvalues near the threshold.
                equilibrium, modes, degree=args.degree, zero_floor=1e-8,
                angular_theta=args.angular[0], angular_zeta=args.angular[1],
                vacuum=vacuum,
            ) as problem:
                result = problem.solve_class(1)
            writer.writerow([wall, distance, result.negative_count,
                             result.lowest_eigenvalue, result.certificate])
            sys.stdout.flush()
            return result.negative_count > 0

        for wall in args.walls:
            unstable(wall)
        if args.bisect is None:
            return
        low, high = args.bisect
        if unstable(low) or not unstable(high):
            sys.exit("--bisect needs a stable LOW and an unstable HIGH wall")
        while high - low > args.tolerance:
            middle = 0.5 * (low + high)
            if unstable(middle):
                high = middle
            else:
                low = middle
        writer.writerow(["critical", low, high, "", ""])


if __name__ == "__main__":
    main()
