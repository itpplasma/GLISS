#!/usr/bin/env python3
"""Compare public Solov'ev stability signs with independent DCON counts."""
import argparse
import dataclasses
import hashlib
import json
import os
from pathlib import Path
import resource
import subprocess
import time

import gliss


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--cells", nargs="+", type=int, default=[32, 64, 128])
    parser.add_argument("--degree", type=int, default=3)
    parser.add_argument("--poloidal-max", type=int, default=6)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[2]
    library = Path(gliss._load_library()._name).resolve()
    report = {
        "commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=repo, text=True).strip(),
        "patch_sha256": hashlib.sha256(subprocess.check_output(
            ["git", "diff", "--binary", "HEAD"], cwd=repo)).hexdigest(),
        "runner_sha256": sha256(Path(__file__)),
        "library_sha256": sha256(library),
        "thread_counts": {key: os.environ.get(key) for key in
                          ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS")},
        "reference": "GPEC DCON f5595c0689c624834d720068ddc8b5a7e4028248; "
                     "q0=1.035 has one Newcomb crossing, q0=1.045 has zero; "
                     "see dcon/run.sh for the independent reference runner",
        "observable": "negative inertia; raw eigenvalues have different norms",
        "timing": "one observation, including equilibrium loading; shared host",
        "cases": [],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    for q0, expected in (("1.035", 1), ("1.045", 0)):
        source = repo / "test/data" / f"solovev_q{q0}.nc"
        for cells in args.cells:
            start = time.perf_counter()
            row = {"q0": q0, "radial_cells": cells,
                   "degree": args.degree, "input_sha256": sha256(source),
                   "expected_negative_count": expected}
            try:
                with gliss.Equilibrium(source) as equilibrium:
                    result = gliss.solve_axisymmetric(
                        equilibrium, poloidal_max=args.poloidal_max,
                        degree=args.degree, radial_cells=cells)
                row.update(dataclasses.asdict(result))
                row["accepted"] = result.negative_count == expected
            except (RuntimeError, ValueError) as error:
                row.update(error=str(error), accepted=False)
            row["seconds"] = time.perf_counter() - start
            report["cases"].append(row)
            args.output.write_text(json.dumps(report, indent=2, allow_nan=False)
                                   + "\n")
            print(q0, cells, row.get("negative_count"),
                  row.get("lowest_eigenvalue"), row["seconds"], flush=True)
    report["peak_rss_kib"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    if not all(row["accepted"] for row in report["cases"]):
        raise SystemExit("DCON sign comparison failed; inspect the saved cases")


if __name__ == "__main__":
    main()
