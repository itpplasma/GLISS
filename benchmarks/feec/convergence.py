#!/usr/bin/env python3
"""Fixed-export radial FE convergence, with source and binary provenance."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import time

import gliss


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("equilibrium", type=Path, nargs="+")
    parser.add_argument("--cells", type=int, nargs="+", default=[8, 16, 32, 64, 128])
    parser.add_argument("--degrees", type=int, nargs="+", default=[1, 2, 3, 4])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if any(cell < 2 for cell in args.cells):
        parser.error("cell counts must be at least two")
    if any(degree not in (1, 2, 3, 4) for degree in args.degrees):
        parser.error("degrees must be between one and four")
    source = Path(__file__).resolve().parents[2]
    library = Path(os.environ["GLISS_LIB"]).resolve()
    revision_file = source / "SOURCE_REVISION"
    revision = (
        revision_file.read_text().strip() if revision_file.exists()
        else subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=source, text=True
        ).strip()
    )
    report = {
        "revision": revision,
        "library": str(library),
        "library_sha256": digest(library),
        "space": "regular quotient, physical perpendicular inertia",
        "scope": "fixed export; empirical rates are not continuum certificates",
        "rows": [],
        "rates": [],
    }

    def save():
        args.output.write_text(json.dumps(report, indent=2) + "\n")

    save()
    for path in args.equilibrium:
        path = path.resolve()
        export_digest = digest(path)
        with gliss.Equilibrium(path) as equilibrium:
            for degree in args.degrees:
                rows = []
                for cells in args.cells:
                    start = time.monotonic()
                    result = gliss.solve_axisymmetric(
                        equilibrium, poloidal_max=6, degree=degree, radial_cells=cells
                    )
                    row = dict(
                        equilibrium=str(path), export_sha256=export_digest,
                        degree=degree, cells=cells,
                        negative_count=result.negative_count,
                        lowest_eigenvalue=result.lowest_eigenvalue,
                        certificate=result.certificate,
                        seconds=time.monotonic() - start,
                    )
                    rows.append(row)
                    report["rows"].append(row)
                    save()
                    print(json.dumps(row), flush=True)
                for a, b, c in zip(rows, rows[1:], rows[2:]):
                    if b["cells"] != 2*a["cells"] or c["cells"] != 2*b["cells"]:
                        continue
                    d1 = a["lowest_eigenvalue"] - b["lowest_eigenvalue"]
                    d2 = b["lowest_eigenvalue"] - c["lowest_eigenvalue"]
                    rate = math.log2(abs(d1/d2)) if d1 and d2 else None
                    report["rates"].append(dict(
                        equilibrium=str(path), degree=degree,
                        cells=[a["cells"], b["cells"], c["cells"]],
                        empirical_order=rate,
                        differences_resolved=(
                            abs(d1) > a["certificate"] + b["certificate"]
                            and abs(d2) > b["certificate"] + c["certificate"]
                        ),
                    ))
                save()


if __name__ == "__main__":
    main()
