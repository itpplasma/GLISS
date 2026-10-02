#!/usr/bin/env python3
"""Refine DCON wall distances after the public free-boundary runner finishes."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--walls", nargs="+", type=float,
                        default=[0.152422, 0.1525, 0.1535])
    args = parser.parse_args()
    args.output = args.output.resolve()
    binary = args.output / "dcon/GPEC/bin/dcon"
    template = args.output / "dcon_free/a0.15"
    rows, evidence = [], []
    for wall in args.walls:
        run = args.output / "dcon_refined" / f"a{wall:.9g}"
        run.mkdir(parents=True, exist_ok=True)
        for name in ("equil.in", "dcon.in", "sol.in", "vac.in"):
            shutil.copy2(template / name, run / name)
        path = run / "vac.in"
        text, replacements = re.subn(r"(?m)^(\s*a\s*=)\s*[^\n]+",
                                     rf"\g<1> {wall:.9g}", path.read_text())
        if replacements != 1:
            raise RuntimeError("expected exactly one conformal-wall parameter")
        path.write_text(text)
        with (run / "stdout.log").open("w") as log:
            subprocess.run([str(binary)], cwd=run, stdout=log,
                           stderr=subprocess.STDOUT, check=True)
        values = re.findall(r"Energies:.*real\s*=\s*([^,]+),",
                            (run / "stdout.log").read_text())
        if len(values) != 1:
            raise RuntimeError("expected exactly one DCON total energy")
        rows.append((wall, float(values[0])))
        evidence.append({"wall_half_widths": wall,
                         "input_hashes": {name: hashlib.sha256(
                             (run / name).read_bytes()).hexdigest() for name in
                             ("equil.in", "dcon.in", "sol.in", "vac.in")}})
    with (args.output / "dcon_walls_refined.csv").open("w") as result:
        writer = csv.writer(result)
        writer.writerow(("wall_half_widths", "total_energy"))
        writer.writerows(rows)
    (args.output / "dcon_refined_provenance.json").write_text(json.dumps({
        "binary_sha256": hashlib.sha256(binary.read_bytes()).hexdigest(),
        "runner_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "runs": evidence}, indent=2) + "\n")


if __name__ == "__main__":
    main()
