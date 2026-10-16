#!/usr/bin/env python3
"""Compare the Cartesian geometry value caller with its published baseline."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shutil
import statistics
import subprocess
import tempfile
import time


BASELINE = "f2da7f4301e9207db38d3fd5fc3c1488b344fb56"
FLAGS = ["-O3", "-ffree-line-length-none", "-Igenerated"]


def digest(data):
    return hashlib.sha256(data).hexdigest()


def host_cpu():
    if platform.system() == "Darwin":
        return subprocess.check_output(
            ["sysctl", "-n", "machdep.cpu.brand_string"], text=True
        ).strip()
    cpuinfo = Path("/proc/cpuinfo")
    if cpuinfo.exists():
        for line in cpuinfo.read_text().splitlines():
            if line.startswith("model name"):
                return line.split(":", 1)[1].strip()
    return platform.processor()


def parse_timing(output):
    rows = [line.split() for line in output.splitlines()
            if re.match(r"^\s+\d+\s+", line)]
    if len(rows) != 11:
        raise RuntimeError("expected eleven completed timing pairs")
    original = [float(row[1]) for row in rows]
    current = [float(row[2]) for row in rows]
    if min(original + current) < 0.1:
        raise RuntimeError("timing blocks shorter than 100 ms; increase harness passes")
    checksums = [[float(row[3]), float(row[4])] for row in rows]
    scalar = {}
    for line in output.splitlines():
        if ":" in line:
            key, value = line.split(":", 1)
            scalar[key] = float(value)
    calls = int(scalar["calls per timed block"])
    return {
        "original_seconds": {"minimum": min(original),
                             "median": statistics.median(original),
                             "maximum": max(original)},
        "current_seconds": {"minimum": min(current),
                            "median": statistics.median(current),
                            "maximum": max(current)},
        "ratio_current_to_original": statistics.median(current) / statistics.median(original),
        "original_ns_per_call": statistics.median(original) * 1.e9 / calls,
        "current_ns_per_call": statistics.median(current) * 1.e9 / calls,
        "maximum_normalized_output_difference": scalar["maximum normalized output difference"],
        "minimum_normalized_frame_determinant": scalar["minimum normalized frame determinant"],
        "negative_orientation_samples": int(scalar["negative orientation samples"]),
        "calls_per_timed_block": calls,
        "trials": len(rows),
        "original_times_seconds": original,
        "current_times_seconds": current,
        "checksums": checksums,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compiler", default=os.environ.get("FC", "gfortran"))
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--output", type=Path, help="write JSON here; otherwise print it")
    args = parser.parse_args()
    compiler = shutil.which(args.compiler)
    if compiler is None:
        parser.error(f"compiler not found: {args.compiler}")
    version = subprocess.check_output([compiler, "--version"], text=True).splitlines()[0]
    if "GNU Fortran" not in version:
        parser.error("this recipe uses GNU Fortran flags")
    repo = args.repo.resolve()
    work = Path(tempfile.mkdtemp(prefix="gliss-geometry-"))
    (work / "generated").mkdir()
    original = subprocess.check_output(
        ["git", "show", f"{BASELINE}:src/cartesian_primitive_geometry.f90"], cwd=repo
    )
    old_name = b"module cartesian_primitive_geometry"
    if original.count(old_name) != 2:
        raise RuntimeError("unexpected baseline module declaration")
    renamed = original.replace(old_name, b"module cartesian_primitive_geometry_original")
    files = {
        "cartesian_primitive_geometry_original.f90": renamed,
        "cartesian_primitive_geometry.f90": (repo / "src/cartesian_primitive_geometry.f90").read_bytes(),
        "generated/gliss_geometry_generated.f90": (repo / "src/generated/gliss_geometry_generated.f90").read_bytes(),
        "generated/gliss_geometry_primal.inc": (repo / "src/generated/gliss_geometry_primal.inc").read_bytes(),
        "benchmark.f90": Path(__file__).with_name("benchmark.f90").read_bytes(),
    }
    for name, content in files.items():
        (work / name).write_bytes(content)
    commands = []

    def run(command):
        commands.append(command)
        result = subprocess.run(command, cwd=work, text=True, capture_output=True)
        with (work / "build.log").open("a") as log:
            log.write(result.stdout + result.stderr)
        result.check_returncode()
        return result

    sources = ["cartesian_primitive_geometry_original.f90",
               "generated/gliss_geometry_generated.f90",
               "cartesian_primitive_geometry.f90", "benchmark.f90"]
    for source in sources:
        run([compiler, *FLAGS, "-c", source])
    run([compiler, *FLAGS, *(Path(source).stem + ".o" for source in sources),
         "-o", "geometry_benchmark"])
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    started = time.time()
    timed = subprocess.run([str(work / "geometry_benchmark")], cwd=work, env=env,
                           text=True, capture_output=True)
    (work / "timing.log").write_text(timed.stdout + timed.stderr)
    timed.check_returncode()
    result = {
        "schema_version": 1,
        "baseline_revision": BASELINE,
        "baseline_git_file_sha256": digest(original),
        "compiler": version,
        "compiler_path": compiler,
        "flags": FLAGS,
        "host": {"cpu": host_cpu(), "architecture": platform.machine(),
                 "system": platform.system(), "release": platform.release()},
        "source_sha256": {name: digest(content) for name, content in files.items()},
        "commands": commands,
        "OMP_NUM_THREADS": "1",
        "started_unix": started,
        "elapsed_seconds": time.time() - started,
        "scratch_directory": str(work),
        "scope": "Serial regular-point value kernel; not application throughput or a universal performance claim.",
        "measurement": parse_timing(timed.stdout),
    }
    serialized = json.dumps(result, indent=2) + "\n"
    (work / "results.json").write_text(serialized)
    if args.output:
        args.output.write_text(serialized)
    else:
        print(serialized, end="")


if __name__ == "__main__":
    main()
