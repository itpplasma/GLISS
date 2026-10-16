"""Compare pressure reverse assembly with public cardinal JVPs on one basis."""

import argparse
import ctypes
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import statistics
import subprocess
import sys
import time

import numpy as np

import gliss
from gliss.pressure_derivatives import _bind, _pointer


REPOSITORY = Path(__file__).resolve().parents[2]
THREAD_VARIABLES = (
    "OMP_NUM_THREADS", "OMP_DYNAMIC", "OMP_PROC_BIND", "OMP_PLACES",
    "OPENBLAS_NUM_THREADS", "GOTO_NUM_THREADS", "MKL_NUM_THREADS",
    "MKL_DYNAMIC", "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
)


def command_output(command):
    try:
        result = subprocess.run(command, capture_output=True, text=True, timeout=10)
    except (OSError, subprocess.TimeoutExpired):
        return None
    return result.stdout.strip() if result.returncode == 0 else None


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def toolchain(library, build_dir):
    cache_path = build_dir / "CMakeCache.txt"
    cache = {}
    if cache_path.is_file():
        for line in cache_path.read_text().splitlines():
            if line.startswith(("#", "//")) or "=" not in line or ":" not in line:
                continue
            key, value = line.split("=", 1)
            cache[key.split(":", 1)[0]] = value
    keys = (
        "CMAKE_BUILD_TYPE", "CMAKE_Fortran_COMPILER", "CMAKE_Fortran_FLAGS",
        "CMAKE_Fortran_FLAGS_RELEASE", "CMAKE_Fortran_FLAGS_RELWITHDEBINFO",
        "BLAS_LIBRARIES", "BLAS_blas_LIBRARY", "BLAS_openblas_LIBRARY",
        "OpenMP_Fortran_FLAGS", "CMAKE_EXE_LINKER_FLAGS", "CMAKE_SHARED_LINKER_FLAGS",
        "GLISS_OPENMP_RUNTIME", "GLISS_OPENMP_SELECTED_LIBRARY",
        "OpenMP_Fortran_LIB_NAMES", "OpenMP_omp_LIBRARY",
    )
    compiler = cache.get("CMAKE_Fortran_COMPILER")
    inspection = ["otool", "-L", str(library)] if sys.platform == "darwin" else [
        "ldd", str(library),
    ]
    dependencies = command_output(inspection)
    dependency_hashes = {}
    if dependencies:
        for line in dependencies.splitlines():
            candidate = line.strip().split(" => ")[-1].split(" (", 1)[0]
            path = Path(candidate)
            if path.is_file():
                dependency_hashes[str(path.resolve())] = sha256(path)
    commands_path = build_dir / "compile_commands.json"
    compile_commands = []
    if commands_path.is_file():
        affected = {"compatible_compressible_stiffness_assembly.f90",
                    "compatible_three_component_problem.f90", "radial_cubic_spline.f90"}
        compile_commands = [entry for entry in json.loads(commands_path.read_text())
                            if Path(entry["file"]).name in affected]
    return {
        "cmake_cache": str(cache_path) if cache_path.is_file() else None,
        "cmake_cache_sha256": sha256(cache_path) if cache_path.is_file() else None,
        "cache_values": {key: cache[key] for key in keys if key in cache},
        "fortran_compiler_version": command_output([compiler, "--version"]) if compiler else None,
        "linked_dependencies": dependencies,
        "linked_dependency_sha256": dependency_hashes,
        "compile_commands_sha256": sha256(commands_path) if commands_path.is_file() else None,
        "affected_compile_commands": compile_commands,
        "metadata_limit": None if cache else "Build cache unavailable; compiler flags unknown.",
    }


def blas_runtime(library):
    """Query the linked OpenBLAS without changing its process-wide settings."""
    result = {}
    for name in ("openblas_get_parallel", "openblas_get_num_threads", "openblas_get_config"):
        try:
            function = getattr(library, name)
        except AttributeError:
            continue
        function.argtypes = []
        function.restype = ctypes.c_char_p if name.endswith("config") else ctypes.c_int
        value = function()
        result[name] = value.decode("ascii") if isinstance(value, bytes) else value
    result["scope"] = "Read outside a parallel region; actual BLAS team size can differ inside assembly."
    return result


def check_status(status, error):
    if status:
        raise RuntimeError(f"GLISS status {status}: {error.value.decode('utf-8')}")


def actions(problem, parity, basis, sample_count):
    """Call the public C ABI with a fixed, already certified cluster basis."""
    _bind(problem._library)
    error = ctypes.create_string_buffer(512)
    gradient = np.empty(sample_count)
    direction = np.zeros(sample_count)
    output = ctypes.c_double()

    def reverse():
        status = problem._library.gliss_stability_problem_pressure_trace_vjp(
            problem._handle, parity, basis.shape[1], basis.shape[0], _pointer(basis),
            1.0, sample_count, _pointer(gradient), error, len(error),
        )
        check_status(status, error)
        return gradient.copy()

    def cardinal():
        result = np.empty(sample_count)
        for sample in range(sample_count):
            direction[sample] = 1.0
            status = problem._library.gliss_stability_problem_pressure_trace_jvp(
                problem._handle, parity, basis.shape[1], basis.shape[0], _pointer(basis),
                sample_count, _pointer(direction), ctypes.byref(output), error, len(error),
            )
            check_status(status, error)
            result[sample] = output.value
            direction[sample] = 0.0
        return result

    return reverse, cardinal


def timed(action):
    start = time.perf_counter()
    result = action()
    return result, time.perf_counter() - start


def distribution(samples):
    return {
        "samples_seconds": samples, "median_seconds": statistics.median(samples),
        "min_seconds": min(samples), "max_seconds": max(samples),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fixture", nargs="?", type=Path,
                        default=REPOSITORY / "test/data/solovev_q1.045.nc")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--build-dir", type=Path)
    parser.add_argument("--mode", type=int, nargs=2, action="append")
    parser.add_argument("--parity", type=int, choices=(0, 1, 2))
    parser.add_argument("--degree", type=int, default=1)
    parser.add_argument("--radial-cells", type=int, default=2)
    parser.add_argument("--angular", type=int, nargs=2, default=(64, 16))
    parser.add_argument("--cluster-size", type=int, default=1)
    parser.add_argument("--gap", type=float, default=1e-4)
    parser.add_argument("--warmup", type=int, default=1)
    parser.add_argument("--repeat", type=int, default=3)
    args = parser.parse_args()
    if args.warmup < 1 or args.repeat < 1 or args.cluster_size < 1:
        parser.error("warmup, repeat and cluster-size must be positive")
    stages = {}
    equilibrium, stages["read_equilibrium"] = timed(lambda: gliss.Equilibrium(args.fixture))
    with equilibrium:
        problem, stages["assemble_problem"] = timed(lambda: gliss.StabilityProblem(
            equilibrium, args.mode or [(1, 1)], degree=args.degree,
            radial_cells=args.radial_cells, angular_theta=args.angular[0],
            angular_zeta=args.angular[1],
        ))
        with problem:
            library = Path(problem._library._name).resolve()
            if not library.is_file():
                raise RuntimeError("Set GLISS_LIB to an actual library file to record its hash")
            library_hash = sha256(library)
            parity = args.parity if args.parity is not None else problem.parity_classes[0]
            spectrum, stages["solve_certified_spectrum"] = timed(
                lambda: problem.solve_full_spectrum_class(parity)
            )
            stop = spectrum.eigenvalues.size
            start = stop - args.cluster_size
            # Public selection certifies the exterior gap and full pressure
            # tangent domain. Timed C actions reuse precisely this basis.
            sensitivity, stages["admit_cluster_and_initial_sensitivity"] = timed(
                lambda: problem.spectral_pressure_sensitivity(parity, start, stop, gap=args.gap)
            )
            basis = np.ascontiguousarray(spectrum.eigenvectors[start:stop])
            nodes, pressure = problem.pressure_samples()
            if not np.all(np.isfinite(basis)) or not np.all(pressure > 0):
                raise RuntimeError("Benchmark requires finite geometry and positive pressure admission")
            reverse, cardinal = actions(problem, parity, basis, pressure.size)
            warmup_start = time.perf_counter()
            for _ in range(args.warmup):
                np.testing.assert_allclose(reverse(), cardinal(), rtol=2e-9, atol=1e-11)
            stages["warmup"] = time.perf_counter() - warmup_start
            times = {"reverse": [], "cardinal_jvps": []}
            agreement = []
            for repeat in range(args.repeat):
                results = {}
                order = ("reverse", "cardinal_jvps") if repeat % 2 == 0 else (
                    "cardinal_jvps", "reverse",
                )
                for name in order:
                    results[name], elapsed = timed(reverse if name == "reverse" else cardinal)
                    times[name].append(elapsed)
                np.testing.assert_allclose(results["reverse"], results["cardinal_jvps"],
                                           rtol=2e-9, atol=1e-11)
                np.testing.assert_allclose(results["reverse"], sensitivity.gradient,
                                           rtol=2e-9, atol=1e-11)
                difference = results["reverse"] - results["cardinal_jvps"]
                scale = float(np.linalg.norm(results["cardinal_jvps"]))
                agreement.append({"max_absolute": float(np.max(np.abs(difference))),
                                  "relative_l2": float(np.linalg.norm(difference) / max(scale, 1e-300))})
            if sha256(library) != library_hash:
                raise RuntimeError("Native library changed during the benchmark")
            manifest = {
                "schema_version": 1, "benchmark": "fixed_geometry_pressure_vjp",
                "created_utc": datetime.now(timezone.utc).isoformat(),
                "command": sys.argv,
                "scope": "Same admitted fixed-boundary imported operator and cluster basis; no force-balance derivative.",
                "fixture": {"path": str(args.fixture.resolve()), "sha256": sha256(args.fixture)},
                "library": {"path": str(library), "sha256": library_hash, "version": gliss.version()},
                "source": {"commit": command_output(["git", "-C", str(REPOSITORY), "rev-parse", "HEAD"]),
                           "status": command_output(["git", "-C", str(REPOSITORY), "status", "--short"]),
                           "runner_sha256": sha256(__file__)},
                "host": {"platform": platform.platform(), "machine": platform.machine(),
                         "processor": platform.processor(), "logical_cpus": os.cpu_count(),
                         "cpu_identity": command_output(["sysctl", "-n", "machdep.cpu.brand_string"])
                             if sys.platform == "darwin" else command_output(["lscpu"])},
                "toolchain": toolchain(library, args.build_dir or library.parent),
                "thread_environment": {name: os.environ.get(name) for name in THREAD_VARIABLES},
                "blas_runtime": blas_runtime(problem._library),
                "python": {"version": platform.python_version(), "numpy": np.__version__},
                "configuration": {"modes": args.mode or [(1, 1)], "parity": parity,
                                  "degree": args.degree, "radial_cells": args.radial_cells,
                                  "angular": args.angular, "sample_count": pressure.size,
                                  "unknowns": basis.shape[1], "cluster": [start, stop],
                                  "exterior_gap": sensitivity.gap, "warmup": args.warmup,
                                  "repeats": args.repeat},
                "admission": {"minimum_sample_pressure_pa": float(np.min(pressure)),
                              "sample_coordinates": nodes.tolist(), "spectral_value": sensitivity.value,
                              "max_selected_residual": float(np.max(spectrum.residuals[start:stop])),
                              "max_selected_resolution": float(np.max(spectrum.resolutions[start:stop]))},
                "stage_seconds": stages,
                "timings": {name: distribution(samples) for name, samples in times.items()},
                "cardinal_over_reverse_median_ratio": statistics.median(times["cardinal_jvps"])
                    / statistics.median(times["reverse"]),
                "agreement": {"rtol": 2e-9, "atol": 1e-11, "runs": agreement},
                "gradient": results["reverse"].tolist(),
                "operation_counts_per_run": {"reverse_quadrature_sweeps": 1,
                                             "cardinal_tangent_assemblies": pressure.size},
            }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"manifest": str(args.output), "timings": manifest["timings"],
                      "cardinal_over_reverse_median_ratio": manifest["cardinal_over_reverse_median_ratio"]}))


if __name__ == "__main__":
    main()
