"""CAS3D2MN coefficient-envelope experiments with frozen conventions.

The default synthetic envelope measures sparse scaling, not agreement with
a published CAS3D deck. An ordered CSV enables explicitly declared literature
prefixes; unspecified historical settings remain hypotheses in the manifest.
"""

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import resource
import subprocess
import sys
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def envelope_table(entries):
    """Retain the original synthetic order, with an explicit finite limit."""
    candidates = [
        (m, n)
        for m in range(16)
        for n in range(-4, 5)
        if (m, n) != (0, 0) and not (m == 0 and n < 0)
    ]
    candidates.sort(key=lambda mode: (mode[0] + abs(mode[1]), mode[0], mode[1]))
    if entries < 1 or entries > len(candidates) + 1:
        raise ValueError(f"synthetic entries must be between 1 and {len(candidates) + 1}")
    return [(0, 0), *candidates[: entries - 1]]


def read_envelope_csv(path, n_sign):
    """Read m,n in file order; extra provenance columns do not alter order."""
    if n_sign not in ("gliss", "opposite"):
        raise ValueError("CSV input requires an explicit envelope n sign")
    with open(path, newline="", encoding="utf-8-sig") as handle:
        rows = csv.DictReader(line for line in handle if not line.lstrip().startswith("#"))
        if rows.fieldnames is None or not {"m", "n"}.issubset(rows.fieldnames):
            raise ValueError("envelope CSV must have m,n header columns")
        if len(set(rows.fieldnames)) != len(rows.fieldnames):
            raise ValueError("envelope CSV must not repeat header columns")
        table = []
        for index, row in enumerate(rows, 1):
            try:
                mode = (int(row["m"]), int(row["n"]))
            except (ValueError, TypeError) as error:
                raise ValueError(f"CSV envelope row {index} needs integer m,n") from error
            if mode[0] < 0 or any(abs(value) > 2**31 - 1 for value in mode):
                raise ValueError(f"CSV envelope row {index} has an invalid mode")
            if n_sign == "opposite":
                mode = (mode[0], -mode[1])
            if mode in table:
                raise ValueError(f"CSV envelope row {index} duplicates {mode}")
            table.append(mode)
    if not table or table[0] != (0, 0):
        raise ValueError("ordered envelope CSV must start with (0,0)")
    return table


def peak_resident_mib(raw_rss, system=None):
    """Darwin reports ru_maxrss in bytes; Linux reports KiB."""
    system = platform.system() if system is None else system
    if system == "Darwin":
        return raw_rss / (1024.0 * 1024.0)
    if system == "Linux":
        return raw_rss / 1024.0
    raise ValueError(f"peak RSS units are not declared for {system}")


def positive_integer(value):
    result = int(value)
    if result < 1:
        raise argparse.ArgumentTypeError("must be positive")
    return result


def positive_real(value):
    result = float(value)
    if not math.isfinite(result) or result <= 0.0:
        raise argparse.ArgumentTypeError("must be finite and positive")
    return result


def argument_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("export")
    parser.add_argument("output")
    prefixes = parser.add_mutually_exclusive_group()
    prefixes.add_argument("--entries", type=positive_integer, nargs="+")
    prefixes.add_argument("--literature-l", type=positive_integer, nargs="+",
                          help="published L values; requires CSV and prefix convention")
    parser.add_argument("--envelope-csv", type=Path,
                        help="ordered m,n CSV, starting with (0,0)")
    parser.add_argument("--envelope-n-sign", choices=("gliss", "opposite"),
                        help="CSV n already follows GLISS, or must be negated")
    parser.add_argument("--prefix-convention", choices=("caption_l", "terminal_l_plus_one"),
                        help="first L CSV rows, or first L+1 rows")
    parser.add_argument("--base", type=int, nargs=2, default=[1, 1])
    parser.add_argument("--parity", type=int, choices=(1, 2), default=1,
                        help="single parity class; coefficient normalization supports 1 or 2")
    parser.add_argument("--angular", type=positive_integer, nargs=2, default=[64, 64],
                        help="stiffness assembly theta,zeta grid")
    parser.add_argument("--coefficient-angular", type=positive_integer, nargs=2,
                        help="coefficient normalization grid; default: assembly grid")
    parser.add_argument("--reference-length", type=positive_real, default=1.0)
    parser.add_argument("--radial-cells", type=positive_integer,
                        help="native FE radial cells; default: equilibrium grid")
    parser.add_argument("--radial-quadrature", choices=("gauss5", "cas3d_midpoint"),
                        default="gauss5")
    parser.add_argument("--assumption", action="append", default=[],
                        help="record an explicit historical hypothesis; repeatable")
    parser.add_argument("--single", action="store_true", help=argparse.SUPPRESS)
    return parser


def configure(arguments, parser):
    if arguments.radial_cells is not None and arguments.radial_cells < 2:
        parser.error("--radial-cells must be at least 2")
    if arguments.envelope_csv is None:
        if arguments.envelope_n_sign or arguments.literature_l or arguments.prefix_convention:
            parser.error("CSV sign and literature prefixes require --envelope-csv")
        runs = [{"envelope_entries": count} for count in
                (arguments.entries or [3, 10, 20, 35, 70])]
        table = envelope_table(max(run["envelope_entries"] for run in runs))
        scope = "synthetic_sparse_scaling_control"
    else:
        if arguments.envelope_n_sign is None:
            parser.error("--envelope-csv requires --envelope-n-sign")
        table = read_envelope_csv(arguments.envelope_csv, arguments.envelope_n_sign)
        if arguments.literature_l:
            if arguments.prefix_convention is None:
                parser.error("--literature-l requires --prefix-convention")
            offset = int(arguments.prefix_convention == "terminal_l_plus_one")
            runs = [{"literature_l": label, "envelope_entries": label + offset}
                    for label in arguments.literature_l]
        else:
            if arguments.prefix_convention:
                parser.error("--prefix-convention requires --literature-l")
            runs = [{"envelope_entries": count} for count in
                    (arguments.entries or [len(table)])]
        if max(run["envelope_entries"] for run in runs) > len(table):
            parser.error("requested prefix exceeds the ordered CSV table")
        scope = "ordered_literature_envelope_hypotheses"
    if arguments.single and len(runs) != 1:
        parser.error("--single requires one prefix")
    configuration = {
        "base_mode": arguments.base,
        "parity_class": arguments.parity,
        "degree": 1,
        "normalization": "cas3d2mn_coefficient",
        "assembly_angular_resolution": arguments.angular,
        "coefficient_angular_resolution": arguments.coefficient_angular or arguments.angular,
        "reference_length": arguments.reference_length,
        "radial_cells": arguments.radial_cells,
        "radial_quadrature": arguments.radial_quadrature,
        "envelope_n_sign": arguments.envelope_n_sign or "gliss",
        "prefix_convention": arguments.prefix_convention or "explicit_entry_count",
        "scope": scope,
        "assumptions": arguments.assumption,
        "boundary_condition": "fixed",
        "base_fourier_convention": "2*pi*(M*theta-N*zeta/N_T)",
        "envelope_fourier_convention": "2*pi*(m*theta-n*zeta)",
        "interpretation": "coefficient Ritz value; not an SI frequency or growth rate",
    }
    return configuration, table, runs


def native_library_identity(library):
    path = Path(library._name).resolve()
    if not path.is_file():
        raise ValueError("set GLISS_LIB to a file path so the loaded library can be hashed")
    return {"path": str(path), "sha256": sha256_file(path)}


def solve(export, table, configuration):
    import gliss

    start = time.perf_counter()
    with gliss.Equilibrium(export) as equilibrium:
        library = native_library_identity(equilibrium._library)
        result = gliss.solve_cas3d_phase_envelope(
            equilibrium,
            base_mode=tuple(configuration["base_mode"]),
            envelope_modes=table,
            parity_class=configuration["parity_class"],
            degree=1,
            angular_theta=configuration["assembly_angular_resolution"][0],
            angular_zeta=configuration["assembly_angular_resolution"][1],
            normalization="cas3d2mn_coefficient",
            coefficient_angular_resolution=tuple(configuration["coefficient_angular_resolution"]),
            reference_length=configuration["reference_length"],
            radial_cells=configuration["radial_cells"],
            radial_quadrature=configuration["radial_quadrature"],
        )
    fields = ("field_periods", "radial_surfaces", "quotient_rank", "labeled_nullity",
              "peak_block_width", "negative_count", "lowest_eigenvalue", "certificate",
              "eigenpair_residual", "force_balance_residual", "inertia_zero_floor",
              "coordinate_handedness", "has_eigenpair")
    run = {field: getattr(result, field) for field in fields}
    run.update({"envelope_entries": len(table),
                "labeled_sidebands": result.labeled_sideband_count,
                "envelope_modes": table,
                "library": library,
                "peak_resident_mib": peak_resident_mib(
                    resource.getrusage(resource.RUSAGE_SELF).ru_maxrss),
                "peak_resident_source_units": "bytes" if sys.platform == "darwin" else "KiB",
                "wall_seconds": time.perf_counter() - start})
    return run


def command_output(command):
    result = subprocess.run(command, capture_output=True, text=True)
    return result.stdout.strip() if result.returncode == 0 else None


def source_identity(repository):
    paths = command_output(["git", "-C", str(repository), "ls-files"])
    if paths is None:
        raise ValueError("source provenance requires a Git checkout")
    source_paths = [path for path in paths.splitlines()
                    if path == "CMakeLists.txt" or path.startswith(
                        ("src/", "include/", "python/gliss/", "cmake/"))]
    digest = hashlib.sha256()
    for relative in sorted(source_paths):
        digest.update(f"{relative}\0{sha256_file(repository / relative)}\n".encode())
    return {"commit": command_output(["git", "-C", str(repository), "rev-parse", "HEAD"]),
            "source_tree_sha256": digest.hexdigest(),
            "source_tree_hash_scheme": "SHA256 of sorted path NUL file-SHA256 LF records",
            "source_file_count": len(source_paths),
            "tracked_changes": command_output(
                ["git", "-C", str(repository), "diff", "HEAD", "--name-only"]),
            "runner_sha256": sha256_file(__file__)}


def toolchain_identity(library):
    cache = Path(library["path"]).parent / "CMakeCache.txt"
    result = {"python": sys.version, "python_executable": str(Path(sys.executable).resolve()),
              "python_executable_sha256": sha256_file(sys.executable),
              "platform": platform.platform(), "machine": platform.processor() or platform.machine()}
    if cache.is_file():
        result["adjacent_cmake_cache_sha256"] = sha256_file(cache)
        result["cmake_cache_scope"] = "informational build settings; not proof of binary provenance"
        names = ("CMAKE_Fortran_COMPILER", "CMAKE_Fortran_FLAGS", "CMAKE_Fortran_FLAGS_RELEASE",
                 "CMAKE_BUILD_TYPE", "BLAS_LIBRARIES", "LAPACK_LIBRARIES")
        settings = {}
        for line in cache.read_text().splitlines():
            if "=" in line and ":" in line.split("=", 1)[0]:
                key, value = line.split("=", 1)
                key = key.split(":", 1)[0]
                if key in names:
                    settings[key] = value
        result["cmake_settings"] = settings
        compiler = settings.get("CMAKE_Fortran_COMPILER")
        if compiler and Path(compiler).is_file():
            result["fortran_compiler_sha256"] = sha256_file(compiler)
            result["fortran_compiler_version"] = command_output([compiler, "--version"])
    return result


def write_manifest(path, manifest):
    path = Path(path)
    with tempfile.NamedTemporaryFile(mode="w", dir=path.parent, prefix=f".{path.name}.",
                                     encoding="utf-8", delete=False) as handle:
        temporary = Path(handle.name)
        try:
            json.dump(manifest, handle, indent=2, allow_nan=False)
            handle.write("\n")
        except BaseException:
            temporary.unlink(missing_ok=True)
            raise
    temporary.replace(path)


def child_command(arguments, entries):
    command = [sys.executable, str(Path(__file__).resolve()), arguments.export, "-", "--single",
               "--entries", str(entries), "--base", *map(str, arguments.base),
               "--parity", str(arguments.parity), "--angular", *map(str, arguments.angular),
               "--coefficient-angular", *map(str, arguments.coefficient_angular or arguments.angular),
               "--reference-length", str(arguments.reference_length),
               "--radial-quadrature", arguments.radial_quadrature]
    if arguments.radial_cells is not None:
        command.extend(["--radial-cells", str(arguments.radial_cells)])
    if arguments.envelope_csv is not None:
        command.extend(["--envelope-csv", str(arguments.envelope_csv),
                        "--envelope-n-sign", arguments.envelope_n_sign])
    return command


def require_frozen_inputs(arguments, manifest):
    identities = [("equilibrium input", manifest["input"]),
                  ("loaded library", manifest["library"])]
    if arguments.envelope_csv is not None:
        identities.append(("envelope CSV", manifest["envelope_csv"]))
    for label, identity in identities:
        if sha256_file(identity["path"]) != identity["sha256"]:
            raise ValueError(f"{label} changed after assumptions were frozen")
    if sha256_file(__file__) != manifest["source"]["runner_sha256"]:
        raise ValueError("runner changed after assumptions were frozen")


def main(argv=None):
    parser = argument_parser()
    arguments = parser.parse_args(argv)
    try:
        configuration, table, runs = configure(arguments, parser)
        if arguments.single:
            count = runs[0]["envelope_entries"]
            print(json.dumps(solve(arguments.export, table[:count], configuration), allow_nan=False))
            return
        import gliss

        repository = Path(__file__).resolve().parents[2]
        library = native_library_identity(gliss._load_library())
        export_identity = {"path": str(Path(arguments.export).resolve()),
                           "sha256": sha256_file(arguments.export)}
        manifest = {
            "schema_version": 2,
            "status": "running",
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "configuration": configuration,
            "export": Path(arguments.export).name,
            "export_sha256": export_identity["sha256"],
            "input": export_identity,
            "source": source_identity(repository),
            "library": library,
            "toolchain": toolchain_identity(library),
            "environment": {key: os.environ.get(key) for key in
                            ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")},
            "ordered_gliss_envelope_modes": table,
            "ordered_envelope_sha256": hashlib.sha256(
                json.dumps(table, separators=(",", ":")).encode()).hexdigest(),
            "planned_runs": runs,
            "runs": [],
        }
        # Preserve the original scaling manifest's top-level metadata keys.
        manifest.update({"gliss_commit": manifest["source"]["commit"],
                         "omp_num_threads": os.environ.get("OMP_NUM_THREADS"),
                         "machine": platform.processor() or platform.machine(),
                         "base_mode": arguments.base,
                         "angular_resolution": arguments.angular})
        if arguments.envelope_csv is not None:
            manifest["envelope_csv"] = {"path": str(arguments.envelope_csv.resolve()),
                                        "sha256": sha256_file(arguments.envelope_csv)}
        # Freeze all assumptions and input identities before observing an eigenvalue.
        write_manifest(arguments.output, manifest)
        for planned in runs:
            require_frozen_inputs(arguments, manifest)
            output = subprocess.run(child_command(arguments, planned["envelope_entries"]),
                                    capture_output=True, text=True)
            if output.returncode != 0:
                manifest["status"] = "failed"
                manifest["failure"] = {**planned, "returncode": output.returncode,
                                       "stderr": output.stderr[-16000:]}
                write_manifest(arguments.output, manifest)
                parser.exit(1, f"{planned['envelope_entries']} entries failed:\n{output.stderr}")
            require_frozen_inputs(arguments, manifest)
            lines = output.stdout.strip().splitlines()
            if not lines:
                raise ValueError("child solve did not produce a result")
            run = json.loads(lines[-1])
            if run["library"] != library:
                raise ValueError("loaded library changed after assumptions were frozen")
            run.update(planned)
            manifest["runs"].append(run)
            print(json.dumps(run, allow_nan=False), flush=True)
            write_manifest(arguments.output, manifest)
        manifest["status"] = "complete"
        write_manifest(arguments.output, manifest)
    except (OSError, ValueError) as error:
        if "manifest" in locals():
            manifest["status"] = "failed"
            manifest["failure"] = str(error)
            write_manifest(arguments.output, manifest)
        parser.exit(1, f"{error}\n")


if __name__ == "__main__":
    main()
