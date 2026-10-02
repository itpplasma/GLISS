"""Audit an installed macOS wheel and compile a consumer of its bundled header.

Run with the fresh environment's Python, outside the source checkout, with
GLISS_LIB, PYTHONPATH and DYLD_* unset. System libraries are permitted; every
other native dependency must resolve inside the installed GLISS package.
"""

import argparse
from importlib.metadata import version
import os
from pathlib import Path
import platform
import shlex
import subprocess
import sys
import tempfile


def _output(*command):
    return subprocess.check_output(command, text=True).strip()


def _load_commands(library):
    """Read dependencies and run paths, excluding the library's own ID."""
    command = None
    for line in _output("otool", "-l", str(library)).splitlines():
        fields = line.strip().split()
        if fields[:1] == ["cmd"]:
            command = fields[1]
        elif fields[:1] in (["name"], ["path"]) and command in {
            "LC_LOAD_DYLIB", "LC_LOAD_WEAK_DYLIB", "LC_REEXPORT_DYLIB",
            "LC_LOAD_UPWARD_DYLIB", "LC_RPATH",
        }:
            yield command, line.strip().split(" (offset ")[0].split(" ", 1)[1]


def audit_bundle(package, architecture):
    libraries = sorted(package.rglob("*.dylib"))
    if not libraries:
        raise RuntimeError("the wheel contains no native libraries")
    for library in libraries:
        architectures = _output("lipo", "-archs", str(library)).split()
        if architecture not in architectures:
            raise RuntimeError(f"{library} has architectures {architectures}")
        for command, name in _load_commands(library):
            if command == "LC_RPATH":
                if not name.startswith(("@loader_path", "@executable_path")):
                    raise RuntimeError(f"{library} has external rpath {name}")
                continue
            if name.startswith(("/usr/lib/", "/System/Library/")):
                continue
            if not name.startswith("@loader_path/"):
                raise RuntimeError(f"{library} has unbundled dependency {name}")
            resolved = (library.parent / name[len("@loader_path/"):]).resolve()
            if not resolved.is_relative_to(package) or not resolved.is_file():
                raise RuntimeError(f"{library} has missing dependency {name}")
        print(f"Bundled {library.relative_to(package)}: {architecture}", flush=True)


def check_c_consumer(gliss, library, source):
    """Use only the installed header and shared library for this C program."""
    with tempfile.TemporaryDirectory(prefix="gliss-c-consumer-") as scratch:
        executable = Path(scratch) / "consumer"
        subprocess.run(
            shlex.split(os.environ.get("CC", "cc")) + [
                "-std=c11", "-Wall", "-Wextra", "-Werror",
                "-I", gliss.get_include(), str(source), str(library),
                f"-Wl,-rpath,{library.parent}", "-o", str(executable),
            ], check=True,
        )
        subprocess.run([str(executable)], cwd=scratch, check=True)
    print("Installed C-header consumer passed", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--architecture", required=True, choices=("arm64", "x86_64"))
    parser.add_argument("--c-consumer", required=True, type=Path)
    arguments = parser.parse_args()
    if platform.system() != "Darwin" or platform.machine() != arguments.architecture:
        raise RuntimeError("the runner architecture does not match this wheel")
    forbidden = [key for key in os.environ if key.startswith("DYLD_")]
    forbidden += [key for key in ("GLISS_LIB", "PYTHONPATH") if os.environ.get(key)]
    if forbidden:
        raise RuntimeError(f"wheel validation requires unset variables: {forbidden}")

    import gliss

    package = Path(gliss.__file__).resolve().parent
    if not package.is_relative_to(Path(sys.prefix).resolve()):
        raise RuntimeError(f"GLISS imported outside the fresh environment: {package}")
    library = Path(gliss._load_library()._name).resolve()
    if not library.is_relative_to(package):
        raise RuntimeError(f"GLISS loaded an external library: {library}")
    if gliss.version() != gliss.__version__ or version("gliss") != gliss.__version__:
        raise RuntimeError("the metadata, Python and native package versions differ")
    if os.environ.get("GITHUB_REF_TYPE") == "tag":
        if os.environ.get("GITHUB_REF_NAME") != f"v{gliss.__version__}":
            raise RuntimeError("the release tag does not match the package version")
    audit_bundle(package, arguments.architecture)
    check_c_consumer(gliss, library, arguments.c_consumer.resolve(strict=True))


if __name__ == "__main__":
    main()
