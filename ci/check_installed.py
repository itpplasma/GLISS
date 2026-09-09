"""Check an installed GLISS distribution using its own native library."""

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tests", type=Path, required=True)
    parser.add_argument("--equilibrium", type=Path, required=True)
    args = parser.parse_args()
    tests = args.tests.resolve(strict=True)
    equilibrium_path = args.equilibrium.resolve(strict=True)

    os.environ.pop("GLISS_LIB", None)
    import gliss

    package = Path(gliss.__file__).resolve().parent
    if not package.is_relative_to(Path(sys.prefix).resolve()):
        raise RuntimeError(f"GLISS was imported outside this environment: {package}")
    library = Path(gliss._load_library()._name).resolve()
    if not library.is_relative_to(package):
        raise RuntimeError(f"GLISS loaded an unbundled native library: {library}")
    print(f"Installed package: {package}", flush=True)
    print(f"Bundled native library: {library}", flush=True)
    assert gliss.version() == gliss.__version__
    assert (Path(gliss.get_include()) / "gliss.h").is_file()

    # The physical generalized spectrum scales as omega^2 ~ 1/rho for
    # fixed equilibrium and gamma. Exercise the complete installed ABI.
    spectra = []
    with gliss.Equilibrium(equilibrium_path) as equilibrium:
        for density, floor in ((2.0, 1.0), (8.0, 0.25)):
            with gliss.StabilityProblem(
                equilibrium,
                modes=[(1, 1), (2, 1)],
                adiabatic_index=5.0 / 3.0,
                density_kg_m3=density,
                zero_floor=floor,
                degree=1,
            ) as problem:
                spectra.append(problem.solve().classes)
    assert len(spectra[0]) == len(spectra[1]) == 2
    for original, scaled in zip(*spectra):
        assert original.has_eigenvector and scaled.has_eigenvector
        assert np.isfinite(original.lowest_eigenvalue)
        assert np.isfinite(scaled.lowest_eigenvalue)
        np.testing.assert_allclose(
            4.0 * scaled.lowest_eigenvalue,
            original.lowest_eigenvalue,
            rtol=1.0e-8,
            atol=1.0e-7,
        )
        assert original.negative_count == scaled.negative_count
    print("Installed production spectrum: inverse-density scaling passed", flush=True)

    with tempfile.TemporaryDirectory(prefix="gliss-installed-tests-") as scratch:
        destination = Path(scratch) / "tests"
        shutil.copytree(tests, destination, ignore=shutil.ignore_patterns("__pycache__"))
        environment = os.environ.copy()
        environment.pop("PYTHONPATH", None)
        environment["GLISS_LIB"] = str(library)
        subprocess.run(
            [
                sys.executable,
                "-m",
                "pytest",
                "tests",
                # This module checks source registry files, not installed behavior.
                "--ignore=tests/test_closure_registry.py",
                "-q",
                "-rs",
            ],
            cwd=scratch,
            env=environment,
            check=True,
        )


if __name__ == "__main__":
    main()
