"""Check an installed GLISS distribution using its own native library."""

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import numpy as np


def check_production_spectrum(gliss, equilibrium_path):
    # The smooth synthetic torus resolves at these angular grids. M is linear
    # in density, while K is independent of density. Check convergence and this
    # exact scaling identity through the real C ABI with non-default grids.
    spectra = []
    with gliss.Equilibrium(equilibrium_path) as equilibrium:
        for theta, zeta, density in (
            (64, 64, 2.0),
            (64, 64, 8.0),
            (128, 128, 2.0),
            (96, 128, 2.0),
            (96, 128, 8.0),
        ):
            with gliss.StabilityProblem(
                equilibrium,
                modes=[(1, 1), (2, 1)],
                adiabatic_index=5.0 / 3.0,
                density_kg_m3=density,
                zero_floor=2.0 / density,
                degree=1,
                angular_theta=theta,
                angular_zeta=zeta,
            ) as problem:
                result = problem.solve()
                assert len(result.classes) == 2
                assert problem.configuration.angular_theta == theta
                assert problem.configuration.angular_zeta == zeta
                assert gliss.StabilityConfiguration.from_dict(
                    problem.configuration.to_dict()
                ) == problem.configuration
                restored = gliss.StabilityResult.read_dict(result.to_dict())
                reference = spectra[0] if spectra else result.classes
                for original, current, persisted in zip(
                    reference, result.classes, restored.classes
                ):
                    assert current.angular_resolution == (theta, zeta)
                    assert persisted.angular_resolution == (theta, zeta)
                    assert current.has_eigenvector
                    assert np.isfinite(current.lowest_eigenvalue)
                    np.testing.assert_allclose(
                        (density / 2.0) * current.lowest_eigenvalue,
                        original.lowest_eigenvalue,
                        rtol=1.0e-8,
                        atol=1.0e-7,
                    )
                    assert original.negative_count == current.negative_count
                    # Use the same baseline vector, normalized to x^T M_2 x=1,
                    # instead of comparing vectors with different normalization.
                    energy = problem.energy(original.parity_class, original.eigenvector)
                    assert np.isfinite(energy.kinetic_energy)
                    np.testing.assert_allclose(
                        energy.kinetic_energy, density / 2.0,
                        rtol=1.0e-10, atol=1.0e-10,
                    )
                spectra.append(result.classes)
    print("Installed production spectrum: angular convergence and inverse-density "
          "scaling passed", flush=True)


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

    check_production_spectrum(gliss, equilibrium_path)

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
