# GLISS

GLISS (Global Linear Ideal Stability Solver) computes the linear ideal-MHD
stability of three-dimensional toroidal equilibria with nested flux surfaces.
It solves the energy-principle eigenvalue problem `K x = omega^2 M x` with
Fourier harmonics in the angles and spline finite elements in the radius, and
is designed for use in differentiable stellarator optimization.

Status: research code under active development; the Python API may change
between 0.0.x releases.

## Features

- Fixed- and free-boundary spectra and energies. The free-boundary problem
  couples the plasma to a current-free vacuum with an optional ideal
  conducting wall.
- Lowest-eigenpair solves with a residual-based error bound, optional full
  spectra, and eigenvalue counting by matrix inertia.
- Mercier interchange stability profiles.
- Equilibrium input from the [GVEC](https://gitlab.mpcdf.mpg.de/gvec-group/gvec)
  CAS3D export, from VMEC `wout` files and from precomputed BOOZ_XFORM files.
- Pressure derivatives of eigenvalues at fixed equilibrium geometry; assembly
  kernels with Enzyme-generated derivatives. The full equilibrium-to-spectrum
  derivative chain is not yet complete.
- Compatibility paths that reproduce TERPSICHORE's stored discretization from
  FORT.23/FORT.24 files, for validation only.

The free-boundary model assumes zero edge pressure and a tangential
equilibrium field that matches a current-free exterior, with no equilibrium
surface current. Users must check that their equilibrium satisfies this.

## Installation

The Python package is the main interface. Release 0.0.2 is on PyPI:

```sh
python -m pip install gliss
```

The development version (0.0.3, documented on the `latest` docs site) installs
from a checkout with `python -m pip install .`, or `python -m pip install
".[vmec]"` for VMEC conversion. Building from source needs CMake, Ninja, a
Fortran compiler, BLAS/LAPACK, PkgConfig and the NetCDF C library. To build the
native library and run the tests:

```sh
cmake -S . -B build -G Ninja
cmake --build build
ctest --test-dir build --output-on-failure
```

See [docs/building.md](docs/building.md) for the release audit, Enzyme
derivative checks and the documentation build.

## Example

Count unstable n = 1 directions of the Solov'ev test equilibria in `test/data`:

```python
import gliss

for path in ("test/data/solovev_q1.035.nc", "test/data/solovev_q1.045.nc"):
    with gliss.Equilibrium(path) as equilibrium:
        result = gliss.solve_axisymmetric(equilibrium, poloidal_max=6, degree=3)
    print(path, result.negative_count, result.lowest_eigenvalue)
```

The first equilibrium (q0 = 1.035) has one unstable direction, the second
(q0 = 1.045) none.

## Documentation

The manual is at <https://itpplasma.github.io/GLISS/>, with a version menu
for releases. Sources are in [docs/](docs/index.md):

- [Quickstart](docs/quickstart.md) and [VMEC input](docs/vmec.md)
- [Free-boundary problem](docs/free_boundary.md) and
  [TERPSICHORE compatibility](docs/compatibility.md)
- [Conventions](docs/conventions.md) and [API reference](docs/api.md)
- [Mode diagnostics](docs/mode_diagnostics.md), including command-line
  eigenvalue counts and eigenfunction profiles

The [Python guide](python/README.md) covers input and output formats, mode
tables and norms, and the optional SIMSOPT adapter. Planned work is in
[ROADMAP.md](ROADMAP.md), changes in [CHANGELOG.md](CHANGELOG.md), and the
merge process in [CONTRIBUTING.md](CONTRIBUTING.md).

## Formulation

The formulation follows the CAS3D energy-principle work of Carolin Schwab
(later Nührenberg). Further methods derive from Bernstein et al. (1958) for the
energy principle, Newcomb (1960) and Suydam (1958) for cylindrical tests,
Mercier (1960) and Landreman and Jorge (2020) for the interchange criterion,
and Anderson et al. (1990) for eigenvalue counting by matrix inertia.
[PROVENANCE.md](PROVENANCE.md) maps each module to its sources.

## Citing

Citation metadata is in [CITATION.cff](CITATION.cff).

## License

MIT. See [LICENSE](LICENSE).
