# VMEC and BOOZ_XFORM input

Install the `vmec` extra with `pip install gliss[vmec]`. `convert_vmec` runs
BOOZ_XFORM on a converged VMEC `wout` file; `convert_boozer` reads a stored
`boozmn` transform. Both produce the Boozer equilibrium export GLISS reads.

The executable examples use the redistributable fixtures in `test/data/vmec`
of the GLISS source distribution or checkout. Set `GLISS_TEST_DATA` to the
`test/data` directory when running elsewhere. The sources are pinned to
SIMSOPT commit `9e027eac38028d57aa23777be52a781aa860e347`; their SHA-256
hashes, transform parameters and MIT license are stored beside the fixtures.
CI runs every Python block below against an installed wheel without downloads.

## Symmetric conversion and precomputed replay

The quasi-axisymmetric vacuum configuration of Landreman and Paul has two
field periods. Its stored transform was generated with BOOZ_XFORM 0.0.9,
including the VMEC flux profiles. Direct conversion and replay must agree.

```python
import os
import tempfile
from pathlib import Path

import numpy as np
import gliss

data = Path(os.environ.get("GLISS_TEST_DATA", "test/data")) / "vmec"
with tempfile.TemporaryDirectory() as directory:
    directory = Path(directory)
    wout = data / "wout_qa_lowres.nc"
    direct = gliss.convert_vmec(
        wout, directory / "direct.nc", poloidal_max=4, toroidal_max=3,
        radial_surfaces=7,
    )
    replay = gliss.convert_boozer(
        data / "boozmn_qa.nc", directory / "replay.nc", wout_path=wout,
        poloidal_max=4, toroidal_max=3,
    )
    spectra = []
    for path in (direct, replay):
        with gliss.Equilibrium(path) as equilibrium:
            assert equilibrium.coordinate_handedness == "left-handed"
            with gliss.StabilityProblem(
                equilibrium, [(1, 1), (2, 1)], degree=1,
                angular_theta=32, angular_zeta=32,
            ) as problem:
                assert problem.parity_classes == (1, 2)
                spectra.append(problem.solve())
assert spectra[0].lowest.lowest_eigenvalue > 0.0
np.testing.assert_allclose(
    [item.lowest_eigenvalue for item in spectra[0].classes],
    [item.lowest_eigenvalue for item in spectra[1].classes], rtol=1e-10,
)
```

## Asymmetric conversion

The Landreman–Sengupta–Plunk section 5.3 vacuum equilibrium has three field
periods and lacks stellarator symmetry. Conversion retains both Fourier
parities; the native solve uses the coupled operator, parity class 0.

```python
with tempfile.TemporaryDirectory() as directory:
    converted = gliss.convert_vmec(
        data / "wout_lsp_asymmetric.nc", Path(directory) / "asymmetric.nc",
        poloidal_max=4, toroidal_max=4, radial_surfaces=8,
    )
    with gliss.Equilibrium(converted) as equilibrium:
        with gliss.StabilityProblem(
            equilibrium, [(1, 1), (2, 1)], degree=1,
            angular_theta=32, angular_zeta=32,
        ) as problem:
            assert problem.parity_classes == (0,)
            asymmetric = problem.solve()
assert asymmetric.lowest.parity_class == 0
assert asymmetric.lowest.lowest_eigenvalue > 0.0
```

These coarse vacuum examples exercise conversion, parity and the native solve.
Physics studies require convergence in equilibrium, harmonics and radial space.

The converter rejects position harmonics whose truncated Jacobian exceeds
`truncation_tolerance`, nonpositive metric harmonics, Boozer Jacobian, field
strength or current identity residuals above 3e-2, and flux-surface averaged
radial force-balance residuals above 1e-2. It stores the residuals in the
export. Precomputed files must contain a centered uniform subset of the
VMEC half grid and all required asymmetric harmonics. The native tests also
compare BOOZ_XFORM with an independent upstream LI383 field-harmonic reference.
