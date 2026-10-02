# Native VMEC conversion fixtures

These files let CI run BOOZ_XFORM and GLISS's native readers and spectra without
network access. Upstream files come from SIMSOPT commit
`9e027eac38028d57aa23777be52a781aa860e347`; `provenance.json` records the
original paths, complete SHA-256 hashes and generated transform parameters.
The upstream repository's MIT license is retained as `LICENSE.simsopt`.

`wout_qa_lowres.nc` is the two-field-period Landreman–Paul quasi-axisymmetric
vacuum equilibrium. `wout_lsp_asymmetric.nc` is the three-field-period,
non-stellarator-symmetric vacuum equilibrium from Landreman, Sengupta and
Plunk, JPP 85, 905850103 (2019), section 5.3. Both pass GLISS's default
Jacobian, current, truncation and averaged force-balance admission gates at
the resolutions recorded in the manifest. The native tests compare direct
conversion, stored BOOZ_XFORM replay and native NetCDF round trips, including
the coupled asymmetric spectrum.

The two generated `boozmn` fixtures use BOOZ_XFORM 0.0.9. Reproduce them by
reading the corresponding `wout` with `read_wout(path, True)`, setting
`mboz`, `nboz` and `compute_surfs` from the manifest, then calling `run()` and
`write_boozmn(path)`. The second `read_wout` argument reads the flux profiles
needed by GLISS. Export QA with GLISS harmonics (4, 3) and seven radial
surfaces; export the asymmetric case with (4, 4) and eight surfaces.

`wout_li383_lowres.nc` and `boozmn_li383_reference.nc` are unchanged upstream
files. The reference predates GLISS; its first and last surfaces provide an
independent BOOZ_XFORM field-harmonic oracle at absolute and relative tolerance
1e-12. LI383's low-resolution geometry does not pass GLISS's Jacobian gate
(relative maximum about 0.086 with all 15 half-grid surfaces), so it supplies
the transform oracle rather than an admitted GLISS equilibrium.

Original names and physical descriptions are in the
[pinned upstream fixture README](https://github.com/hiddenSymmetries/simsopt/blob/9e027eac38028d57aa23777be52a781aa860e347/tests/test_files/README.md).
The [upstream reference test](https://github.com/hiddenSymmetries/simsopt/blob/9e027eac38028d57aa23777be52a781aa860e347/tests/mhd/test_boozer.py)
specifies the LI383 coefficient comparison. These fixtures establish input
conversion and representation consistency; they are not converged cross-code
stability benchmarks.
