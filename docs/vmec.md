# VMEC and BOOZ_XFORM input

`convert_vmec` runs `booz_xform` on a converged VMEC `wout` file and writes
the Boozer export GLISS reads; `convert_boozer` reads a precomputed `boozmn`
file instead and produces the same export for the same transform. Both need
the `vmec` extra (`pip install gliss[vmec]`).

```text
converted = gliss.convert_vmec("wout_W7X.nc", "W7X_gliss.nc",
                               poloidal_max=12, toroidal_max=12)
converted = gliss.convert_boozer("boozmn_W7X.nc", "W7X_gliss.nc",
                                 wout_path="wout_W7X.nc")
```

The conversion rejects exports whose truncated position harmonics reproduce
the Jacobian worse than `truncation_tolerance`, whose metric harmonics are
not positive definite, whose Boozer Jacobian, |B| or current identities
disagree with the geometry by more than 3e-2, or whose flux-surface-averaged
radial force balance fails by more than 1e-2. The residuals are stored as
attributes of the export. Asymmetric (`lasym`) files are converted with both
parities and solved with the coupled parity operator (parity class 0). See
the Python README for the full contract.
