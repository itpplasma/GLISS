# Conventions

- Radial coordinate: normalized toroidal flux `s`, fixed boundary at `s = 1`.
- Fourier phases: `2*pi*(m*theta - n*zeta/N_T)` for spectra.
- Boozer exports are left-handed with the pyGVEC CAS3D frame (winding -1);
  `Equilibrium.coordinate_handedness` reports the reconstructed chart.
- Eigenvalues of `StabilityProblem` are `omega^2` in s^-2 with mass
  normalization `x.T @ M @ x = 1`; negative values are unstable.
- Mercier: positive `D_Mercier` is stable.
- Regular displacements behave like `xi^s ~ s^(|m|/2)` at the axis; the
  FEEC trial space carries this factor for the normal and tangential
  components.
