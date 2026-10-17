# QAS3 current-regularization comparison

The [source audit](source-audit.md) traces the pinned public TERPSICHORE and
GLISS current-reconstruction paths into stiffness. Their present QAS3 recipe
uses different resonance regularizations: TERPSICHORE's fixed `DELTAJP=0.04`
and GLISS's width from the equilibrium-data rotational-transform variation.
This identifies a comparison that needs isolation, without establishing a
defect or explaining the historical unstable-mode counts. Historical executed
deck hashes have not been recovered.

The [derivation](derivation.md) gives both responses in a common chart,
their domains, force-balance residuals, exact difference and singular limits.
The standalone native CMake/CTest suite proves 27 rational identities with
both native FortSym and SymEngine: 54 strict checks, zero probes. The generated
kernel CI runs it using the repository's pinned FortSym.

```sh
cmake -S benchmarks/qas3/symbolic_current/src \
  -B /tmp/gliss-current-response-proof -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/path/to/fortsym-prefix
cmake --build /tmp/gliss-current-response-proof --parallel 2
ctest --test-dir /tmp/gliss-current-response-proof --output-on-failure
```

The frozen [source receipt](source-receipt.json) and [proof receipt](receipt.json)
bind original sources and results. Build-log paths refer to the isolated
execution; the source report's `receipt.json` refers to `source-receipt.json`
here. No numerical spectra were run. The next control must qualify a common
equilibrium, then compare the two explicit responses and trial energies with
all other settings fixed. Radial stability-mesh refinement alone does not
make these regularizations equivalent.
