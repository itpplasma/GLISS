# GVEC figure-eight qualification protocol

The primary input is **N2.12-v3.1-hi**, the folded GVEC N2-12 example. Its
50,810-byte frame is identical to the [public research dataset](https://doi.org/10.5281/zenodo.14714598)
(CC-BY-4.0; G-frame and boundary attributed to G. G. Plunk in its README).
`protocol.json` pins the GVEC revision, input hashes and archive member.
`geometry-evidence.json` records an independent geometry audit: the guiding-curve
inclination is approximately 50.47 degrees, versus 13.43 degrees for the optional
QUASR189686 weak control. This is the available N2 geometry, not the approximately
76-degree maximum-fold paper case. The sampled boundary azimuth derivatives stay
positive on the sampled grids. These examples do not test toroidal-angle reversal
or prove validity of a cylindrical angular chart everywhere.

The native stored-state axes also retain the strong fold: `axis-evidence.json`
measures 49.474 degrees at M4/N10 and 49.493 degrees at M6/N20. Axis position is
independent of poloidal angle within 6e-17 m. The branches at the bowtie's
projected crossing are approximately 2.146 m apart in three dimensions. These
are stored rho=0 axis measurements, not just guiding-frame measurements;
physical equilibrium and stability remain unqualified.

**No physical acceptance fixture is promoted here.** `qualification-evidence.json`
retains failed attempts without their large state files. The fresh N2-12 attempt
at M4/N10, four degree-5 radial cells, reached coefficient force
3.49775e-9 and preserved the prescribed edge to 1.02e-14 m, but its independent
physical force RMS was 0.0165203 against the frozen 1e-4 goal. The QUASR weak
control and archived Brezellator also fail physical force balance. Increasing
the physical diagnostic's quadrature is not equilibrium resolution convergence.

A bounded continuation from the correct fresh boundary at M6/N20, six degree-5
cells and 27-by-99 angular quadrature stopped after 540 seconds/14,500 steps.
Coefficient force 8.59e-9 passes the frozen goal, but physical RMS 0.00964139
still fails. It improves on the medium resolution but worsens from the 5,000-step
checkpoint's 0.00840334 as coefficient force decreases. This does not establish
continuous-force convergence. The final state has not been exported or promoted.
On that identical state, doubling angular quadrature to 55-by-199 changes the
coefficient norm by only 1.02e-5 relative and preserves its pass. The remaining
physical discrepancy is unresolved.

`force-norm-evidence.json` records a frozen-state source and numerical check.
GVEC reports the largest Euclidean norm of constrained, preconditioned
Fourier/B-spline energy-gradient coefficients. The independent diagnostic
reports normalized volume-RMS Cartesian force. Four selected admissible
variations give physical virtual-work projections consistent with native raw
gradients within 4.98e-13 absolute; the small projections involve substantial
cancellation. The independent magnetic-energy integral agrees within 5.92e-15
relative. These checks establish neither a full physical-residual decomposition
nor a continuous-force error bound, and establish no GVEC force bug. Their probe
sources and full reports remain explicitly identified scratch artifacts.

## Prepare and attempt a solve

The executed environment was Python 3.12.13, GNU Fortran 16.1.0, OpenBLAS and
netCDF-Fortran; `requirements.txt` records its Python package versions. Install
the native development libraries first. Supply their CMake prefixes and the
compiler appropriate to the machine. From the GLISS repository root:

```sh
work=$(mktemp -d "${TMPDIR:-/tmp}/gliss-figure8.XXXXXX")
git clone https://gitlab.mpcdf.mpg.de/gvec-group/gvec.git "$work/gvec"
git -C "$work/gvec" checkout c0dc66fe2b9faa147c76a728e5cd053ce693d472
python3.12 -m venv "$work/venv"
py="$work/venv/bin/python"
"$py" -m pip install -r benchmarks/figure8/requirements.txt
export FC=gfortran-16
# Set CMAKE_PREFIX_PATH to your OpenBLAS, netCDF and netCDF-Fortran prefixes.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 CMAKE_BUILD_PARALLEL_LEVEL=2
"$py" -m pip install --no-deps "$work/gvec" \
  --config-setting="build-dir=$work/build" \
  --config-setting=cmake.define.BLA_VENDOR=OpenBLAS \
  --config-setting=cmake.define.COMPILE_GVEC_TO_GENE=On \
  --config-setting=cmake.define.COMPILE_GVEC_TO_CASTOR3D=Off \
  --config-setting=cmake.define.COMPILE_GVEC_TO_JOREK=Off \
  --config-setting=cmake.define.COMPILE_GVEC_TO_HOPR=Off
"$py" benchmarks/figure8/prepare.py "$work/gvec" "$work/prepared"
"$py" benchmarks/figure8/solve.py "$work/prepared/parameters.toml" "$work/run" \
  --iterations 20000
```

The GENE option remains enabled because this revision's configuration references
its test target when disabled. Build parallelism is bounded to two processes;
every execution script forces **one OpenMP thread before loading native code**.
GNU `libgomp` combined with Homebrew OpenBLAS `libomp` failed a separate concurrent
BLAS probe at two OpenMP threads, including with `OPENBLAS_NUM_THREADS=1`.
That variable does not constrain an OpenMP-built OpenBLAS. The GLISS runtime fix
does not fix GVEC's global `-fopenmp` configuration. Qualify a compatible GVEC
runtime before parallel solves or performance comparisons. The stored physical
failure evidence was repeated serially; it is not a backend performance claim.

Preparation defaults to a **fresh M2/N10, two degree-5-cell** initialization. This
restores the publication's boundary mode support: the current CI input's N5
truncates the edge. It retains the current revision's explicit profile convention
of pressure zero and iota +0.70; this is not a literal reproduction of an old
parameter reader. Fixed iota permits current and describes a force-free target,
not a prescribed zero-current vacuum. The archive's 400-iteration log explicitly
did not converge and supplies no native restart state or acceptance oracle.

For the executed fresh medium resolution, prepare into a different fresh directory
with `--resolution 4 10 4 --quadrature 13 49`. Each custom resolution requires an
explicit quadrature. Preserve boundary, flux and profiles across the refinement
family; measure all qualification goals independently.

## Boundary and physical checks

GVEC restart transfers boundary degrees of freedom. Raising N5 to N10 during a
restart retains the old edge; mode capacity is not evidence of preservation.
Use fresh initialization when restoring the prescribed target. `solve.py` checks
a supplied restart against the target frame **before solving**, and checks the
final edge again. A restart requires both `--restart-parameter` and
`--restart-state`; `--target-frame` defaults to `boundary_filename`.

`boundary.py FRAME PARAMETER STATE OUTPUT` independently DFT-interpolates the
source samples on shifted off-knot grids. It implements the pinned native
contracts in `hmap_axisNB.F90:195-291` and `boundaryFromFile.F90:207-240`: rotate
frame vectors into field-periodic hat components before interpolation, and apply
the input's declared `lasym=0` X-cosine/Y-sine projection. For N2 the removed
opposite-parity noise is approximately 4.46e-8/9.86e-8 m. This is not an additional
symmetry imposed by the harness. The frozen boundary tolerance is 1e-12 m;
the N5 control fails by 4.66433e-4 m. Sampled edge agreement does not establish
global embedding or nestedness.

After any completed attempt, including a reported failure, use the exact final
paths from its report:

```sh
parameterfile=$("$py" -c 'import json,sys; print(json.load(open(sys.argv[1]))["parameterfile"])' "$work/run/minimizer-report.json")
statefile=$("$py" -c 'import json,sys; print(json.load(open(sys.argv[1]))["statefile"])' "$work/run/minimizer-report.json")
"$py" benchmarks/figure8/physical_residual.py "$parameterfile" "$statefile" \
  "$work/physical.json" --toroidal-multiplier 2 --force-free
"$py" benchmarks/figure8/compare_geometry.py "$work/geometry.png" \
  "$work/prepared/N2.12-v3.1-hi.nc"
```

The physical diagnostic forms Cartesian curl from field derivatives and checks
`J × B − ∇p` independently of the coefficient minimizer norm. It integrates three
off-knot Gauss/angular refinements and checks an independent central-difference
step sweep. Normalization is the volume mean of `|B||∇|B||/mu0`;
current normalization is `volume^(1/3) mu0 sqrt(∫J²/∫B²)`. The N2 frame has high
toroidal support, so multiplier 2 produces 99/131/163 toroidal samples per period.
Use `--force-free` for N2; omit it for a prescribed zero-current vacuum.

`acceptance-goals.json` freezes coefficient/Picard goals of 1e-8 and physical RMS
goals of 1e-4. Solve, boundary and residual scripts write reports and return
status 2 on failed applicable goals. Status 0 still does not qualify a fixture:
three equilibrium resolutions, global embedding/self-contact/nestedness, a fresh
supported Cartesian export with producer pin, and GLISS parity spectra remain
required. The pinned Python current/Picard controller also has an independently
reproduced stage-loop defect with a local repair prepared, not yet published
upstream; the wrapper rejects its
unmet final goals rather than assuming that a completed loop converged.

## Bounded native GLISS compatibility

The actual fresh N2-12 state passed the native GLISS reader and spline, nine
off-knot radial surfaces on 37-by-73 angular grids, and sampled piecewise-linear
boundary admission at 25-by-64, 33-by-96 and 49-by-128 full-torus points.
`compatibility-evidence.json` pins producer patch, exporter, state, export and
native-library hashes. This is geometric ingestion evidence, not a stability
spectrum or physical-equilibrium certificate.

Apply the existing reviewed consumer exporter patch to the private producer
checkout, then export the exact final state. The wrapper checks a pointwise
Boozer-frame identity before writing the diagnostic file. For the executed N2
state the old native-angle reconstruction fails by 7.91 mm, while the corrected
Boozer-angle reconstruction agrees to 4.52e-16 m. The boundary's 1.02e-14 m
agreement describes native GVEC before finite-mode export, not export accuracy.

```sh
repo_root=$PWD
git -C "$work/gvec" apply "$repo_root/benchmarks/solovev/public/gvec-1.5.0-cas3d.patch"
"$py" benchmarks/figure8/export_compatibility.py \
  "$work/gvec/python/gvec/scripts/cas3d.py" "$parameterfile" "$statefile" \
  "$work/N2-12-UNQUALIFIED.nc"
cmake -S benchmarks/figure8/native -B "$work/native" -G Ninja \
  -DGLISS_BUILD_DIR="$repo_root/build" \
  -DFIGURE8_EXPORT="$work/N2-12-UNQUALIFIED.nc"
cmake --build "$work/native" -j2
ctest --test-dir "$work/native" --output-on-failure
```

Supply an existing verified GLISS build made with the same Fortran compiler.
The standalone probe builds only its small program and uses GLISS's compatible
OpenMP runtime selection. It requires the N2 field periods, true-edge schema,
finite regular geometry with the expected negative CAS3D Jacobian orientation,
positive field magnitude, zero pressure and consistent native edge flux slopes.
Its sampled boundary admission does not prove smooth global embedding.

## Optional controls and verification

Download the exact small QUASR JSON from the URL in `protocol.json`, then pass
`prepare.py SOURCE OUTPUT --quasr-json FILE`. Preparation verifies the SHA256,
disables cleaning and imposed symmetry, and preserves its derived M7/N11 frame.
The QUASR v4 [license record](https://doi.org/10.5281/zenodo.13717741) is MIT;
membership of this online JSON in the large archive has not been independently
verified. It is a weak-shape control, not strong-eight acceptance. Its original
OMP2 solve backend is unqualified; its serial physical failure is repeatable.
Compare frames with `compare_geometry.py OUTPUT FRAME... --labels LABEL...`.
The plot records hashes, units and shared Cartesian limits alongside the image;
the centreline is a frame guiding curve, not a certified equilibrium axis.

The portable physical diagnostic reproduces every original QUASR serial
quadrature and finite-difference value exactly. Fresh N2 preparation, a bounded
one-iteration solver smoke check, positive/negative boundary cases, and rejection
of a truncated-edge restart were executed. The geometry comparison was rendered
and inspected. The portable export matches the original NetCDF data arrays and
the standalone native CTest passes. Syntax and Ruff checks pass. No long solve was run to validate
these wrappers, and no failed states or generated images are bundled.
