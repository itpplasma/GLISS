# CAS3D2MN coefficient solve: scaling through L139

`run.py` solves the CAS3D2MN labeled phase envelope with the coefficient
normalization (`gliss.solve_cas3d_phase_envelope(...,
normalization="cas3d2mn_coefficient")`) for envelope tables of increasing
size, each in a fresh process, and writes a manifest of labeled and physical
dimensions, quotient rank, labeled nullity, widest block, inertia, lowest
eigenvalue, certificate, peak resident memory and wall time.

```sh
GLISS_LIB=build/libgliss_c.so PYTHONPATH=python OMP_NUM_THREADS=2 \
    python benchmarks/cas3d_envelope/run.py cd_w7x.nc l139_manifest.json
```

The labeled pencil is the congruence of the physical block-tridiagonal
pencil by a block-diagonal map, applied block by block
(`apply_cas3d_phase_envelope_block_congruence`): no dense matrix is formed
at physical or labeled size. Coincident labels give an exact null space of
the labeled pencil (`labeled_nullity`); the physical inertia is the quotient
inertia. `test_cas3d_block_congruence` checks the blockwise congruence
entrywise against the retained dense congruence on a colliding table and
the quotient rank and nullity exactly.

## Result

`l139_manifest.json`: base mode (1, 1), angular grid 64 x 64, 100 radial
surfaces, two OpenMP threads on a 4-vCPU Xeon at 2.1 GHz.

| entries | labels | quotient rank | nullity | widest block | negative | lowest | certificate | peak MiB | seconds |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 3 | 5 | 995 | 0 | 10 | 0 | 0.279558 | 8.8e-06 | 85 | 7 |
| 10 | 19 | 3781 | 0 | 38 | 0 | 0.271981 | 2.8e-05 | 91 | 21 |
| 20 | 39 | 7761 | 0 | 78 | 0 | 0.261682 | 5.5e-05 | 107 | 57 |
| 35 | 69 | 13731 | 0 | 138 | 0 | 0.258997 | 9.4e-05 | 169 | 132 |
| 70 | 139 | 27661 | 0 | 278 | 0 | 0.254767 | 1.9e-04 | 523 | 382 |

Memory grows with the widest block, not with the square of the unknown
count: a dense pencil at L139 would need 6.1 GB per matrix. The eigenvalue
decreases monotonically as sidebands are added, and the inertia stays
stable.

The equilibrium is a fixed-boundary W7-X VMEC solution (NFP 5, beta 4.5%,
101 surfaces; `wout` sha256 `f8abe4dc...c957e`) converted with
`gliss.convert_vmec` (BOOZ_XFORM 52 x 52, 100 surfaces); the manifest
records the export hash and the GLISS commit. It is not the deck of the
published CAS3D W7-X L139 study, whose angular grid, radial form functions
and reference length are not public, so the digitized CAS3D eigenvalue is
not a target for these numbers: the table measures scaling, quotient
dimensions and certified inertia, not same-deck agreement.
