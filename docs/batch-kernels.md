# Shared MHD batch kernels

The native `two_component_energy_batch`, `three_component_energy_batch` and
`perpendicular_kinetic_energy_batch` interfaces evaluate the same definitions
as their scalar counterparts. The generated `ideal_mhd_batch` module also
provides the KS3 frozen-in Hessian contraction. All variants come from
`gen/ideal_mhd/gen_ideal_mhd_batches.f90` and `ideal_mhd_symbolic.f90`.

A call is `call NAME(n, inputs, outputs)`, with contiguous structure-of-arrays
storage `inputs(n,2*m+1)` and `outputs(n,m+2)`. Columns contain the `m` physical
inputs, their `m` tangents, and one cotangent; outputs are value, JVP, and the
`m` VJP entries. FP64 input/output is the reference. FP32 inputs may have FP32
or FP64 outputs; the latter boundary preserves FP32 arithmetic and permits
FP64 caller reductions. It does not restore precision lost inside the leaf.

| Family | m | Physical input order |
|---|---:|---|
| two | 20 | Scalar two-component density argument order |
| three | 27 | Scalar three-component density argument order, followed by mu0 |
| mass | 7 | J, B, grad_s2, sigma_tilde, density, xi_s, eta |
| xi | 23 | B(3), p, mu0, A(3,3), C(3,3), matrices column-major |

The coordinate domain requires J nonzero, B and grad_s2 positive, and nonzero
flux norm; mass requires positive density. The signed-J branches are tested
separately; derivatives cannot cross J=0. No new enclosure claim is made.
The xi form is the gamma=0 plasma-volume Hessian at the identity. It supplies
no vacuum, equilibrium-boundary or complete global spectral operator.

CPU loops use OpenMP SIMD and static thread scheduling without allocation.
The CUDA header has equivalent `mhd_FAMILY_PRECISION_batch(n,x,y)` global
entry points and zero-based packed SoA storage. The caller owns streams,
launch geometry, residency, reductions and precision admission. Ordinary
GLISS builds do not require a CUDA compiler or CAS.

The standalone native oracle and timing harness is in
`benchmarks/ideal_mhd_batches`. Configure normally for CPU or add
`-DMHD_BATCH_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=120` for the RTX 5060 Ti.
`ctest` runs analytic slab/mass limits, directional differences, adjoint
pairings and precision checks. The CUDA oracle compares every output against
the CPU reference and checks tail sentinels. Timing arguments are family
(0=all, 1=two, 2=three, 3=mass, 4=xi), batch size and repetitions; the CUDA
executable takes batch size and repetitions. Each uses two warmups.

[V] GNU 12.2 cluster measurements at 262144 points, 11 samples: FP32/FP64
speedups at 1/8/24 threads are two 2.01/2.19/3.82, three 2.00/1.99/2.40,
mass 1.44/1.61/1.61, xi 2.14/2.13/2.70. FP64 24-thread scaling is
8.03/9.31/14.02/14.18 respectively. These shared-host kernel results include
value/JVP/VJP output, with FP64 reductions outside the timed kernel. They
establish neither a solver speedup nor general FP32 admission.

The caller's accepted residuals and certificates remain FP64. The new entry
points do not switch solver precision automatically. A candidate stage must
be admitted against its own sensitivity/error budget, with escalation and an
independent final residual. BF16/FP16 are storage formats supplied by the
shared precision runtime, not native arithmetic kinds in these leaves.
