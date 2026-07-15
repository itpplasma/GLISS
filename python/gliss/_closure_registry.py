"""Release registry for numerical closures and supported compositions."""

from dataclasses import dataclass
from typing import Dict, Tuple


CATEGORIES = (
    "equilibrium",
    "basis",
    "topology",
    "physics",
    "assembly",
    "norm",
    "boundary",
    "solver",
    "derivative",
)


@dataclass(frozen=True)
class Evidence:
    path: str
    token: str


@dataclass(frozen=True)
class Closure:
    identifier: str
    category: str
    implementations: Tuple[str, ...]
    evidence: Tuple[Evidence, ...]


@dataclass(frozen=True)
class Profile:
    identifier: str
    interface: str
    entrypoint: Evidence
    selection: Tuple[Tuple[str, str], ...]

    @property
    def closures(self) -> Dict[str, str]:
        return dict(self.selection)


def _evidence(path: str, token: str) -> Evidence:
    return Evidence(path, token)


CLOSURES = (
    Closure(
        "gvec-cas3d-v1",
        "equilibrium",
        ("src/gvec_cas3d_reader.f90", "src/gvec_cas3d_writer.f90"),
        (_evidence("test/test_gvec_cas3d_reader.f90", "reader_schema_error"),),
    ),
    Closure(
        "vmec-symmetric",
        "equilibrium",
        ("python/gliss/vmec.py", "python/gliss/_vmec_geometry.py"),
        (_evidence("python/tests/test_vmec.py", "test_convert_vmec"),),
    ),
    Closure(
        "terpsichore-fort23",
        "equilibrium",
        ("src/terpsichore_matrix_fixture.f90",),
        (
            _evidence(
                "test/test_terpsichore_reduced_mass_adapter.f90",
                "terpsichore_matrix_fixture_ok",
            ),
        ),
    ),
    Closure(
        "radial-feec-p1-p4",
        "basis",
        ("src/radial_feec_complex.f90",),
        (
            _evidence(
                "test/test_radial_feec_complex.f90",
                "commuting derivative identity failed",
            ),
        ),
    ),
    Closure(
        "terpsichore-p1-p0",
        "basis",
        ("src/terpsichore_reduced_layout.f90",),
        (
            _evidence(
                "test/test_terpsichore_reduced_mass_family_assembly.f90",
                "reduced family normal boundary elimination is wrong",
            ),
        ),
    ),
    Closure(
        "stellarator-parity",
        "topology",
        ("src/mode_topology.f90", "src/trial_space_topology.f90"),
        (_evidence("test/test_mode_topology.f90", "distinct families are coupled"),),
    ),
    Closure(
        "axisymmetric-family",
        "topology",
        ("src/axisymmetric_spectrum.f90",),
        (_evidence("test/test_gliss_marginality_capi.f90", "require_same_result"),),
    ),
    Closure(
        "terpsichore-mask",
        "topology",
        ("src/terpsichore_topology.f90",),
        (_evidence("test/test_terpsichore_topology.f90", "PARFAC"),),
    ),
    Closure(
        "compressible-ideal-mhd",
        "physics",
        ("src/three_component_kernel.f90",),
        (_evidence("test/test_three_component_kernel.f90", "compressible_divergence"),),
    ),
    Closure(
        "incompressible-two-component",
        "physics",
        ("src/two_component_kernel.f90",),
        (_evidence("test/test_two_component_kernel.f90", "two_component_density"),),
    ),
    Closure(
        "terpsichore-modelk0",
        "physics",
        ("src/terpsichore_noninteracting_stiffness.f90",),
        (
            _evidence(
                "test/test_terpsichore_noninteracting_stiffness.f90",
                "assemble_terpsichore_noninteracting_fixed_boundary_stiffness",
            ),
        ),
    ),
    Closure(
        "compatible-primitive-feec",
        "assembly",
        (
            "src/compatible_two_component_problem.f90",
            "src/compatible_three_component_problem.f90",
        ),
        (
            _evidence(
                "test/test_cartesian_harmonic_spline.f90",
                "compatible problem assembly failed",
            ),
        ),
    ),
    Closure(
        "terpsichore-preassembled",
        "assembly",
        ("src/terpsichore_matrix_fixture.f90",),
        (
            _evidence(
                "test/test_terpsichore_reduced_mass_adapter.f90",
                "raw TERPSICHORE MODELK is wrong",
            ),
        ),
    ),
    Closure(
        "physical-kinetic",
        "norm",
        ("src/compatible_physical_mass_assembly.f90",),
        (
            _evidence(
                "test/test_cartesian_harmonic_spline.f90",
                "compatible mass is not positive definite",
            ),
        ),
    ),
    Closure(
        "perpendicular-marginality",
        "norm",
        ("src/compatible_physical_mass_assembly.f90",),
        (
            _evidence(
                "test/test_gliss_marginality_capi.f90",
                "forced-general solve failed",
            ),
        ),
    ),
    Closure(
        "terpsichore-reduced",
        "norm",
        ("src/terpsichore_reduced_mass.f90",),
        (
            _evidence(
                "test/test_terpsichore_reduced_mass.f90",
                "assemble_terpsichore_reduced_mass_element",
            ),
        ),
    ),
    Closure(
        "feec-fixed-edge",
        "boundary",
        ("src/radial_feec_complex.f90",),
        (
            _evidence(
                "test/test_radial_feec_complex.f90",
                "right trace was not eliminated",
            ),
        ),
    ),
    Closure(
        "terpsichore-fixed-edge",
        "boundary",
        ("src/dynamic_family_layout.f90",),
        (
            _evidence(
                "test/test_dynamic_family_layout.f90",
                "edge normal coefficient was retained",
            ),
        ),
    ),
    Closure(
        "pseudoplasma-vacuum",
        "boundary",
        ("src/terpsichore_pseudoplasma_coupling.f90",),
        (
            _evidence(
                "test/test_terpsichore_pseudoplasma_stiffness.f90",
                "add_terpsichore_pseudoplasma_schur",
            ),
        ),
    ),
    Closure(
        "certified-feec-inertia",
        "solver",
        (
            "src/fixed_boundary_eigen_bracket.f90",
            "src/variable_generalized_solver.f90",
        ),
        (_evidence("test/test_fixed_boundary_spectrum.f90", "certificate"),),
    ),
    Closure(
        "dense-lapack-complete",
        "solver",
        ("src/dense_spectrum_support.f90",),
        (_evidence("test/test_fixed_boundary_spectrum.f90", "check_full_spectrum"),),
    ),
    Closure(
        "rayleigh-displacement",
        "derivative",
        ("python/gliss/derivatives.py",),
        (_evidence("python/tests/test_stability.py", "rayleigh_jvp_and_vjp"),),
    ),
)


PROFILES = (
    Profile(
        "physical-fixed-certified",
        "python",
        _evidence("python/gliss/stability.py", "def solve_class"),
        (
            ("equilibrium", "gvec-cas3d-v1"),
            ("basis", "radial-feec-p1-p4"),
            ("topology", "stellarator-parity"),
            ("physics", "compressible-ideal-mhd"),
            ("assembly", "compatible-primitive-feec"),
            ("norm", "physical-kinetic"),
            ("boundary", "feec-fixed-edge"),
            ("solver", "certified-feec-inertia"),
            ("derivative", "rayleigh-displacement"),
        ),
    ),
    Profile(
        "physical-fixed-dense-vmec",
        "python",
        _evidence("python/gliss/stability.py", "def solve_full_spectrum_class"),
        (
            ("equilibrium", "vmec-symmetric"),
            ("basis", "radial-feec-p1-p4"),
            ("topology", "stellarator-parity"),
            ("physics", "compressible-ideal-mhd"),
            ("assembly", "compatible-primitive-feec"),
            ("norm", "physical-kinetic"),
            ("boundary", "feec-fixed-edge"),
            ("solver", "dense-lapack-complete"),
            ("derivative", "rayleigh-displacement"),
        ),
    ),
    Profile(
        "axisymmetric-marginality",
        "python",
        _evidence("python/gliss/axisymmetric.py", "def solve_axisymmetric"),
        (
            ("equilibrium", "gvec-cas3d-v1"),
            ("basis", "radial-feec-p1-p4"),
            ("topology", "axisymmetric-family"),
            ("physics", "incompressible-two-component"),
            ("assembly", "compatible-primitive-feec"),
            ("norm", "perpendicular-marginality"),
            ("boundary", "feec-fixed-edge"),
            ("solver", "certified-feec-inertia"),
        ),
    ),
    Profile(
        "cas3d-marginality",
        "python",
        _evidence("python/gliss/marginality.py", "def solve_cas3d_marginality"),
        (
            ("equilibrium", "gvec-cas3d-v1"),
            ("basis", "radial-feec-p1-p4"),
            ("topology", "stellarator-parity"),
            ("physics", "incompressible-two-component"),
            ("assembly", "compatible-primitive-feec"),
            ("norm", "perpendicular-marginality"),
            ("boundary", "feec-fixed-edge"),
            ("solver", "certified-feec-inertia"),
        ),
    ),
    Profile(
        "cas3d-phase-envelope",
        "python",
        _evidence("python/gliss/marginality.py", "def solve_cas3d_phase_envelope"),
        (
            ("equilibrium", "gvec-cas3d-v1"),
            ("basis", "radial-feec-p1-p4"),
            ("topology", "stellarator-parity"),
            ("physics", "incompressible-two-component"),
            ("assembly", "compatible-primitive-feec"),
            ("norm", "perpendicular-marginality"),
            ("boundary", "feec-fixed-edge"),
            ("solver", "certified-feec-inertia"),
        ),
    ),
    Profile(
        "terpsichore-fixed",
        "python",
        _evidence(
            "python/gliss/terpsichore.py", "def solve_terpsichore_fixed_boundary"
        ),
        (
            ("equilibrium", "terpsichore-fort23"),
            ("basis", "terpsichore-p1-p0"),
            ("topology", "terpsichore-mask"),
            ("physics", "terpsichore-modelk0"),
            ("assembly", "terpsichore-preassembled"),
            ("norm", "terpsichore-reduced"),
            ("boundary", "terpsichore-fixed-edge"),
            ("solver", "certified-feec-inertia"),
        ),
    ),
    Profile(
        "terpsichore-pseudoplasma",
        "python",
        _evidence(
            "python/gliss/terpsichore.py",
            "def solve_terpsichore_pseudoplasma",
        ),
        (
            ("equilibrium", "terpsichore-fort23"),
            ("basis", "terpsichore-p1-p0"),
            ("topology", "terpsichore-mask"),
            ("physics", "terpsichore-modelk0"),
            ("assembly", "terpsichore-preassembled"),
            ("norm", "terpsichore-reduced"),
            ("boundary", "pseudoplasma-vacuum"),
            ("solver", "certified-feec-inertia"),
        ),
    ),
)


def manifest_document() -> dict:
    closures = []
    for closure in CLOSURES:
        closures.append(
            {
                "id": closure.identifier,
                "category": closure.category,
                "implementations": list(closure.implementations),
                "evidence": [
                    {"path": item.path, "token": item.token}
                    for item in closure.evidence
                ],
            }
        )
    profiles = []
    for profile in PROFILES:
        profiles.append(
            {
                "id": profile.identifier,
                "interface": profile.interface,
                "entrypoint": {
                    "path": profile.entrypoint.path,
                    "token": profile.entrypoint.token,
                },
                "closures": profile.closures,
            }
        )
    return {"schema_version": 1, "closures": closures, "profiles": profiles}
