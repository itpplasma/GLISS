"""Native checks of libgliss_c against independent oracles.

Every test here loads the real library; none substitutes a fake. The
oracles are the GPEC/DCON Newcomb results for the public Solov'ev fixtures
(benchmarks/solovev) and exact identities of the Rayleigh quotient.
"""

import math

import numpy as np
import pytest

import gliss

pytestmark = pytest.mark.native


def test_native_solovev_mercier_sign_matches_dcon(native_library, test_data):
    # Independent GPEC/DCON Newcomb runs of this analytic Solov'ev family near
    # q0=1.04 find the core Mercier-unstable (D_I > 0) for s below about 0.16
    # and stable outside.  Positive D_Mercier is stable in GLISS.
    path = test_data / "solovev_q1.035.nc"
    s, d_mercier = gliss.mercier_profile(path)
    assert d_mercier[0] < 0.0
    assert np.all(d_mercier[s > 0.25] > 0.0)
    crossing = s[np.argmax(d_mercier > 0.0)]
    assert 0.1 < crossing < 0.3
    assert gliss.mercier_objective(path) == pytest.approx(-d_mercier.min())
    assert gliss.mercier_objective(path) > 0.0


def test_native_gvec_export_reports_left_handed_chart(native_library, test_data):
    with gliss.Equilibrium(test_data / "solovev_q1.035.nc") as equilibrium:
        assert equilibrium.coordinate_handedness == "left-handed"


@pytest.mark.parametrize(("name", "unstable"), [("1.035", True), ("1.045", False)])
def test_native_solovev_n1_stability_matches_dcon(
    native_library, test_data, name, unstable
):
    # GPEC/DCON (tolerance 1e-8) puts the n=1 fixed-boundary marginal point of
    # this family between q0=1.0391 (unstable) and q0=1.0398 (stable). The
    # conforming axis space bounds eigenvalues from above; degree 3 resolves
    # the near-marginal kink on these 16-surface exports.
    with gliss.Equilibrium(test_data / f"solovev_q{name}.nc") as equilibrium:
        result = gliss.solve_axisymmetric(equilibrium, poloidal_max=6, degree=3)
    assert result.negative_count == int(unstable)
    assert (result.lowest_eigenvalue < 0.0) == unstable
    assert result.certificate < abs(result.lowest_eigenvalue)
    assert 0.0 < result.inertia_zero_floor < abs(result.lowest_eigenvalue)


def test_native_rayleigh_gradient_matches_finite_differences(
    native_library, test_data
):
    # The Rayleigh quotient R(x) = x.K.x / x.M.x has gradient
    # 2 (K x - R M x) / x.M.x; it matches central differences of the native
    # energy for any x and vanishes at the lowest eigenvector.
    generator = np.random.default_rng(7)
    with gliss.Equilibrium(test_data / "solovev_q1.045.nc") as equilibrium:
        with gliss.StabilityProblem(
            equilibrium,
            modes=[(0, 1), (1, -1), (1, 1)],
            degree=1,
            angular_theta=32,
            angular_zeta=8,
        ) as problem:
            size = problem._unknown_count(1)
            vector = generator.standard_normal(size)
            tangent = generator.standard_normal(size)
            gradient = problem.rayleigh_vjp(1, vector)
            step = 1.0e-6 * np.linalg.norm(vector) / np.linalg.norm(tangent)
            plus = problem.energy(1, vector + step * tangent).rayleigh_quotient
            minus = problem.energy(1, vector - step * tangent).rayleigh_quotient
            difference = (plus - minus) / (2.0 * step)
            assert problem.rayleigh_jvp(1, vector, tangent) == pytest.approx(
                difference, rel=1.0e-6
            )
            assert np.dot(gradient, tangent) == pytest.approx(difference, rel=1.0e-6)

            lowest = problem.solve().classes[0]
            stationary = problem.rayleigh_vjp(1, lowest.eigenvector)
            scale = np.linalg.norm(gradient) * np.linalg.norm(vector)
            assert np.linalg.norm(stationary) * np.linalg.norm(
                lowest.eigenvector
            ) < 1.0e-6 * scale


def test_native_coupled_operator_reproduces_parity_classes(native_library, test_data):
    # solovev_q1.035_shifted.nc is the symmetric export with the poloidal
    # angle origin moved by 3/32 (test/data/shift_poloidal_origin.py): the
    # same equilibrium stored with both parities. Its coupled operator must
    # reproduce the union of the two parity classes of the original.
    modes = [(0, 1), (1, 1), (2, 1)]
    options = {"degree": 1, "angular_theta": 32, "angular_zeta": 8}
    with gliss.Equilibrium(test_data / "solovev_q1.035.nc") as equilibrium:
        with gliss.StabilityProblem(equilibrium, modes, **options) as problem:
            assert not problem.coupled and problem.parity_classes == (1, 2)
            symmetric = problem.solve()
            union = np.sort(
                np.concatenate(
                    [problem.solve_full_spectrum_class(item).eigenvalues
                     for item in (1, 2)]
                )
            )
    with gliss.Equilibrium(test_data / "solovev_q1.035_shifted.nc") as equilibrium:
        with gliss.StabilityProblem(equilibrium, modes, **options) as problem:
            assert problem.coupled and problem.parity_classes == (0,)
            with pytest.raises(ValueError, match="parity_class must be 0"):
                problem.solve_class(1)
            coupled = problem.solve()
            full = problem.solve_full_spectrum()
            energy = problem.energy(0, coupled.classes[0].eigenvector)
    (item,) = coupled.classes
    assert item.parity_class == 0
    assert item.negative_count == sum(c.negative_count for c in symmetric.classes)
    assert item.lowest_eigenvalue == pytest.approx(
        symmetric.lowest.lowest_eigenvalue, rel=1e-9
    )
    assert energy.rayleigh_quotient == pytest.approx(item.lowest_eigenvalue, rel=1e-9)
    assert np.allclose(full.classes[0].eigenvalues, union, rtol=0.0,
                       atol=1e-9 * np.max(np.abs(union)))
    restored = gliss.StabilityResult.read_dict(coupled.to_dict())
    assert restored.classes[0].parity_class == 0


def test_native_coupled_axisymmetric_family(native_library, test_data):
    with gliss.Equilibrium(test_data / "solovev_q1.035.nc") as equilibrium:
        symmetric = gliss.solve_axisymmetric(equilibrium, poloidal_max=6, degree=3)
    with gliss.Equilibrium(test_data / "solovev_q1.035_shifted.nc") as equilibrium:
        coupled = gliss.solve_axisymmetric(equilibrium, poloidal_max=6, degree=3)
    assert (symmetric.parity_class, coupled.parity_class) == (1, 0)
    # The n = 1 kink of an axisymmetric equilibrium appears in both parities.
    assert coupled.negative_count == 2 * symmetric.negative_count == 2
    assert coupled.lowest_eigenvalue == pytest.approx(
        symmetric.lowest_eigenvalue, rel=1e-7
    )


def test_native_free_boundary_kink_and_wall(native_library, test_data):
    # q0 = 1.045 is fixed-boundary stable (DCON), but its edge carries current
    # and pressure gradient: without a wall the n=1 external kink is unstable,
    # and an ideal wall close to the edge restores stability. The fixed
    # boundary is the free space with a vanishing edge and the vacuum energy
    # grows as the wall approaches, so min-max orders the lowest eigenvalues.
    modes = [(0, 1), (1, -1), (1, 1), (2, -1), (2, 1)]
    options = dict(degree=2, angular_theta=24, angular_zeta=8)
    theta = np.linspace(0.0, 2.0 * np.pi, 24, endpoint=False)
    phi = np.linspace(0.0, 2.0 * np.pi, 12, endpoint=False)
    # A circular shell of radius 0.75 m about R = 0.935 m encloses the
    # R in [0.583, 1.288] m, |Z| < 0.565 m edge.
    radius = np.broadcast_to(0.935 + 0.75 * np.cos(theta)[:, None], (24, 12))
    shell = np.stack(
        [
            radius * np.cos(phi)[None, :],
            radius * np.sin(phi)[None, :],
            np.broadcast_to(0.75 * np.sin(theta)[:, None], radius.shape),
        ]
    )
    with gliss.Equilibrium(test_data / "solovev_q1.045.nc") as equilibrium:
        with gliss.StabilityProblem(equilibrium, modes, **options) as problem:
            fixed = problem.solve_class(1)
        lowest = {}
        for name, wall in (("none", None), ("shell", shell), ("close", 0.03)):
            vacuum = gliss.VacuumModel((24, 12), wall)
            with gliss.StabilityProblem(
                equilibrium, modes, vacuum=vacuum, **options
            ) as problem:
                assert problem.boundary_condition == "free"
                result = problem.solve_class(1)
                lowest[name] = result
                if name == "none":
                    energy = problem.energy(1, result.eigenvector)
                    configuration = problem.configuration
        document = configuration.to_dict()
        assert document["boundary_condition"] == "free"
        replayed = gliss.StabilityConfiguration.from_dict(document)
        assert replayed == configuration
        with replayed.create_problem(equilibrium) as problem:
            again = problem.solve_class(1)
    assert fixed.boundary_condition == "fixed"
    assert lowest["none"].boundary_condition == "free"
    assert fixed.negative_count == 0
    assert lowest["none"].negative_count >= 1
    assert lowest["close"].negative_count == 0
    assert (
        lowest["none"].lowest_eigenvalue
        <= lowest["shell"].lowest_eigenvalue
        <= lowest["close"].lowest_eigenvalue
        <= fixed.lowest_eigenvalue * (1.0 + 1e-10)
    )
    assert lowest["none"].normal_unknowns > fixed.normal_unknowns
    assert energy.vacuum_energy > 0.0
    assert energy.pressure_drive < 0.0
    assert energy.potential_energy == pytest.approx(math.fsum(energy.components))
    assert again.lowest_eigenvalue == lowest["none"].lowest_eigenvalue


def test_native_free_boundary_rejects_bad_vacuum(native_library, test_data):
    modes = [(1, 1), (2, 1)]
    with gliss.Equilibrium(test_data / "solovev_q1.045.nc") as equilibrium:
        with pytest.raises(gliss.GlissArgumentError, match="edge mesh"):
            gliss.StabilityProblem(
                equilibrium, modes, degree=1, angular_theta=24,
                angular_zeta=8, vacuum=gliss.VacuumModel((4, 12)),
            )
        inside = np.zeros((3, 8, 8))
        inside[0] = 0.9
        inside[1] = np.linspace(-0.01, 0.01, 8)[None, :]
        inside[2] = np.linspace(-0.01, 0.01, 8)[:, None]
        with pytest.raises(gliss.GlissError):
            gliss.StabilityProblem(
                equilibrium, modes, degree=1, angular_theta=24,
                angular_zeta=8,
                vacuum=gliss.VacuumModel((24, 12), inside),
            )
    with pytest.raises(ValueError, match="positive"):
        gliss.VacuumModel((24, 12), -0.1)
    with pytest.raises(ValueError, match="shape"):
        gliss.VacuumModel((24, 12), np.zeros((2, 4, 4)))


def test_native_radial_refinement_converges(native_library, test_data):
    # The n=1 kink of the 16-surface q0=1.035 export, solved on finite-element
    # meshes refined independently of the equilibrium surfaces. One cell per
    # surface reproduces the default mesh exactly; doubling the cells
    # converges monotonically at high order (the successive differences fall
    # by more than eight).
    with gliss.Equilibrium(test_data / "solovev_q1.035.nc") as equilibrium:
        default = gliss.solve_axisymmetric(equilibrium, poloidal_max=6, degree=3)
        lowest = []
        for cells in (16, 32, 64):
            result = gliss.solve_axisymmetric(
                equilibrium, poloidal_max=6, degree=3, radial_cells=cells
            )
            assert result.radial_surfaces == cells
            assert result.negative_count == 1
            lowest.append(result.lowest_eigenvalue)
        assert default.radial_surfaces == 16
        assert lowest[0] == default.lowest_eigenvalue
        assert lowest[0] > lowest[1] > lowest[2]
        assert lowest[1] - lowest[2] < (lowest[0] - lowest[1]) / 8.0
        marginal = gliss.cas3d_marginality_inertia(
            equilibrium, [(1, 1), (2, 1)], degree=1, radial_cells=24
        )
        assert marginal.radial_surfaces == 24
        envelope = gliss.cas3d_phase_envelope_inertia(
            equilibrium, (1, 1), [(0, 0)], degree=1, radial_cells=24
        )
        assert envelope.radial_surfaces == 24
        coefficient = gliss.cas3d_phase_envelope_inertia(
            equilibrium, (1, 1), [(0, 0)], degree=1, radial_cells=24,
            normalization="cas3d2mn_coefficient",
            coefficient_angular_resolution=(16, 16), reference_length=1.0,
        )
        assert coefficient.radial_surfaces == 24


def test_native_stability_problem_radial_cells(native_library, test_data, tmp_path):
    modes = [(1, 1), (2, 1)]
    with gliss.Equilibrium(test_data / "solovev_q1.045.nc") as equilibrium:
        with gliss.StabilityProblem(
            equilibrium, modes, degree=1, angular_theta=24, angular_zeta=8
        ) as coarse:
            coarse_results = coarse.solve()
            coarse_result = coarse_results.classes[0]
        with pytest.raises(ValueError, match="equilibrium SHA-256"):
            gliss.write_run_manifest(
                tmp_path / "wrong-equilibrium.json",
                test_data / "solovev_q1.035.nc", coarse.configuration, coarse_results,
            )
        with gliss.StabilityProblem(
            equilibrium, modes, degree=1, angular_theta=24, angular_zeta=8,
            radial_cells=32,
        ) as fine:
            fine_result = fine.solve_class(1)
            configuration = fine.configuration
            with pytest.raises(ValueError, match="configuration SHA-256"):
                fine.write_manifest(tmp_path / "wrong-mesh.json", coarse_results)
        with gliss.StabilityProblem(
            equilibrium, modes, degree=1, angular_theta=24, angular_zeta=8,
            radial_cells=32, vacuum=gliss.VacuumModel((24, 12)),
        ) as free:
            assert free.boundary_condition == "free"
            assert free.configuration.radial_cells == 32
            free_results = free.solve()
            free_result = free_results.classes[0]
            manifest = free.write_manifest(tmp_path / "free.json", free_results)
            assert manifest.configuration_verified
            assert manifest.equilibrium_verified
            tampered = manifest.to_dict()
            tampered["equilibrium"]["sha256"] = "0" * 64
            with pytest.raises(ValueError, match="equilibrium SHA-256"):
                gliss.RunManifest.from_dict(tampered)
            changed_wall = gliss.StabilityConfiguration(
                modes, degree=1, angular_theta=24, angular_zeta=8,
                radial_cells=32, vacuum=gliss.VacuumModel((24, 12), 0.03),
            )
            with pytest.raises(ValueError, match="configuration SHA-256"):
                gliss.write_run_manifest(
                    tmp_path / "wrong-wall.json", equilibrium.path,
                    changed_wall, free_results,
                )
        replayed = configuration.create_problem(equilibrium)
        with replayed:
            assert replayed.solve_class(1).lowest_eigenvalue == (
                fine_result.lowest_eigenvalue
            )
        for cells in (0, 1, -4, True, 2.0):
            with pytest.raises((TypeError, ValueError), match="radial_cells"):
                gliss.StabilityProblem(equilibrium, modes, radial_cells=cells)
            with pytest.raises((TypeError, ValueError), match="radial_cells"):
                gliss.solve_axisymmetric(equilibrium, radial_cells=cells)
    # Twice the cells: twice the normal unknowns away from the axis.
    assert fine_result.normal_unknowns > 1.9 * coarse_result.normal_unknowns
    # The edge is free: the vacuum lowers the energy of the fixed edge.
    assert free_result.lowest_eigenvalue < fine_result.lowest_eigenvalue
    path = tmp_path / "configuration.json"
    configuration.write(path)
    assert gliss.StabilityConfiguration.read(path) == configuration
    assert configuration.to_dict()["radial_cells"] == 32
