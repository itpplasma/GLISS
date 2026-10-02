import hashlib
import json
from dataclasses import replace

import numpy as np
import pytest

import gliss


def spectrum(parity_class, lowest):
    vector = np.arange(6, dtype=np.float64) + parity_class
    vector.setflags(write=False)
    return gliss.SpectrumResult(
        parity_class=parity_class,
        field_periods=3,
        modes=((1, 1), (2, -1)),
        degree=2,
        angular_resolution=(64, 64),
        adiabatic_index=5.0 / 3.0,
        density_kg_m3=2.0,
        zero_floor=1.0,
        negative_count=1,
        floor_count=2,
        lowest_eigenvalue=lowest,
        certificate=6.0e-10,
        eigenpair_residual=1.0e-10,
        eigenpair_resolution=2.0e-10,
        inertia_interval=3.0e-10,
        eigenvector=vector,
        normal_unknowns=2,
        eta_unknowns=2,
        mu_unknowns=2,
        has_chart_metric=True,
        has_eigenvector=True,
        configuration_sha256=gliss.StabilityConfiguration(
            modes=((1, 1), (2, -1)), density_kg_m3=2.0
        ).sha256,
    )


@pytest.fixture
def configuration():
    return gliss.StabilityConfiguration(
        modes=[(1, 1), (2, -1)],
        adiabatic_index=5.0 / 3.0,
        density_kg_m3=2.0,
        zero_floor=1.0,
        degree=2,
    )


@pytest.fixture
def result():
    return gliss.StabilityResult((spectrum(1, -5.0), spectrum(2, -4.0)))


def test_configuration_round_trip_is_deterministic(configuration, tmp_path):
    first = tmp_path / "configuration.json"
    second = tmp_path / "configuration-copy.json"
    configuration.write(first)
    loaded = gliss.StabilityConfiguration.read(first)
    loaded.write(second)

    assert loaded == configuration
    assert first.read_bytes() == second.read_bytes()
    document = json.loads(first.read_text(encoding="utf-8"))
    assert document["schema"] == "gliss.stability.configuration"
    assert document["schema_version"] == 8
    assert document["discretization_revision"] == 4
    assert document["boundary_condition"] == "fixed"


def test_schema_three_round_trip_preserves_solver_tolerances(configuration, result):
    tolerances = gliss.SolverTolerances(
        eigenvalue_relative=2.0e-13,
        residual_relative=3.0e-12,
        negative_bracket_relative=4.0e-9,
        negative_bracket_floor=5.0e-3,
        inverse_iteration_limit=321,
        bracket_iteration_limit=123,
    )
    configured = replace(configuration, solver_tolerances=tolerances)
    configured_document = configured.to_dict()
    assert configured_document["schema_version"] == 8
    assert gliss.StabilityConfiguration.from_dict(configured_document) == configured

    controlled_result = gliss.StabilityResult(
        tuple(replace(item, solver_tolerances=tolerances) for item in result.classes)
    )
    result_document = controlled_result.to_dict()
    assert result_document["schema_version"] == 8
    loaded = gliss.StabilityResult.read_dict(result_document)
    assert all(item.solver_tolerances == tolerances for item in loaded.classes)


def test_schema_one_reads_historical_solver_tolerances(configuration, result):
    configuration_document = configuration.to_dict()
    configuration_document["schema_version"] = 1
    configuration_document.pop("angular_theta")
    configuration_document.pop("angular_zeta")
    configuration_document["radial_quadrature"] = "midpoint"
    configuration_document.pop("degree")
    configuration_document.pop("solver_tolerances")
    configuration_document.pop("discretization_revision")
    configuration_document.pop("vacuum")
    configuration_document.pop("radial_cells")
    result_document = result.to_dict()
    result_document["schema_version"] = 1
    for item in result_document["classes"]:
        item["radial_quadrature"] = "midpoint"
        item.pop("degree")
        item.pop("solver_tolerances")
        item.pop("discretization_revision")
        item.pop("configuration_sha256")
        item.pop("equilibrium_sha256")
    loaded_configuration = gliss.StabilityConfiguration.from_dict(
        configuration_document
    )
    loaded_result = gliss.StabilityResult.read_dict(result_document)
    defaults = gliss.SolverTolerances.historical_defaults()
    assert loaded_configuration.degree == 1
    assert all(item.degree == 1 for item in loaded_result.classes)
    assert loaded_configuration.solver_tolerances == defaults
    assert all(item.solver_tolerances == defaults for item in loaded_result.classes)
    # Midpoint records are readable but cannot be replayed on today's operator.
    assert loaded_configuration.discretization_revision == 0
    assert all(item.discretization_revision == 0 for item in loaded_result.classes)
    with pytest.raises(ValueError, match="operator changed.*revision 0"):
        loaded_configuration.create_problem(object())


@pytest.mark.parametrize(
    ("change", "match"),
    [
        (lambda value: value.update(extra=1), "unknown field 'extra'"),
        (lambda value: value.pop("modes"), "missing field 'modes'"),
        (
            lambda value: value.update(schema_version=9),
            "schema_version.*expected 1 or 2 or 3",
        ),
        (lambda value: value.update(boundary_condition="free"), "fixed"),
        (lambda value: value.update(density_kg_m3=float("nan")), "finite"),
        (lambda value: value.update(modes=[[1, 1], [1, 1]]), "duplicate"),
    ],
)
def test_configuration_rejects_invalid_documents(
    configuration, tmp_path, change, match
):
    path = tmp_path / "configuration.json"
    configuration.write(path)
    document = json.loads(path.read_text(encoding="utf-8"))
    change(document)
    path.write_text(json.dumps(document), encoding="utf-8")

    with pytest.raises(ValueError, match=match):
        gliss.StabilityConfiguration.read(path)


def test_configuration_rejects_truncated_json(tmp_path):
    path = tmp_path / "configuration.json"
    path.write_text('{"schema":', encoding="utf-8")
    with pytest.raises(ValueError, match="invalid or truncated JSON.*line 1"):
        gliss.StabilityConfiguration.read(path)


@pytest.mark.parametrize("name", ["schema", "schema_version"])
def test_configuration_rejects_missing_schema_fields(configuration, name):
    document = configuration.to_dict()
    del document[name]
    with pytest.raises(ValueError, match=f"missing field '{name}'"):
        gliss.StabilityConfiguration.from_dict(document)


def test_configuration_rejects_duplicate_json_fields(configuration, tmp_path):
    path = tmp_path / "configuration.json"
    configuration.write(path)
    text = path.read_text(encoding="utf-8")
    text = text.replace(
        '"schema": "gliss.stability.configuration",',
        '"schema": "gliss.stability.configuration",\n'
        '  "schema": "gliss.stability.configuration",',
        1,
    )
    path.write_text(text, encoding="utf-8")
    with pytest.raises(ValueError, match="duplicate field 'schema'"):
        gliss.StabilityConfiguration.read(path)


def test_result_round_trip_preserves_metadata_and_vectors(result, tmp_path):
    first = tmp_path / "result.json"
    second = tmp_path / "result-copy.json"
    result.write(first)
    loaded = gliss.StabilityResult.read(first)
    loaded.write(second)

    assert first.read_bytes() == second.read_bytes()
    assert loaded.lowest.parity_class == 1
    assert loaded.classes[0].modes == result.classes[0].modes
    assert loaded.classes[0].certificate == result.classes[0].certificate
    np.testing.assert_array_equal(
        loaded.classes[0].eigenvector, result.classes[0].eigenvector
    )
    assert not loaded.classes[0].eigenvector.flags.writeable


@pytest.mark.parametrize(
    ("change", "match"),
    [
        (lambda value: value.update(extra=1), "unknown field 'extra'"),
        (lambda value: value.pop("classes"), "missing field 'classes'"),
        (lambda value: value.update(schema_version=9), "schema_version.*expected 1"),
        (
            lambda value: value["classes"][0].update(lowest_eigenvalue=float("inf")),
            "lowest_eigenvalue.*finite",
        ),
        (
            lambda value: value["classes"][0].update(eigenvector=[1.0]),
            "eigenvector.*expected 6",
        ),
        (
            lambda value: value["classes"][1].update(parity_class=1),
            "parity classes.*1 then 2",
        ),
    ],
)
def test_result_rejects_invalid_documents(result, tmp_path, change, match):
    path = tmp_path / "result.json"
    result.write(path)
    document = json.loads(path.read_text(encoding="utf-8"))
    change(document)
    path.write_text(json.dumps(document), encoding="utf-8")

    with pytest.raises(ValueError, match=match):
        gliss.StabilityResult.read(path)


def test_result_writer_rejects_inconsistent_objects(result, tmp_path):
    invalid_class = replace(result.classes[0], certificate=float("nan"))
    invalid = gliss.StabilityResult((invalid_class, result.classes[1]))
    with pytest.raises(ValueError, match="certificate.*finite"):
        invalid.write(tmp_path / "result.json")


def test_run_manifest_is_portable_and_round_trips(
    configuration, result, tmp_path, monkeypatch
):
    equilibrium = tmp_path / "private" / "equilibrium.nc"
    equilibrium.parent.mkdir()
    equilibrium.write_bytes(b"equilibrium fixture")
    path = tmp_path / "run.json"
    monkeypatch.setattr("gliss.schema._native_version", lambda: "0.0.1")
    monkeypatch.setattr("gliss.schema._equilibrium_schema_version", lambda path: 0)

    manifest = gliss.write_run_manifest(path, equilibrium, configuration, result)
    loaded = gliss.RunManifest.read(path)
    text = path.read_text(encoding="utf-8")

    assert loaded == manifest
    assert loaded.equilibrium_filename == "equilibrium.nc"
    assert loaded.equilibrium_schema_version == 0
    assert (
        loaded.equilibrium_sha256
        == hashlib.sha256(equilibrium.read_bytes()).hexdigest()
    )
    assert str(equilibrium.parent) not in text
    assert loaded.configuration == configuration
    np.testing.assert_array_equal(
        loaded.result.classes[1].eigenvector, result.classes[1].eigenvector
    )


def test_manifest_rejects_result_configuration_mismatch(
    configuration, result, tmp_path, monkeypatch
):
    equilibrium = tmp_path / "equilibrium.nc"
    equilibrium.touch()
    monkeypatch.setattr("gliss.schema._native_version", lambda: "0.0.1")
    monkeypatch.setattr("gliss.schema._equilibrium_schema_version", lambda path: 0)
    mismatch = replace(configuration, modes=((3, 2),))
    with pytest.raises(ValueError, match="result modes.*configuration modes"):
        gliss.write_run_manifest(tmp_path / "run.json", equilibrium, mismatch, result)


@pytest.mark.parametrize("configuration_is_free", [False, True])
def test_manifest_rejects_boundary_mismatch(
    configuration, result, tmp_path, configuration_is_free
):
    # The boundary changes the admissible displacement space and operator.
    # A persisted run must never advertise the opposite boundary on replay.
    if configuration_is_free:
        configuration = replace(configuration, vacuum=gliss.VacuumModel())
    else:
        result = gliss.StabilityResult(
            tuple(replace(item, boundary_condition="free") for item in result.classes)
        )
    destination = tmp_path / "run.json"
    with pytest.raises(ValueError, match="result boundary_condition.*configuration"):
        gliss.write_run_manifest(
            destination, tmp_path / "equilibrium.nc", configuration, result
        )
    assert not destination.exists()


def test_historical_manifest_remains_explicitly_unverified(
    configuration, result, tmp_path, monkeypatch
):
    equilibrium = tmp_path / "equilibrium.nc"
    equilibrium.write_bytes(b"fixture")
    monkeypatch.setattr("gliss.schema._native_version", lambda: "0.0.2")
    monkeypatch.setattr("gliss.schema._equilibrium_schema_version", lambda path: 1)
    manifest = gliss.write_run_manifest(
        tmp_path / "run.json", equilibrium, configuration, result
    )
    assert manifest.configuration_verified
    document = manifest.to_dict()
    document["schema_version"] = 7
    document["configuration"]["schema_version"] = 7
    document["result"]["schema_version"] = 7
    for item in document["result"]["classes"]:
        del item["configuration_sha256"]
        del item["equilibrium_sha256"]
    historical = gliss.RunManifest.from_dict(document)
    assert not historical.configuration_verified
    assert not historical.equilibrium_verified
    assert all(item.configuration_sha256 is None for item in historical.result.classes)
    historical.write(tmp_path / "historical.json")
    assert not gliss.RunManifest.read(tmp_path / "historical.json").configuration_verified
    with pytest.raises(ValueError, match="provenance is unavailable"):
        gliss.write_run_manifest(
            tmp_path / "invented.json", equilibrium, configuration, historical.result
        )
    assert not (tmp_path / "invented.json").exists()


def test_manifest_rejects_changed_wall_nodes(configuration, result, tmp_path):
    wall = np.zeros((3, 4, 4))
    original = replace(configuration, vacuum=gliss.VacuumModel((8, 8), wall))
    recorded = gliss.StabilityResult(
        tuple(
            replace(item, boundary_condition="free", configuration_sha256=original.sha256)
            for item in result.classes
        )
    )
    wall[0, 0, 0] = np.nextafter(0.0, 1.0)
    changed = replace(original, vacuum=gliss.VacuumModel((8, 8), wall))
    with pytest.raises(ValueError, match="configuration SHA-256"):
        gliss.write_run_manifest(
            tmp_path / "wrong-nodes.json", tmp_path / "equilibrium.nc", changed, recorded
        )


def test_manifest_rejects_equilibrium_changed_during_metadata_collection(
    configuration, result, tmp_path, monkeypatch
):
    equilibrium = tmp_path / "equilibrium.nc"
    equilibrium.write_bytes(b"before")
    monkeypatch.setattr("gliss.schema._native_version", lambda: "0.0.1")

    def change_source(path):
        path.write_bytes(b"after")
        return 0

    monkeypatch.setattr("gliss.schema._equilibrium_schema_version", change_source)
    with pytest.raises(gliss.GlissIOError, match="changed while collecting"):
        gliss.write_run_manifest(
            tmp_path / "run.json", equilibrium, configuration, result
        )
    assert not (tmp_path / "run.json").exists()


def test_manifest_rejects_equilibrium_deleted_during_metadata_collection(
    configuration, result, tmp_path, monkeypatch
):
    equilibrium = tmp_path / "equilibrium.nc"
    equilibrium.write_bytes(b"before")
    monkeypatch.setattr("gliss.schema._native_version", lambda: "0.0.1")

    def remove_source(path):
        path.unlink()
        return 0

    monkeypatch.setattr("gliss.schema._equilibrium_schema_version", remove_source)
    with pytest.raises(gliss.GlissIOError, match="changed while collecting"):
        gliss.write_run_manifest(
            tmp_path / "run.json", equilibrium, configuration, result
        )
    assert not (tmp_path / "run.json").exists()


def test_manifest_rejects_changed_equilibrium(configuration, result, tmp_path):
    equilibrium = tmp_path / "equilibrium.nc"
    equilibrium.write_bytes(b"before")
    path = tmp_path / "run.json"
    document = {
        "schema": "gliss.stability.run",
        "schema_version": 8,
        "equilibrium": {
            "format": "gvec-cas3d-netcdf",
            "schema_version": 1,
            "filename": equilibrium.name,
            "size_bytes": 6,
            "sha256": "0" * 64,
        },
        "software": {
            "gliss_python": "0.0.1",
            "gliss_native": "0.0.1",
            "gliss_abi": 1,
            "numpy": np.__version__,
            "python": "3.9.0",
        },
        "configuration": configuration.to_dict(),
        "result": result.to_dict(),
    }
    path.write_text(json.dumps(document), encoding="utf-8")
    loaded = gliss.RunManifest.read(path)
    with pytest.raises(ValueError, match="SHA-256 does not match"):
        loaded.verify_equilibrium(equilibrium)


def test_manifest_constructor_rejects_private_or_invalid_metadata(
    configuration, result
):
    with pytest.raises(ValueError, match="filename.*base name"):
        gliss.RunManifest(
            equilibrium_filename="/private/equilibrium.nc",
            equilibrium_size_bytes=1,
            equilibrium_sha256="0" * 64,
            equilibrium_schema_version=0,
            configuration=configuration,
            result=result,
            gliss_python_version="0.0.1",
            gliss_native_version="0.0.1",
            gliss_abi_version=1,
            numpy_version=np.__version__,
            python_version="3.9.0",
        )


@pytest.mark.parametrize("version", [-1, 2, True, "1"])
def test_manifest_rejects_invalid_equilibrium_schema_version(
    configuration, result, version
):
    with pytest.raises((TypeError, ValueError), match="equilibrium_schema_version"):
        gliss.RunManifest(
            equilibrium_filename="equilibrium.nc",
            equilibrium_size_bytes=1,
            equilibrium_sha256="0" * 64,
            equilibrium_schema_version=version,
            configuration=configuration,
            result=result,
            gliss_python_version="0.0.1",
            gliss_native_version="0.0.1",
            gliss_abi_version=1,
            numpy_version=np.__version__,
            python_version="3.9.0",
        )


def test_manifest_reader_accepts_legacy_and_rejects_unknown_equilibrium_schema(
    configuration, result, tmp_path
):
    path = tmp_path / "run.json"
    document = {
        "schema": "gliss.stability.run",
        "schema_version": 8,
        "equilibrium": {
            "format": "gvec-cas3d-netcdf",
            "schema_version": 0,
            "filename": "equilibrium.nc",
            "size_bytes": 0,
            "sha256": "0" * 64,
        },
        "software": {
            "gliss_python": "0.0.1",
            "gliss_native": "0.0.1",
            "gliss_abi": 1,
            "numpy": np.__version__,
            "python": "3.9.0",
        },
        "configuration": configuration.to_dict(),
        "result": result.to_dict(),
    }
    path.write_text(json.dumps(document), encoding="utf-8")
    assert gliss.RunManifest.read(path).equilibrium_schema_version == 0

    document["equilibrium"]["schema_version"] = 2
    path.write_text(json.dumps(document), encoding="utf-8")
    with pytest.raises(ValueError, match="equilibrium.schema_version.*0 or 1"):
        gliss.RunManifest.read(path)


def test_custom_angular_configuration_round_trip(tmp_path):
    config = gliss.StabilityConfiguration(modes=((1, 0),), angular_theta=128, angular_zeta=96)
    path = tmp_path / "quadrature.json"
    config.write(path)
    loaded = gliss.StabilityConfiguration.read(path)
    assert loaded == config
    legacy = config.to_dict()
    legacy["schema_version"] = 3
    del legacy["angular_theta"], legacy["angular_zeta"]
    del legacy["discretization_revision"], legacy["vacuum"]
    del legacy["radial_cells"]
    restored = gliss.StabilityConfiguration.from_dict(legacy)
    assert restored.angular_theta == 64
    # Schema 3 and 4 records predate the axis-conforming FEEC space (#16).
    assert restored.discretization_revision == 1
    with pytest.raises(ValueError, match="operator changed.*revision 1"):
        restored.create_problem(object())
    accepted = replace(restored, discretization_revision=4)
    assert accepted.to_dict()["discretization_revision"] == 4


@pytest.mark.parametrize("theta,zeta", [(0, 64), (-1, 64), (64, 0), (2**31, 1), (65536, 65536), (True, 64), (64.0, 64)])
def test_invalid_angular_configuration(theta, zeta):
    with pytest.raises((ValueError, TypeError)):
        gliss.StabilityConfiguration(modes=((1, 0),), angular_theta=theta, angular_zeta=zeta)
