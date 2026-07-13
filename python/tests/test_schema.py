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
        radial_quadrature="midpoint",
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
    )


@pytest.fixture
def configuration():
    return gliss.StabilityConfiguration(
        modes=[(1, 1), (2, -1)],
        adiabatic_index=5.0 / 3.0,
        density_kg_m3=2.0,
        zero_floor=1.0,
        radial_quadrature="midpoint",
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
    assert document["schema_version"] == 1
    assert document["boundary_condition"] == "fixed"


@pytest.mark.parametrize(
    ("change", "match"),
    [
        (lambda value: value.update(extra=1), "unknown field 'extra'"),
        (lambda value: value.pop("modes"), "missing field 'modes'"),
        (lambda value: value.update(schema_version=2), "schema_version.*expected 1"),
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

    manifest = gliss.write_run_manifest(path, equilibrium, configuration, result)
    loaded = gliss.RunManifest.read(path)
    text = path.read_text(encoding="utf-8")

    assert loaded == manifest
    assert loaded.equilibrium_filename == "equilibrium.nc"
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
    mismatch = replace(configuration, modes=((3, 2),))
    with pytest.raises(ValueError, match="result modes.*configuration modes"):
        gliss.write_run_manifest(tmp_path / "run.json", equilibrium, mismatch, result)


def test_manifest_rejects_changed_equilibrium(configuration, result, tmp_path):
    equilibrium = tmp_path / "equilibrium.nc"
    equilibrium.write_bytes(b"before")
    path = tmp_path / "run.json"
    document = {
        "schema": "gliss.stability.run",
        "schema_version": 1,
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
            configuration=configuration,
            result=result,
            gliss_python_version="0.0.1",
            gliss_native_version="0.0.1",
            gliss_abi_version=1,
            numpy_version=np.__version__,
            python_version="3.9.0",
        )
