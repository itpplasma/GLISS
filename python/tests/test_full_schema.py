import hashlib
import json
import zipfile
from dataclasses import replace

import numpy as np
import pytest

import gliss


def _certified(parity_class, eigenvalues):
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
        floor_count=1,
        lowest_eigenvalue=float(eigenvalues[0]),
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


def _full_class(parity_class):
    eigenvalues = np.array(
        [-6.0 + parity_class, 0.0, 2.0, 4.0, 6.0, 8.0], dtype=np.float64
    )
    arrays = (
        eigenvalues,
        eigenvalues.copy(),
        np.linspace(1.0e-12, 6.0e-12, 6),
        np.linspace(2.0e-12, 7.0e-12, 6),
        np.arange(36, dtype=np.float64).reshape(6, 6) + 100 * parity_class,
    )
    for array in arrays:
        array.setflags(write=False)
    return gliss.FullSpectrumResult(_certified(parity_class, eigenvalues), *arrays)


@pytest.fixture
def full_result():
    return gliss.FullStabilityResult((_full_class(1), _full_class(2)))


@pytest.fixture
def configuration():
    return gliss.StabilityConfiguration(modes=((1, 1), (2, -1)), density_kg_m3=2.0)


def _archive_entries(path):
    with zipfile.ZipFile(path) as archive:
        return {name: archive.read(name) for name in archive.namelist()}


def _rewrite_archive(path, entries, compression=zipfile.ZIP_STORED):
    with zipfile.ZipFile(path, "w", compression=compression) as archive:
        for name, data in entries.items():
            archive.writestr(name, data)


def test_full_result_round_trip_is_exact_read_only_and_deterministic(
    full_result, tmp_path
):
    first = tmp_path / "full-result.gliss"
    second = tmp_path / "full-result-copy.gliss"
    full_result.write(first)
    loaded = gliss.FullStabilityResult.read(first)
    loaded.write(second)

    assert first.read_bytes() == second.read_bytes()
    with zipfile.ZipFile(first) as archive:
        metadata = json.loads(archive.read("metadata.json"))
        assert metadata["schema"] == "gliss.stability.full-result"
        assert metadata["schema_version"] == 1
        assert metadata["storage"]["array_order"] == "eigenpair-component"
        assert all(
            item.compress_type == zipfile.ZIP_STORED for item in archive.infolist()
        )
    for expected, actual in zip(full_result.classes, loaded.classes):
        assert (
            actual.certified_lowest.lowest_eigenvalue
            == expected.certified_lowest.lowest_eigenvalue
        )
        assert actual.certified_lowest.modes == expected.certified_lowest.modes
        np.testing.assert_array_equal(
            actual.certified_lowest.eigenvector,
            expected.certified_lowest.eigenvector,
        )
        for name in (
            "eigenvalues",
            "rayleigh_quotients",
            "residuals",
            "resolutions",
            "eigenvectors",
        ):
            array = getattr(actual, name)
            np.testing.assert_array_equal(array, getattr(expected, name))
            assert array.dtype == np.float64
            assert not array.flags.writeable


@pytest.mark.parametrize(
    ("mutate", "match"),
    [
        (lambda entries: entries.pop("class-2-residuals.npy"), "missing entry"),
        (lambda entries: entries.__setitem__("extra.npy", b"x"), "unknown entry"),
        (
            lambda entries: entries.__setitem__("metadata.json", b'{"schema":'),
            "invalid or truncated JSON",
        ),
        (
            lambda entries: entries.__setitem__(
                "class-1-eigenvalues.npy", entries["class-1-eigenvalues.npy"][:-1]
            ),
            "invalid or truncated NumPy array",
        ),
    ],
)
def test_full_result_rejects_malformed_containers(full_result, tmp_path, mutate, match):
    path = tmp_path / "full-result.gliss"
    full_result.write(path)
    entries = _archive_entries(path)
    mutate(entries)
    _rewrite_archive(path, entries)

    with pytest.raises(ValueError, match=match):
        gliss.FullStabilityResult.read(path)


def test_full_result_rejects_wrong_array_dtype_and_shape(full_result, tmp_path):
    path = tmp_path / "full-result.gliss"
    full_result.write(path)
    entries = _archive_entries(path)
    replacement = tmp_path / "wrong.npy"
    np.save(replacement, np.arange(5, dtype=np.float32), allow_pickle=False)
    entries["class-1-eigenvalues.npy"] = replacement.read_bytes()
    _rewrite_archive(path, entries)

    with pytest.raises(ValueError, match="class-1-eigenvalues.*dtype|shape"):
        gliss.FullStabilityResult.read(path)


def test_full_result_writer_rejects_inconsistent_data_without_replacing_destination(
    full_result, tmp_path
):
    path = tmp_path / "full-result.gliss"
    path.write_bytes(b"existing")
    invalid_class = replace(
        full_result.classes[0], eigenvalues=np.array([0.0, -2.0, 3.0])
    )
    invalid = gliss.FullStabilityResult((invalid_class, full_result.classes[1]))

    with pytest.raises(ValueError, match="shape|sorted|negative counts"):
        invalid.write(path)
    assert path.read_bytes() == b"existing"


def test_full_result_rejects_truncated_archive(full_result, tmp_path):
    path = tmp_path / "full-result.gliss"
    full_result.write(path)
    path.write_bytes(path.read_bytes()[:100])

    with pytest.raises(
        ValueError, match="invalid or truncated full-spectrum container"
    ):
        gliss.FullStabilityResult.read(path)


def test_full_result_rejects_boolean_schema_version(full_result, tmp_path):
    path = tmp_path / "full-result.gliss"
    full_result.write(path)
    entries = _archive_entries(path)
    metadata = json.loads(entries["metadata.json"])
    metadata["schema_version"] = True
    entries["metadata.json"] = json.dumps(metadata).encode()
    _rewrite_archive(path, entries)

    with pytest.raises(ValueError, match="schema_version.*expected 1"):
        gliss.FullStabilityResult.read(path)


def test_full_result_rejects_oversized_declared_shape(full_result, tmp_path):
    path = tmp_path / "full-result.gliss"
    full_result.write(path)
    entries = _archive_entries(path)
    metadata = json.loads(entries["metadata.json"])
    item = metadata["certified_result"]["classes"][0]
    item["normal_unknowns"] = 2**63
    item["eta_unknowns"] = 0
    item["mu_unknowns"] = 0
    item["has_eigenvector"] = False
    item["eigenvector"] = []
    entries["metadata.json"] = json.dumps(metadata).encode()
    _rewrite_archive(path, entries)

    with pytest.raises(ValueError, match="unknown count.*NumPy.*limit"):
        gliss.FullStabilityResult.read(path)


@pytest.mark.parametrize("compression", [zipfile.ZIP_DEFLATED, zipfile.ZIP_BZIP2])
def test_full_result_rejects_compressed_entries(full_result, tmp_path, compression):
    path = tmp_path / "full-result.gliss"
    full_result.write(path)
    _rewrite_archive(path, _archive_entries(path), compression=compression)

    with pytest.raises(ValueError, match="must be stored"):
        gliss.FullStabilityResult.read(path)


def test_full_result_rejects_duplicate_entries(full_result, tmp_path):
    path = tmp_path / "full-result.gliss"
    full_result.write(path)
    entries = _archive_entries(path)
    with pytest.warns(UserWarning, match="Duplicate name"):
        with zipfile.ZipFile(path, "w") as archive:
            for name, data in entries.items():
                archive.writestr(name, data)
            archive.writestr("class-1-residuals.npy", entries["class-1-residuals.npy"])

    with pytest.raises(ValueError, match="duplicate entries"):
        gliss.FullStabilityResult.read(path)


def test_full_result_reports_array_allocation_failure(
    full_result, tmp_path, monkeypatch
):
    path = tmp_path / "full-result.gliss"
    full_result.write(path)

    def fail_allocation(*args, **kwargs):
        raise MemoryError

    monkeypatch.setattr(np.lib.format, "read_array", fail_allocation)
    with pytest.raises(gliss.GlissAllocationError, match="class-1-eigenvalues"):
        gliss.FullStabilityResult.read(path)


def test_full_result_write_failure_preserves_destination(
    full_result, tmp_path, monkeypatch
):
    path = tmp_path / "full-result.gliss"
    path.write_bytes(b"existing")

    def fail_write(*args, **kwargs):
        raise OSError("disk full")

    monkeypatch.setattr(np.lib.format, "write_array", fail_write)
    with pytest.raises(OSError, match="disk full"):
        full_result.write(path)
    assert path.read_bytes() == b"existing"
    assert list(tmp_path.iterdir()) == [path]


def test_full_result_reports_conversion_allocation_failure(
    full_result, tmp_path, monkeypatch
):
    path = tmp_path / "full-result.gliss"
    path.write_bytes(b"existing")

    def fail_allocation(*args, **kwargs):
        raise MemoryError

    monkeypatch.setattr(np, "ascontiguousarray", fail_allocation)
    with pytest.raises(gliss.GlissAllocationError, match="class 1 eigenvalues"):
        full_result.write(path)
    assert path.read_bytes() == b"existing"


def test_full_result_writer_enforces_reader_metadata_limit(
    full_result, tmp_path, monkeypatch
):
    path = tmp_path / "full-result.gliss"
    monkeypatch.setattr("gliss.full_schema._MAX_METADATA_BYTES", 1)

    with pytest.raises(ValueError, match="metadata exceeds 1 byte"):
        full_result.write(path)
    assert not path.exists()


def test_full_and_certified_formats_are_not_implicitly_interchangeable(
    full_result, tmp_path
):
    full_path = tmp_path / "full-result.gliss"
    certified_path = tmp_path / "certified-result.json"
    full_result.write(full_path)
    certified = gliss.StabilityResult(
        tuple(item.certified_lowest for item in full_result.classes)
    )
    certified.write(certified_path)

    with pytest.raises(ValueError, match="valid UTF-8"):
        gliss.StabilityResult.read(full_path)
    with pytest.raises(ValueError, match="full-spectrum container"):
        gliss.FullStabilityResult.read(certified_path)


def test_full_run_manifest_is_self_contained_and_deterministic(
    full_result, configuration, tmp_path, monkeypatch
):
    equilibrium = tmp_path / "private" / "equilibrium.nc"
    equilibrium.parent.mkdir()
    equilibrium.write_bytes(b"equilibrium fixture")
    first = tmp_path / "full-run.gliss"
    second = tmp_path / "full-run-copy.gliss"
    monkeypatch.setattr("gliss.schema._native_version", lambda: "0.0.1")
    monkeypatch.setattr("gliss.schema._equilibrium_schema_version", lambda path: 0)

    manifest = gliss.write_full_run_manifest(
        first, equilibrium, configuration, full_result
    )
    loaded = gliss.FullRunManifest.read(first)
    loaded.write(second)

    assert loaded == manifest
    assert first.read_bytes() == second.read_bytes()
    assert loaded.equilibrium_filename == "equilibrium.nc"
    assert (
        loaded.equilibrium_sha256
        == hashlib.sha256(equilibrium.read_bytes()).hexdigest()
    )
    loaded.verify_equilibrium(equilibrium)
    assert str(equilibrium.parent).encode() not in first.read_bytes()
    with zipfile.ZipFile(first) as archive:
        metadata = json.loads(archive.read("metadata.json"))
    assert metadata["schema"] == "gliss.stability.full-run"
    assert metadata["configuration"] == configuration.to_dict()


def test_full_run_manifest_rejects_configuration_mismatch(
    full_result, configuration, tmp_path, monkeypatch
):
    equilibrium = tmp_path / "equilibrium.nc"
    equilibrium.touch()
    mismatch = replace(configuration, modes=((3, 2),))
    monkeypatch.setattr("gliss.schema._native_version", lambda: "0.0.1")
    monkeypatch.setattr("gliss.schema._equilibrium_schema_version", lambda path: 0)

    with pytest.raises(ValueError, match="result modes.*configuration modes"):
        gliss.write_full_run_manifest(
            tmp_path / "full-run.gliss", equilibrium, mismatch, full_result
        )
    assert not (tmp_path / "full-run.gliss").exists()


def test_full_run_manifest_rejects_wrong_result_type(configuration):
    with pytest.raises(TypeError, match="gliss.FullStabilityResult"):
        gliss.FullRunManifest(
            equilibrium_filename="equilibrium.nc",
            equilibrium_size_bytes=1,
            equilibrium_sha256="0" * 64,
            equilibrium_schema_version=1,
            configuration=configuration,
            result=object(),
            gliss_python_version="0.0.1",
            gliss_native_version="0.0.1",
            gliss_abi_version=1,
            numpy_version=np.__version__,
            python_version="3.9.0",
        )
