"""Conventions and provenance of the CAS3D experiment runner."""

import importlib.util
import json
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

_PATH = Path(__file__).resolve().parents[2] / "benchmarks/cas3d_envelope/run.py"
_SPEC = importlib.util.spec_from_file_location("cas3d_benchmark_runner", _PATH)
runner = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(runner)


def configuration(*options):
    parser = runner.argument_parser()
    arguments = parser.parse_args(["input.nc", "result.json", *map(str, options)])
    return arguments, runner.configure(arguments, parser)


def test_synthetic_order_and_finite_table():
    assert runner.envelope_table(5) == [(0, 0), (0, 1), (1, 0), (0, 2), (1, -1)]
    assert len(runner.envelope_table(70)) == 70
    with pytest.raises(ValueError, match="between"):
        runner.envelope_table(1000)
    _, (settings, table, runs) = configuration()
    assert settings["scope"] == "synthetic_sparse_scaling_control"
    assert settings["coefficient_angular_resolution"] == [64, 64]
    assert settings["radial_cells"] is None
    assert table == runner.envelope_table(70)
    assert [run["envelope_entries"] for run in runs] == [3, 10, 20, 35, 70]


def test_ordered_csv_preserves_rows_and_declares_sign(tmp_path):
    source = tmp_path / "table.csv"
    source.write_text("# published order\ncell,m,n\n0,0,0\n1,2,-3\n2,1,2\n")
    assert runner.read_envelope_csv(source, "gliss") == [(0, 0), (2, -3), (1, 2)]
    assert runner.read_envelope_csv(source, "opposite") == [(0, 0), (2, 3), (1, -2)]
    with pytest.raises(SystemExit):
        configuration("--envelope-csv", source)


@pytest.mark.parametrize("body,message", [
    ("x,y\n0,0\n", "header"),
    ("m,n,n\n0,0,0\n", "repeat"),
    ("m,n\n1,0\n", "start"),
    ("m,n\n0,0\n0,0\n", "duplicates"),
    ("m,n\n0,0\n1.5,2\n", "integer"),
    ("m,n\n0,0\n-1,2\n", "invalid"),
])
def test_invalid_csv_is_rejected(tmp_path, body, message):
    source = tmp_path / "table.csv"
    source.write_text(body)
    with pytest.raises(ValueError, match=message):
        runner.read_envelope_csv(source, "gliss")


@pytest.mark.parametrize("convention,count", [("caption_l", 2), ("terminal_l_plus_one", 3)])
def test_literature_prefixes_are_explicit(tmp_path, convention, count):
    source = tmp_path / "table.csv"
    source.write_text("m,n\n0,0\n1,1\n2,1\n")
    _, (settings, table, runs) = configuration(
        "--envelope-csv", source, "--envelope-n-sign", "opposite",
        "--literature-l", 2, "--prefix-convention", convention)
    assert settings["prefix_convention"] == convention
    assert runs == [{"literature_l": 2, "envelope_entries": count}]
    assert table[:count] == [(0, 0), (1, -1), (2, -1)][:count]
    with pytest.raises(SystemExit):
        configuration("--envelope-csv", source, "--envelope-n-sign", "gliss",
                      "--literature-l", 2)


def test_split_grids_parity_and_native_radial_options_reach_child(tmp_path):
    source = tmp_path / "table.csv"
    source.write_text("m,n\n0,0\n1,1\n")
    arguments, (settings, _, _) = configuration(
        "--envelope-csv", source, "--envelope-n-sign", "opposite", "--entries", 2,
        "--base", 3, 2, "--parity", 2, "--angular", 96, 64,
        "--coefficient-angular", 36, 24, "--reference-length", 10,
        "--radial-cells", 12, "--radial-quadrature", "cas3d_midpoint")
    child = runner.argument_parser().parse_args(runner.child_command(arguments, 2)[2:])
    child_settings, child_table, child_runs = runner.configure(child, runner.argument_parser())
    assert child_settings == settings
    assert child_table == [(0, 0), (1, -1)]
    assert child_runs == [{"envelope_entries": 2}]


def test_rss_platform_units():
    assert runner.peak_resident_mib(104857600, "Darwin") == 100.0
    assert runner.peak_resident_mib(102400, "Linux") == 100.0
    with pytest.raises(ValueError, match="not declared"):
        runner.peak_resident_mib(100, "Unknown")


def test_coefficient_runner_rejects_unsupported_coupled_parity():
    with pytest.raises(SystemExit):
        configuration("--parity", 0)


def test_solver_receives_declared_options_and_preserves_certificates(tmp_path, monkeypatch):
    library = tmp_path / "libgliss_c.fake"
    library.write_bytes(b"identity")
    calls = []
    result = SimpleNamespace(
        field_periods=5, radial_surfaces=13, quotient_rank=24, labeled_nullity=0,
        peak_block_width=4, negative_count=1, lowest_eigenvalue=-0.2, certificate=1e-9,
        eigenpair_residual=2e-10, force_balance_residual=0.01, inertia_zero_floor=3e-12,
        coordinate_handedness="left-handed", has_eigenpair=True, labeled_sideband_count=3)

    class Equilibrium:
        _library = SimpleNamespace(_name=str(library))

        def __init__(self, path):
            self.path = path

        def __enter__(self):
            return self

        def __exit__(self, *unused):
            return False

    def solve(equilibrium, **kwargs):
        calls.append((equilibrium.path, kwargs))
        return result

    monkeypatch.setitem(sys.modules, "gliss", SimpleNamespace(
        Equilibrium=Equilibrium, solve_cas3d_phase_envelope=solve))
    _, (settings, _, _) = configuration(
        "--entries", 2, "--parity", 2, "--angular", 96, 64,
        "--coefficient-angular", 36, 24, "--reference-length", 10,
        "--radial-cells", 12, "--radial-quadrature", "cas3d_midpoint")
    solved = runner.solve("input.nc", [(0, 0), (1, -2)], settings)
    assert calls[0][1] == {
        "base_mode": (1, 1), "envelope_modes": [(0, 0), (1, -2)], "parity_class": 2,
        "degree": 1, "angular_theta": 96, "angular_zeta": 64,
        "normalization": "cas3d2mn_coefficient", "coefficient_angular_resolution": (36, 24),
        "reference_length": 10.0, "radial_cells": 12, "radial_quadrature": "cas3d_midpoint"}
    assert solved["certificate"] == 1e-9
    assert solved["eigenpair_residual"] == 2e-10
    assert solved["inertia_zero_floor"] == 3e-12
    assert solved["library"]["sha256"] == runner.sha256_file(library)


def test_manifest_is_frozen_before_first_solve_and_retained_on_failure(tmp_path, monkeypatch):
    export = tmp_path / "input.nc"
    export.write_bytes(b"equilibrium identity")
    library = tmp_path / "libgliss_c.fake"
    library.write_bytes(b"library identity")
    output = tmp_path / "manifest.json"
    monkeypatch.setitem(sys.modules, "gliss", SimpleNamespace(
        _load_library=lambda: SimpleNamespace(_name=str(library))))
    original_run = subprocess.run

    def run(command, **kwargs):
        if "--single" not in command:
            return original_run(command, **kwargs)
        frozen = json.loads(output.read_text())
        assert frozen["status"] == "running"
        assert frozen["runs"] == []
        assert frozen["input"]["sha256"] == runner.sha256_file(export)
        assert frozen["configuration"]["assumptions"] == ["Rref is provisional"]
        assert frozen["source"]["runner_sha256"] == runner.sha256_file(_PATH)
        return subprocess.CompletedProcess(command, 9, "", "native solve failed")

    monkeypatch.setattr(runner.subprocess, "run", run)
    with pytest.raises(SystemExit) as error:
        runner.main([str(export), str(output), "--entries", "2",
                     "--assumption", "Rref is provisional"])
    assert error.value.code == 1
    frozen = json.loads(output.read_text())
    assert frozen["status"] == "failed"
    assert frozen["failure"]["returncode"] == 9
    assert frozen["failure"]["stderr"] == "native solve failed"


@pytest.mark.parametrize("changed", [None, "input", "csv", "library"])
def test_completed_result_requires_frozen_input_table_and_library(tmp_path, monkeypatch, changed):
    export = tmp_path / "input.nc"
    export.write_bytes(b"equilibrium identity")
    library = tmp_path / "libgliss_c.fake"
    library.write_bytes(b"library identity")
    csv = tmp_path / "table.csv"
    csv.write_text("m,n\n0,0\n1,1\n")
    output = tmp_path / "manifest.json"
    identity = {"path": str(library), "sha256": runner.sha256_file(library)}
    monkeypatch.setitem(sys.modules, "gliss", SimpleNamespace(
        _load_library=lambda: SimpleNamespace(_name=str(library))))
    original_run = subprocess.run

    def run(command, **kwargs):
        if "--single" not in command:
            return original_run(command, **kwargs)
        if changed is not None:
            {"input": export, "csv": csv, "library": library}[changed].write_bytes(b"changed")
        return subprocess.CompletedProcess(command, 0, json.dumps({"library": identity}), "")

    monkeypatch.setattr(runner.subprocess, "run", run)
    options = [str(export), str(output), "--envelope-csv", str(csv),
               "--envelope-n-sign", "opposite", "--literature-l", "2",
               "--prefix-convention", "caption_l"]
    if changed is None:
        runner.main(options)
        manifest = json.loads(output.read_text())
        assert manifest["status"] == "complete"
        assert manifest["runs"][0]["literature_l"] == 2
        assert manifest["runs"][0]["envelope_entries"] == 2
    else:
        with pytest.raises(SystemExit):
            runner.main(options)
        manifest = json.loads(output.read_text())
        assert manifest["status"] == "failed"
        assert "changed after assumptions were frozen" in manifest["failure"]
        assert manifest["runs"] == []
