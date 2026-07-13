from pathlib import Path

import numpy as np
import pytest
from scipy.io import netcdf_file

import gliss.vmec as vmec
from gliss._vmec_geometry import ConvertedGeometry, _EVEN_FIELDS, _ODD_FIELDS


class _Transform:
    nfp = 3
    asym = False

    def __init__(self):
        self.verbose = 1
        self.mboz = 0
        self.nboz = 0
        self.read_arguments = None
        self.ran = False

    def read_wout(self, path, flux):
        self.read_arguments = (path, flux)

    def run(self):
        self.ran = True


class _BoozModule:
    Booz_xform = _Transform


def _geometry(m_max=2, n_max=1):
    shape = (5, m_max + 1, 2 * n_max + 1)
    harmonics = {}
    for name in _EVEN_FIELDS | _ODD_FIELDS:
        cosine = np.zeros(shape)
        sine = np.zeros(shape)
        target = cosine if name in _EVEN_FIELDS else sine
        target[:, 0, 0] = 1.0
        harmonics[name] = cosine, sine
    profiles = {
        name: np.linspace(0.1, 0.5, 5)
        for name in ("p", "B_theta_avg", "B_zeta_avg", "Phi", "chi", "iota")
    }
    return ConvertedGeometry(
        (np.arange(5) + 0.5) / 5,
        profiles,
        harmonics,
        {"toroidal_flux": 0.0},
    )


def test_convert_vmec_validates_paths_before_optional_dependencies(tmp_path):
    with pytest.raises(FileNotFoundError, match="input_path does not exist"):
        vmec.convert_vmec(tmp_path / "missing.nc", tmp_path / "out.nc")

    source = tmp_path / "wout.nc"
    source.write_bytes(b"not netcdf")
    destination = tmp_path / "out.nc"
    destination.write_bytes(b"preserve")
    with pytest.raises(FileExistsError, match="already exists"):
        vmec.convert_vmec(source, destination)
    assert destination.read_bytes() == b"preserve"


def test_convert_vmec_rejects_overwriting_its_input(tmp_path):
    source = tmp_path / "wout.nc"
    source.write_bytes(b"preserve")
    with pytest.raises(ValueError, match="must name different files"):
        vmec.convert_vmec(source, source, overwrite=True)
    assert source.read_bytes() == b"preserve"


@pytest.mark.parametrize("name", ["poloidal_max", "toroidal_max", "transform_factor"])
def test_convert_vmec_rejects_boolean_integer_options(tmp_path, name):
    source = tmp_path / "wout.nc"
    source.write_bytes(b"input")
    with pytest.raises(TypeError, match=f"{name} must be an integer"):
        vmec.convert_vmec(source, tmp_path / "out.nc", **{name: True})


def test_convert_vmec_is_atomic_and_writes_reader_schema(tmp_path, monkeypatch):
    source = tmp_path / "wout.nc"
    source.write_bytes(b"input")
    destination = tmp_path / "converted.nc"
    transform = _Transform()
    monkeypatch.setattr(
        vmec,
        "_dependencies",
        lambda: (type("Module", (), {"Booz_xform": lambda: transform}), netcdf_file),
    )
    monkeypatch.setattr(vmec, "_metadata", lambda path, reader: 0.02)
    monkeypatch.setattr(vmec, "convert_geometry", lambda *args: _geometry())

    result = vmec.convert_vmec(source, destination, poloidal_max=2, toroidal_max=1)

    assert result == destination
    assert transform.verbose == 0
    assert transform.read_arguments == (str(source), True)
    assert transform.mboz == 12
    assert transform.nboz == 8
    assert transform.ran
    assert not list(tmp_path.glob(".converted.nc.*.tmp"))
    with netcdf_file(destination, "r", mmap=False) as file:
        assert file.gliss_schema == b"gvec-cas3d-export"
        assert file.gliss_schema_version == b"1"
        assert file.stellarator_symmetry == b"True"
        assert file.creator == b"gliss.convert_vmec"
        assert float(file.conversion_residual_toroidal_flux) == 0.0
        assert int(file.variables["N_FP"].data) == 3
        assert file.variables["Jac_mnc"].data.shape == (5, 3, 3)
        assert "Jac_mns" not in file.variables
        assert file.variables["g_st_mns"].data.shape == (5, 3, 3)


def test_convert_vmec_removes_partial_output_after_write_failure(tmp_path, monkeypatch):
    source = tmp_path / "wout.nc"
    source.write_bytes(b"input")
    destination = tmp_path / "converted.nc"
    monkeypatch.setattr(vmec, "_dependencies", lambda: (_BoozModule, netcdf_file))
    monkeypatch.setattr(vmec, "_metadata", lambda path, reader: 0.02)
    monkeypatch.setattr(vmec, "convert_geometry", lambda *args: _geometry())

    def fail(*args):
        Path(args[0]).write_bytes(b"partial")
        raise OSError("disk failure")

    monkeypatch.setattr(vmec, "_write", fail)
    with pytest.raises(OSError, match="disk failure"):
        vmec.convert_vmec(source, destination)
    assert not destination.exists()
    assert not list(tmp_path.glob(".converted.nc.*.tmp"))
