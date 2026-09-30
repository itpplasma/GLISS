from pathlib import Path

import numpy as np
import pytest
from scipy.io import netcdf_file

import gliss._vmec_geometry as _vmec_geometry
import gliss.vmec as vmec
from gliss._vmec_geometry import (
    ConvertedGeometry,
    _EVEN_FIELDS,
    _ODD_FIELDS,
    _position_frame,
)


class _Transform:
    nfp = 3
    asym = False
    ns_in = 15

    def __init__(self):
        self.verbose = 1
        self.mboz = 0
        self.nboz = 0
        self.read_arguments = None
        self.ran = False
        self.compute_surfs = np.arange(self.ns_in, dtype=np.int32)

    def read_wout(self, path, flux):
        self.read_arguments = (path, flux)

    def run(self):
        self.ran = True


class _BoozModule:
    Booz_xform = _Transform


def _geometry(m_max=2, n_max=1, force_balance=0.0):
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
        {"toroidal_flux": 0.0, "force_balance": force_balance},
    )


@pytest.mark.parametrize(
    ("radius", "boozer_shift", "boozer_zeta"),
    [(6.2, 0.0, 0.37), (5.9, -0.23, 1.41), (6.6, 0.31, -0.72)],
)
def test_vmec_position_frame_reconstructs_cylindrical_position(
    radius, boozer_shift, boozer_zeta
):
    xhat, yhat = _position_frame(np.asarray(radius), np.asarray(boozer_shift))
    nfp = 5
    zeta_period = -nfp * boozer_zeta / (2.0 * np.pi)
    winding = -1
    rotation = 2.0 * np.pi * winding * zeta_period / nfp
    physical_x = np.cos(rotation) * xhat - np.sin(rotation) * yhat
    physical_y = np.sin(rotation) * xhat + np.cos(rotation) * yhat
    expected_angle = boozer_zeta - boozer_shift
    np.testing.assert_allclose(
        [physical_x, physical_y],
        radius * np.asarray([np.cos(expected_angle), np.sin(expected_angle)]),
        rtol=0.0,
        atol=2.0e-15,
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


def test_convert_vmec_wraps_truncated_netcdf_error(tmp_path, monkeypatch):
    source = tmp_path / "wout.nc"
    source.write_bytes(b"truncated")

    def broken_reader(*args, **kwargs):
        raise ValueError("cannot reshape truncated variable")

    monkeypatch.setattr(vmec, "_dependencies", lambda: (_BoozModule, broken_reader))
    with pytest.raises(ValueError, match="cannot read standard VMEC NetCDF file"):
        vmec.convert_vmec(source, tmp_path / "out.nc")


@pytest.mark.parametrize(
    "name", ["poloidal_max", "toroidal_max", "transform_factor", "radial_surfaces"]
)
def test_convert_vmec_rejects_boolean_integer_options(tmp_path, name):
    source = tmp_path / "wout.nc"
    source.write_bytes(b"input")
    with pytest.raises(TypeError, match=f"{name} must be an integer"):
        vmec.convert_vmec(source, tmp_path / "out.nc", **{name: True})


@pytest.mark.parametrize("value", [None, True, 1])
def test_convert_vmec_rejects_nonstring_force_balance_policy(tmp_path, value):
    source = tmp_path / "wout.nc"
    source.write_bytes(b"input")
    with pytest.raises(TypeError, match="force_balance_policy must be a string"):
        vmec.convert_vmec(
            source, tmp_path / "out.nc", force_balance_policy=value
        )


def test_convert_vmec_rejects_unknown_force_balance_policy(tmp_path):
    source = tmp_path / "wout.nc"
    source.write_bytes(b"input")
    with pytest.raises(ValueError, match="must be 'error' or 'warn'"):
        vmec.convert_vmec(
            source, tmp_path / "out.nc", force_balance_policy="ignore"
        )


@pytest.mark.parametrize("policy", ["error", "warn"])
def test_convert_vmec_force_balance_policy_is_explicit(
    tmp_path, monkeypatch, policy
):
    source = tmp_path / "wout.nc"
    source.write_bytes(b"input")
    destination = tmp_path / "converted.nc"
    monkeypatch.setattr(vmec, "_dependencies", lambda: (_BoozModule, netcdf_file))
    monkeypatch.setattr(vmec, "_metadata", lambda path, reader: 0.02)
    monkeypatch.setattr(
        vmec, "convert_geometry", lambda *args: _geometry(force_balance=0.2)
    )
    if policy == "error":
        with pytest.raises(ValueError, match="force_balance"):
            vmec.convert_vmec(source, destination)
        assert not destination.exists()
    else:
        with pytest.warns(RuntimeWarning, match="force_balance"):
            vmec.convert_vmec(
                source,
                destination,
                poloidal_max=2,
                toroidal_max=1,
                force_balance_policy="warn",
            )
        with netcdf_file(destination, "r", mmap=False) as file:
            assert float(file.conversion_residual_force_balance) == 0.2


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

    result = vmec.convert_vmec(
        source,
        destination,
        poloidal_max=2,
        toroidal_max=1,
        radial_surfaces=5,
    )

    assert result == destination
    assert transform.verbose == 0
    assert transform.read_arguments == (str(source), True)
    assert transform.mboz == 12
    assert transform.nboz == 8
    np.testing.assert_array_equal(transform.compute_surfs, [1, 4, 7, 10, 13])
    assert transform.ran
    assert not list(tmp_path.glob(".converted.nc.*.tmp"))
    with netcdf_file(destination, "r", mmap=False) as file:
        assert file.gliss_schema == b"gvec-cas3d-export"
        assert file.gliss_schema_version == b"1"
        assert file.stellarator_symmetry == b"True"
        assert file.position_frame == b"xhat,yhat rotated by winding*zeta_B"
        assert file.creator == b"gliss.convert_vmec"
        assert int(file.vmec_half_grid_surfaces) == 15
        assert int(file.booz_xform_surfaces) == 5
        assert float(file.conversion_residual_toroidal_flux) == 0.0
        assert int(file.variables["N_FP"].data) == 3
        assert int(file.variables["winding"].data) == -1
        assert file.variables["Jac_mnc"].data.shape == (5, 3, 3)
        assert "Jac_mns" not in file.variables
        assert file.variables["g_st_mns"].data.shape == (5, 3, 3)


@pytest.mark.parametrize(("available", "radial_surfaces"), [(15, 6), (30, 5)])
def test_convert_vmec_rejects_noncentered_radial_subsampling(
    tmp_path, monkeypatch, available, radial_surfaces
):
    source = tmp_path / "wout.nc"
    source.write_bytes(b"input")
    destination = tmp_path / "converted.nc"
    transform = _Transform()
    transform.ns_in = available
    monkeypatch.setattr(
        vmec,
        "_dependencies",
        lambda: (type("Module", (), {"Booz_xform": lambda: transform}), netcdf_file),
    )
    monkeypatch.setattr(vmec, "_metadata", lambda path, reader: 0.02)

    with pytest.raises(ValueError, match="radial_surfaces"):
        vmec.convert_vmec(source, destination, radial_surfaces=radial_surfaces)

    assert not transform.ran
    assert not destination.exists()


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


def test_convert_vmec_rejects_inaccurate_position_truncation(tmp_path, monkeypatch):
    source = tmp_path / "wout.nc"
    source.write_bytes(b"not read")
    geometry = _geometry()
    geometry.residuals["truncated_jacobian"] = 0.2
    monkeypatch.setattr(vmec, "_dependencies", lambda: (_BoozModule, netcdf_file))
    monkeypatch.setattr(vmec, "_metadata", lambda path, reader: 0.02)
    monkeypatch.setattr(vmec, "convert_geometry", lambda *args: geometry)
    with pytest.raises(ValueError, match="increase poloidal_max"):
        vmec.convert_vmec(source, tmp_path / "out.nc", poloidal_max=2, toroidal_max=1)
    assert not (tmp_path / "out.nc").exists()
    with pytest.raises(ValueError, match="truncation_tolerance"):
        vmec.convert_vmec(source, tmp_path / "out.nc", truncation_tolerance=0.0)


def test_truncated_jacobian_residual_detects_missing_harmonics():
    # A shaped torus whose radius needs m=2: truncating the export to m<=1
    # must be visible in the Jacobian that GLISS reconstructs.
    s = (np.arange(8) + 0.5) / 8
    theta = np.linspace(0.0, 2.0 * np.pi, 32, endpoint=False)
    zeta = np.linspace(0.0, 2.0 * np.pi, 8, endpoint=False)
    rho = np.sqrt(s)[:, None, None]
    t = theta[None, :, None] + 0.0 * zeta[None, None, :]
    radius = 3.0 + 0.5 * rho * np.cos(t) + 0.2 * rho**2 * np.cos(2 * t)
    radius_t = -0.5 * rho * np.sin(t) - 0.4 * rho**2 * np.sin(2 * t)
    height = 0.6 * rho * np.sin(t)
    height_t = 0.6 * rho * np.cos(t)
    zero = np.zeros_like(radius)
    exact = (
        np.stack((radius, zero, height)),
        np.stack((radius_t, zero, height_t)),
        np.stack((zero, zero, zero)),
    )
    residuals = []
    for m_max in (1, 2):
        harmonics = {
            name: _vmec_geometry._project(values, theta, zeta, 1, m_max, 0)
            for name, values in (("xhat", radius), ("yhat", zero), ("zhat", height))
        }
        residuals.append(
            _vmec_geometry._truncated_jacobian_residual(
                harmonics, exact, s, theta, zeta, 1
            )
        )
    assert residuals[0] > 0.1
    assert residuals[1] < 1e-12


def test_truncated_metric_positivity_detects_indefinite_metric():
    # g_tt = g_zz = 1 and g_tz = a cos(theta): det/(g_tt g_zz) = 1 - a^2 at
    # theta = 0, positive only for |a| < 1.
    theta = np.linspace(0.0, 2.0 * np.pi, 16, endpoint=False)
    zeta = np.linspace(0.0, 2.0 * np.pi, 4, endpoint=False)

    def pair(constant, first):
        cosine = np.zeros((1, 2, 1))
        cosine[0, 0, 0] = constant
        cosine[0, 1, 0] = first
        return cosine, np.zeros_like(cosine)

    for amplitude, expected in ((0.5, 0.75), (1.5, -1.25)):
        harmonics = {
            "g_tt": pair(1.0, 0.0),
            "g_zz": pair(1.0, 0.0),
            "g_tz": pair(0.0, amplitude),
        }
        positivity = _vmec_geometry._truncated_metric_positivity(
            harmonics, theta, zeta, 1
        )
        assert positivity == pytest.approx(expected)


def test_force_balance_gate_uses_the_flux_surface_average():
    # Uniform fields with mu0 p' J + F G' + P I' = 0 close the averaged
    # radial force balance exactly; an angular metric B_s = g_st B^theta
    # (g_st = a sin theta) breaks only the pointwise closure, which GLISS
    # does not use because it solves the magnetic differential equation.
    ns = 12
    s = np.linspace(0.06, 0.88, ns)
    theta = np.linspace(0.0, 2.0 * np.pi, 16, endpoint=False)
    zeta = np.linspace(0.0, 2.0 * np.pi, 4, endpoint=False)

    def pair(cosine0=0.0, sine1=0.0):
        cosine = np.zeros((ns, 2, 1))
        sine = np.zeros((ns, 2, 1))
        cosine[:, 0, 0] = cosine0
        sine[:, 1, 0] = sine1
        return cosine, sine

    phi_slope, chi_slope, pressure_slope = 0.8, 0.3, -2.0e4
    jacobian, g_slope = 1.5, 0.7
    # mu0 p' J + F G' + P I' = 0 with F = -Phi', P = -chi', I' = 0.
    g_zeta = 4.0e-7 * np.pi * pressure_slope * jacobian / phi_slope
    profiles = {
        "p": pressure_slope * s,
        "Phi": phi_slope * s,
        "chi": chi_slope * s,
        "B_theta_avg": np.zeros(ns),
        "B_zeta_avg": g_zeta * s,
        "iota": np.full(ns, chi_slope / phi_slope),
    }
    for amplitude, expected_pointwise in ((0.0, False), (0.2, True)):
        harmonics = {
            "Jac": pair(cosine0=jacobian),
            "B_contra_t": pair(cosine0=g_slope),
            "B_contra_z": pair(cosine0=g_slope),
            "g_st": pair(sine1=amplitude),
            "g_sz": pair(),
        }
        average, pointwise = _vmec_geometry._force_balance_residual(
            harmonics, profiles, s, theta, zeta, 1
        )
        assert average < 1.0e-12
        assert (pointwise > 1.0e-2) == expected_pointwise
    profiles["B_zeta_avg"] = 2.0 * g_zeta * s
    average, _ = _vmec_geometry._force_balance_residual(
        harmonics, profiles, s, theta, zeta, 1
    )
    assert average > 1.0e-1
