"""Conversion of standard VMEC ``wout`` files to GLISS equilibria."""

import operator
import os
import tempfile
from pathlib import Path
from typing import Any, Optional, Union

import numpy as np

from ._vmec_geometry import _EVEN_FIELDS, convert_geometry


PathLike = Union[str, os.PathLike]


def _path(value: PathLike, name: str, must_exist: bool) -> Path:
    try:
        text = os.fspath(value)
    except TypeError as error:
        raise TypeError(f"{name} must be a string or path-like object") from error
    if not isinstance(text, str):
        raise TypeError(f"{name} must resolve to a string")
    if "\0" in text:
        raise ValueError(f"{name} contains a null byte")
    path = Path(text)
    if must_exist:
        if not path.exists():
            raise FileNotFoundError(f"{name} does not exist: {path}")
        if not path.is_file():
            raise ValueError(f"{name} is not a regular file: {path}")
    elif not path.parent.is_dir():
        raise FileNotFoundError(f"{name} directory does not exist: {path.parent}")
    return path


def _integer(value: Any, name: str, minimum: int, maximum: int) -> int:
    if isinstance(value, (bool, np.bool_)):
        raise TypeError(f"{name} must be an integer")
    try:
        result = operator.index(value)
    except TypeError as error:
        raise TypeError(f"{name} must be an integer") from error
    if not minimum <= result <= maximum:
        raise ValueError(f"{name} must be between {minimum} and {maximum}")
    return result


def _dependencies():
    try:
        import booz_xform  # type: ignore[import-not-found]
        from scipy.io import netcdf_file  # type: ignore[import-untyped]
    except ImportError as error:
        raise ImportError(
            "VMEC import requires the optional dependencies; install gliss[vmec]"
        ) from error
    return booz_xform, netcdf_file


def _scalar(file: Any, name: str) -> float:
    if name not in file.variables:
        raise ValueError(f"VMEC wout file is missing {name}")
    value = np.asarray(file.variables[name].data)
    if value.size != 1 or not np.isfinite(value.item()):
        raise ValueError(f"VMEC wout variable {name} is not a finite scalar")
    return float(value.item())


def _scalar_alias(file: Any, *names: str) -> float:
    for name in names:
        if name in file.variables:
            return _scalar(file, name)
    raise ValueError(f"VMEC wout file is missing {' or '.join(names)}")


def _metadata(path: Path, netcdf_file: Any) -> float:
    try:
        file = netcdf_file(path, "r", mmap=False)
    except (OSError, TypeError, ValueError) as error:
        raise ValueError(f"cannot read standard VMEC NetCDF file {path}") from error
    with file:
        if int(_scalar(file, "ier_flag")) != 0:
            raise ValueError("VMEC wout reports a failed equilibrium solve")
        if int(_scalar_alias(file, "lasym__logical__", "lasym")) != 0:
            raise ValueError("GLISS does not yet support asymmetric VMEC equilibria")
        if int(_scalar_alias(file, "lrfp__logical__", "lrfp")) != 0:
            raise ValueError("GLISS does not support reversed-field-pinch VMEC output")
        if int(_scalar(file, "signgs")) != -1:
            raise ValueError("GLISS requires the standard VMEC signgs=-1 convention")
        return _scalar_alias(file, "betatotal", "betatot")


def _select_surfaces(transform: Any, radial_surfaces: Optional[int]) -> None:
    if radial_surfaces is None:
        return
    available = int(transform.ns_in)
    if available < 1:
        raise ValueError("booz_xform reported no VMEC half-grid surfaces")
    if available % radial_surfaces:
        raise ValueError(
            "radial_surfaces must divide the VMEC half-grid surface count"
        )
    stride = available // radial_surfaces
    if stride % 2 == 0:
        raise ValueError(
            "radial_surfaces must give an odd centered-subsampling stride"
        )
    first = (stride - 1) // 2
    transform.compute_surfs = np.arange(
        first, available, stride, dtype=np.int32
    )


def _write(
    path: Path,
    converted: Any,
    beta_average: float,
    nfp: int,
    m_max: int,
    n_max: int,
    netcdf_file: Any,
    source_name: str = "",
    transform_resolution: tuple[int, int] = (0, 0),
    radial_resolution: tuple[int, int] = (0, 0),
    booz_xform_version: str = "unknown",
) -> None:
    n_modes = np.concatenate((np.arange(n_max + 1), np.arange(-n_max, 0)))
    with netcdf_file(path, "w", version=2) as file:
        file.createDimension("s", converted.s.size)
        file.createDimension("m", m_max + 1)
        file.createDimension("n", n_modes.size)
        file.gliss_schema = b"gvec-cas3d-export"
        file.gliss_schema_version = b"1"
        file.stellarator_symmetry = b"True"
        file.creator = b"gliss.convert_vmec"
        file.vmec_source = source_name.encode("utf-8")
        file.booz_xform_mboz = transform_resolution[0]
        file.booz_xform_nboz = transform_resolution[1]
        file.vmec_half_grid_surfaces = radial_resolution[0]
        file.booz_xform_surfaces = radial_resolution[1]
        file.booz_xform_version = booz_xform_version.encode("ascii")
        for name, value in converted.residuals.items():
            setattr(file, f"conversion_residual_{name}", np.float64(value))
        for name, value in (("N_FP", nfp), ("winding", 1)):
            variable = file.createVariable(name, "i", ())
            variable[...] = value
        variable = file.createVariable("beta_avg", "d", ())
        variable[...] = beta_average
        for name, values, code, dimensions in (
            ("m", np.arange(m_max + 1), "i", ("m",)),
            ("n", n_modes, "i", ("n",)),
            ("s", converted.s, "d", ("s",)),
            ("rho", np.sqrt(converted.s), "d", ("s",)),
        ):
            variable = file.createVariable(name, code, dimensions)
            variable[:] = values
        for name, values in converted.profiles.items():
            variable = file.createVariable(name, "d", ("s",))
            variable[:] = values
        for name, (cosine, sine) in converted.harmonics.items():
            suffix, values = ("mnc", cosine) if name in _EVEN_FIELDS else ("mns", sine)
            variable = file.createVariable(f"{name}_{suffix}", "d", ("s", "m", "n"))
            variable[:] = values


def convert_vmec(
    input_path: PathLike,
    output_path: PathLike,
    *,
    poloidal_max: int = 7,
    toroidal_max: int = 7,
    transform_factor: int = 4,
    radial_surfaces: Optional[int] = None,
    overwrite: bool = False,
) -> Path:
    """Convert a stellarator-symmetric VMEC ``wout`` file for GLISS.

    The result uses the left-handed, one-field-period Boozer convention of
    pyGVEC's CAS3D exporter. ``radial_surfaces`` optionally selects an exact
    centered uniform subset of the VMEC half grid. Existing outputs are
    preserved unless ``overwrite=True``.
    """
    source_path = _path(input_path, "input_path", True)
    destination = _path(output_path, "output_path", False)
    poloidal_max = _integer(poloidal_max, "poloidal_max", 0, 64)
    toroidal_max = _integer(toroidal_max, "toroidal_max", 0, 64)
    transform_factor = _integer(transform_factor, "transform_factor", 2, 16)
    if radial_surfaces is not None:
        radial_surfaces = _integer(
            radial_surfaces, "radial_surfaces", 5, 1_000_000
        )
    if not isinstance(overwrite, bool):
        raise TypeError("overwrite must be a bool")
    if destination.exists() and not overwrite:
        raise FileExistsError(f"output_path already exists: {destination}")
    if destination.exists() and not destination.is_file():
        raise ValueError(f"output_path is not a regular file: {destination}")
    if destination.exists() and source_path.samefile(destination):
        raise ValueError("input_path and output_path must name different files")
    if not destination.exists() and source_path.resolve() == destination.resolve():
        raise ValueError("input_path and output_path must name different files")
    booz_xform, netcdf_file = _dependencies()
    beta_average = _metadata(source_path, netcdf_file)
    transform = booz_xform.Booz_xform()
    transform.verbose = 0
    try:
        transform.read_wout(os.fspath(source_path), True)
    except Exception as error:
        raise RuntimeError(f"cannot load VMEC equilibrium {source_path}") from error
    _select_surfaces(transform, radial_surfaces)
    try:
        transform.mboz = transform_factor * (poloidal_max + 1)
        transform.nboz = transform_factor * (toroidal_max + 1)
        transform.run()
    except Exception as error:
        raise RuntimeError(f"Boozer transformation failed for {source_path}") from error
    if bool(transform.asym):
        raise ValueError("booz_xform reported asymmetric geometry")
    converted = convert_geometry(
        transform, beta_average, poloidal_max, toroidal_max
    )
    descriptor, temporary_name = tempfile.mkstemp(
        dir=destination.parent, prefix=f".{destination.name}.", suffix=".tmp"
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        _write(
            temporary,
            converted,
            beta_average,
            int(transform.nfp),
            poloidal_max,
            toroidal_max,
            netcdf_file,
            source_path.name,
            (int(transform.mboz), int(transform.nboz)),
            (int(transform.ns_in), int(converted.s.size)),
            getattr(booz_xform, "__version__", "unknown"),
        )
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)
    return destination


__all__ = ["convert_vmec"]
