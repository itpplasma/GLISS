"""Conversion of standard VMEC ``wout`` files to GLISS equilibria."""

import operator
import os
import tempfile
import warnings
from pathlib import Path
from typing import Any, Optional, Tuple, Union

import numpy as np

from ._vmec_geometry import (
    _EVEN_FIELDS,
    _POSITION_FRAME,
    _VMEC_WINDING,
    convert_geometry,
)

PathLike = Union[str, os.PathLike]
# The only VMEC Jacobian sign GLISS accepts; recorded in every export.
_VMEC_SIGNGS = -1


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


def _force_balance_policy(value: Any) -> str:
    if not isinstance(value, str):
        raise TypeError("force_balance_policy must be a string")
    if value not in {"error", "warn"}:
        raise ValueError("force_balance_policy must be 'error' or 'warn'")
    return value


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


def _metadata(path: Path, netcdf_file: Any) -> Tuple[float, bool]:
    """Return the volume-averaged beta and the stellarator symmetry flag."""
    try:
        file = netcdf_file(path, "r", mmap=False)
    except (OSError, TypeError, ValueError) as error:
        raise ValueError(f"cannot read standard VMEC NetCDF file {path}") from error
    with file:
        if int(_scalar(file, "ier_flag")) != 0:
            raise ValueError("VMEC wout reports a failed equilibrium solve")
        symmetric = int(_scalar_alias(file, "lasym__logical__", "lasym")) == 0
        if int(_scalar_alias(file, "lrfp__logical__", "lrfp")) != 0:
            raise ValueError("GLISS does not support reversed-field-pinch VMEC output")
        if int(_scalar(file, "signgs")) != -1:
            raise ValueError("GLISS requires the standard VMEC signgs=-1 convention")
        return _scalar_alias(file, "betatotal", "betatot"), symmetric


def _select_surfaces(transform: Any, radial_surfaces: Optional[int]) -> None:
    if radial_surfaces is None:
        return
    available = int(transform.ns_in)
    if available < 1:
        raise ValueError("booz_xform reported no VMEC half-grid surfaces")
    if available % radial_surfaces:
        raise ValueError("radial_surfaces must divide the VMEC half-grid surface count")
    stride = available // radial_surfaces
    if stride % 2 == 0:
        raise ValueError("radial_surfaces must give an odd centered-subsampling stride")
    first = (stride - 1) // 2
    transform.compute_surfs = np.arange(first, available, stride, dtype=np.int32)


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
    creator: str = "gliss.convert_vmec",
    boozmn_source: str = "",
) -> None:
    n_modes = np.concatenate((np.arange(n_max + 1), np.arange(-n_max, 0)))
    with netcdf_file(path, "w", version=2) as file:
        file.createDimension("s", converted.s.size)
        file.createDimension("m", m_max + 1)
        file.createDimension("n", n_modes.size)
        file.gliss_schema = b"gvec-cas3d-export"
        file.gliss_schema_version = b"1"
        symmetric = converted.stellarator_symmetric
        file.stellarator_symmetry = b"True" if symmetric else b"False"
        file.position_frame = _POSITION_FRAME.encode("ascii")
        file.creator = creator.encode("ascii")
        file.vmec_signgs = np.int32(_VMEC_SIGNGS)
        file.booz_xform_source = boozmn_source.encode("utf-8")
        file.vmec_source = source_name.encode("utf-8")
        file.booz_xform_mboz = transform_resolution[0]
        file.booz_xform_nboz = transform_resolution[1]
        file.vmec_half_grid_surfaces = radial_resolution[0]
        file.booz_xform_surfaces = radial_resolution[1]
        file.booz_xform_version = booz_xform_version.encode("ascii")
        for name, value in converted.residuals.items():
            setattr(file, f"conversion_residual_{name}", np.float64(value))
        for name, value in (("N_FP", nfp), ("winding", _VMEC_WINDING)):
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
            # A symmetric export keeps only the populated parity of each field;
            # an asymmetric one stores both, as the GLISS reader requires.
            if symmetric:
                parts = [("mnc", cosine) if name in _EVEN_FIELDS else ("mns", sine)]
            else:
                parts = [("mnc", cosine), ("mns", sine)]
            for suffix, values in parts:
                variable = file.createVariable(
                    f"{name}_{suffix}", "d", ("s", "m", "n")
                )
                variable[:] = values


def _check_options(
    source_path: Path,
    destination: Path,
    poloidal_max: int,
    toroidal_max: int,
    force_balance_policy: str,
    truncation_tolerance: float,
    overwrite: bool,
) -> Tuple[int, int, str, float]:
    poloidal_max = _integer(poloidal_max, "poloidal_max", 0, 64)
    toroidal_max = _integer(toroidal_max, "toroidal_max", 0, 64)
    force_balance_policy = _force_balance_policy(force_balance_policy)
    if isinstance(truncation_tolerance, bool) or not isinstance(
        truncation_tolerance, (int, float)
    ):
        raise TypeError("truncation_tolerance must be a real number")
    truncation_tolerance = float(truncation_tolerance)
    if not 0.0 < truncation_tolerance < float("inf"):
        raise ValueError("truncation_tolerance must be positive and finite")
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
    return poloidal_max, toroidal_max, force_balance_policy, truncation_tolerance


def _convert_transform(
    transform: Any,
    booz_xform: Any,
    netcdf_file: Any,
    destination: Path,
    beta_average: float,
    stellarator_symmetric: bool,
    poloidal_max: int,
    toroidal_max: int,
    force_balance_policy: str,
    truncation_tolerance: float,
    source_name: str,
    creator: str,
    boozmn_source: str = "",
) -> Path:
    """Check and atomically write a transform that has Boozer spectra."""
    if bool(transform.asym) == stellarator_symmetric:
        raise ValueError(
            "the Boozer transform and the VMEC equilibrium disagree on "
            "stellarator symmetry"
        )
    converted = convert_geometry(
        transform, beta_average, poloidal_max, toroidal_max, stellarator_symmetric
    )
    truncation = converted.residuals.get("truncated_jacobian", 0.0)
    if truncation > truncation_tolerance:
        raise ValueError(
            "the exported position harmonics reproduce the Jacobian only to "
            f"{truncation:.3g} (relative maximum); increase poloidal_max and "
            f"toroidal_max (now {poloidal_max}, {toroidal_max}) or raise "
            "truncation_tolerance"
        )
    positivity = converted.residuals.get("metric_positivity", 1.0)
    if positivity <= 0.0:
        raise ValueError(
            "the exported metric harmonics are not positive definite "
            f"(min det(g)/(g_tt g_zz) = {positivity:.3g}); increase "
            f"poloidal_max and toroidal_max (now {poloidal_max}, "
            f"{toroidal_max})"
        )
    force_balance = converted.residuals["force_balance"]
    if force_balance > 1.0e-2:
        message = (
            "VMEC conversion failed field-identity checks: "
            f"{converted.residuals}"
        )
        if force_balance_policy == "error":
            raise ValueError(message)
        warnings.warn(message, RuntimeWarning, stacklevel=3)
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
            source_name,
            (int(transform.mboz), int(transform.nboz)),
            (int(transform.ns_in), int(converted.s.size)),
            getattr(booz_xform, "__version__", "unknown"),
            creator,
            boozmn_source,
        )
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)
    return destination


def convert_vmec(
    input_path: PathLike,
    output_path: PathLike,
    *,
    poloidal_max: int = 7,
    toroidal_max: int = 7,
    transform_factor: int = 4,
    radial_surfaces: Optional[int] = None,
    force_balance_policy: str = "error",
    truncation_tolerance: float = 0.05,
    overwrite: bool = False,
) -> Path:
    """Convert a VMEC ``wout`` file for GLISS.

    The result uses the left-handed, one-field-period Boozer convention of
    pyGVEC's CAS3D exporter. A stellarator-symmetric file stores the
    populated parity of each field; an asymmetric (``lasym``) file stores
    both parities and ``stellarator_symmetry="False"``, which GLISS reads
    but its parity-class operators refuse. ``radial_surfaces`` optionally
    selects an exact centered uniform subset of the VMEC half grid. A failed
    flux-surface averaged radial force balance (relative residual above
    1e-2), the solvability condition of the Pfirsch-Schlueter equation GLISS
    solves, rejects the conversion by default; the pointwise closure with the
    metric B_s is reported as ``force_balance_pointwise`` only. The explicit
    ``force_balance_policy="warn"`` option retains the diagnostic export for
    convergence studies. GLISS rebuilds the geometry from the truncated
    ``xhat, yhat, zhat`` harmonics only; the conversion therefore rejects a
    maximum relative Jacobian error of that truncated reconstruction above
    ``truncation_tolerance``, which usually calls for larger ``poloidal_max``
    and ``toroidal_max``. Existing outputs are preserved unless
    ``overwrite=True``.
    """
    source_path = _path(input_path, "input_path", True)
    destination = _path(output_path, "output_path", False)
    transform_factor = _integer(transform_factor, "transform_factor", 2, 16)
    if radial_surfaces is not None:
        radial_surfaces = _integer(radial_surfaces, "radial_surfaces", 5, 1_000_000)
    options = _check_options(
        source_path,
        destination,
        poloidal_max,
        toroidal_max,
        force_balance_policy,
        truncation_tolerance,
        overwrite,
    )
    poloidal_max, toroidal_max = options[0], options[1]
    booz_xform, netcdf_file = _dependencies()
    beta_average, symmetric = _metadata(source_path, netcdf_file)
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
    return _convert_transform(
        transform,
        booz_xform,
        netcdf_file,
        destination,
        beta_average,
        symmetric,
        poloidal_max,
        toroidal_max,
        options[2],
        options[3],
        source_path.name,
        "gliss.convert_vmec",
    )


def convert_boozer(
    input_path: PathLike,
    output_path: PathLike,
    *,
    beta_average: Optional[float] = None,
    wout_path: Optional[PathLike] = None,
    poloidal_max: int = 7,
    toroidal_max: int = 7,
    force_balance_policy: str = "error",
    truncation_tolerance: float = 0.05,
    overwrite: bool = False,
) -> Path:
    """Convert a precomputed BOOZ_XFORM ``boozmn`` file for GLISS.

    The Boozer transform is not re-run: ``mboz``, ``nboz`` and the surface
    list are those of the file, which must be a centered uniform subset of
    the VMEC half grid. A ``boozmn`` file does not store the volume-averaged
    beta; pass it as ``beta_average`` or give the parent ``wout_path``,
    whose metadata (solve status, signgs, symmetry and beta) is then also
    checked. The same geometric and force-balance gates as
    :func:`convert_vmec` apply, and both entry points produce identical
    exports for the same transform.
    """
    source_path = _path(input_path, "input_path", True)
    destination = _path(output_path, "output_path", False)
    options = _check_options(
        source_path,
        destination,
        poloidal_max,
        toroidal_max,
        force_balance_policy,
        truncation_tolerance,
        overwrite,
    )
    poloidal_max, toroidal_max = options[0], options[1]
    if (beta_average is None) == (wout_path is None):
        raise ValueError("give exactly one of beta_average and wout_path")
    booz_xform, netcdf_file = _dependencies()
    symmetric = None
    if wout_path is not None:
        beta_average, symmetric = _metadata(
            _path(wout_path, "wout_path", True), netcdf_file
        )
    elif isinstance(beta_average, bool) or not isinstance(
        beta_average, (int, float)
    ):
        raise TypeError("beta_average must be a real number")
    elif not np.isfinite(beta_average) or beta_average < 0.0:
        raise ValueError("beta_average must be finite and nonnegative")
    transform = booz_xform.Booz_xform()
    transform.verbose = 0
    try:
        transform.read_boozmn(os.fspath(source_path))
    except Exception as error:
        raise RuntimeError(f"cannot load BOOZ_XFORM file {source_path}") from error
    if symmetric is None:
        symmetric = not bool(transform.asym)
    return _convert_transform(
        transform,
        booz_xform,
        netcdf_file,
        destination,
        float(beta_average),
        symmetric,
        poloidal_max,
        toroidal_max,
        options[2],
        options[3],
        "" if wout_path is None else _path(wout_path, "wout_path", True).name,
        "gliss.convert_boozer",
        source_path.name,
    )


__all__ = ["convert_boozer", "convert_vmec"]
