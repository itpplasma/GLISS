"""Versioned JSON interchange for fixed-boundary stability results."""

import re
from typing import Any, Dict, Mapping

import numpy as np

from ._schema_support import (
    SCHEMA_VERSION,
    SCHEMA_VERSIONS,
    discretization_revision,
    fields,
    integer,
    read_json,
    real,
    schema,
    write_json,
)
from .equilibrium import PathLike
from ._stability_input import angular_grid
from ._stability_input import validate_modes as _validate_modes
from .stability import SpectrumResult, StabilityResult
from .solver import SolverTolerances

RESULT_SCHEMA = "gliss.stability.result"
_SPECTRUM_FIELDS = {
    "parity_class",
    "field_periods",
    "modes",
    "degree",
    "angular_resolution",
    "adiabatic_index",
    "density_kg_m3",
    "zero_floor",
    "negative_count",
    "floor_count",
    "lowest_eigenvalue",
    "certificate",
    "eigenpair_residual",
    "eigenpair_resolution",
    "inertia_interval",
    "eigenvector",
    "normal_unknowns",
    "eta_unknowns",
    "mu_unknowns",
    "has_chart_metric",
    "has_eigenvector",
    "eigenvalue_unit",
    "boundary_condition",
    "normalization",
    "coordinate_handedness",
    "fourier_convention",
}
_SPECTRUM_FIELDS_OLD = (_SPECTRUM_FIELDS - {"degree"}) | {
    "radial_quadrature"
}
_SPECTRUM_FIELDS_V2 = _SPECTRUM_FIELDS_OLD | {"solver_tolerances"}
_SPECTRUM_FIELDS_V3 = _SPECTRUM_FIELDS | {"solver_tolerances"}
_SPECTRUM_FIELDS_V5 = _SPECTRUM_FIELDS_V3 | {"discretization_revision"}
_SPECTRUM_FIELDS_V8 = _SPECTRUM_FIELDS_V5 | {
    "configuration_sha256", "equilibrium_sha256"
}


def _spectrum_to_dict(result: SpectrumResult) -> Dict[str, Any]:
    document = {
        "parity_class": result.parity_class,
        "field_periods": result.field_periods,
        "modes": [list(mode) for mode in result.modes],
        "degree": result.degree,
        "angular_resolution": list(result.angular_resolution),
        "adiabatic_index": result.adiabatic_index,
        "density_kg_m3": result.density_kg_m3,
        "zero_floor": result.zero_floor,
        "negative_count": result.negative_count,
        "floor_count": result.floor_count,
        "lowest_eigenvalue": result.lowest_eigenvalue,
        "certificate": result.certificate,
        "eigenpair_residual": result.eigenpair_residual,
        "eigenpair_resolution": result.eigenpair_resolution,
        "inertia_interval": result.inertia_interval,
        "eigenvector": np.asarray(result.eigenvector).tolist(),
        "normal_unknowns": result.normal_unknowns,
        "eta_unknowns": result.eta_unknowns,
        "mu_unknowns": result.mu_unknowns,
        "has_chart_metric": result.has_chart_metric,
        "has_eigenvector": result.has_eigenvector,
        "eigenvalue_unit": result.eigenvalue_unit,
        "boundary_condition": result.boundary_condition,
        "normalization": result.normalization,
        "coordinate_handedness": result.coordinate_handedness,
        "fourier_convention": result.fourier_convention,
    }
    document["solver_tolerances"] = result.solver_tolerances.to_dict()
    document["discretization_revision"] = result.discretization_revision
    document["configuration_sha256"] = result.configuration_sha256
    document["equilibrium_sha256"] = result.equilibrium_sha256
    return document


def _problem_metadata(value: Mapping[str, Any], context: str, version: int) -> tuple:
    try:
        modes = _validate_modes(value["modes"])
    except (TypeError, ValueError) as error:
        raise ValueError(f"{context}.modes: {error}") from error
    if version in (1, 2):
        if value["radial_quadrature"] != "midpoint":
            raise ValueError(f"{context}.radial_quadrature must be 'midpoint'")
        degree = 1
    else:
        degree = integer(value["degree"], f"{context}.degree", 1)
        if degree > 4:
            raise ValueError(f"{context}.degree must be between 1 and 4")
    angular = value["angular_resolution"]
    if not isinstance(angular, list) or len(angular) != 2:
        raise ValueError(f"{context}.angular_resolution must contain two integers")
    resolution = (
        integer(angular[0], f"{context}.angular_resolution[0]", 1),
        integer(angular[1], f"{context}.angular_resolution[1]", 1),
    )
    angular_grid(*resolution)
    return modes, degree, resolution


def _component_vector(value: Mapping[str, Any], context: str) -> tuple:
    counts = tuple(
        integer(value[name], f"{context}.{name}")
        for name in ("normal_unknowns", "eta_unknowns", "mu_unknowns")
    )
    has_vector = value["has_eigenvector"]
    if not isinstance(has_vector, bool):
        raise ValueError(f"{context}.has_eigenvector must be a boolean")
    vector_data = value["eigenvector"]
    if not isinstance(vector_data, list):
        raise ValueError(f"{context}.eigenvector must be an array")
    vector = np.asarray(
        [
            real(item, f"{context}.eigenvector[{i}]")
            for i, item in enumerate(vector_data)
        ],
        dtype=np.float64,
    )
    expected = sum(counts) if has_vector else 0
    if vector.size != expected:
        raise ValueError(
            f"{context}.eigenvector has {vector.size} values; expected {expected}"
        )
    vector.setflags(write=False)
    return counts, has_vector, vector


def _validate_conventions(
    value: Mapping[str, Any], context: str, version: int
) -> None:
    boundaries = ("fixed", "free") if version >= 6 else ("fixed",)
    if value["boundary_condition"] not in boundaries:
        raise ValueError(
            f"{context}.boundary_condition must be one of {boundaries!r}"
        )
    constants = {
        "eigenvalue_unit": "s^-2",
        "normalization": "x.T @ M @ x = 1",
        "fourier_convention": "2*pi*(m*theta - n*zeta/N_T)",
    }
    for name, expected in constants.items():
        if value[name] != expected:
            raise ValueError(f"{context}.{name} must be {expected!r}")
    if value["coordinate_handedness"] not in ("left-handed", "right-handed"):
        raise ValueError(
            f"{context}.coordinate_handedness must be 'left-handed' or "
            "'right-handed'"
        )


def _certificate_components(value: Mapping[str, Any], context: str) -> Dict[str, float]:
    components = {
        name: real(value[name], f"{context}.{name}", 0.0)
        for name in (
            "certificate",
            "eigenpair_residual",
            "eigenpair_resolution",
            "inertia_interval",
        )
    }
    if components["certificate"] != (
        components["inertia_interval"]
        + components["eigenpair_residual"]
        + components["eigenpair_resolution"]
    ):
        raise ValueError(f"{context}.certificate does not equal its components")
    return components


def _spectrum_from_dict(document: Any, index: int, version: int) -> SpectrumResult:
    context = f"result.classes[{index}]"
    if version == 1:
        expected = _SPECTRUM_FIELDS_OLD
    elif version == 2:
        expected = _SPECTRUM_FIELDS_V2
    elif version in (3, 4):
        expected = _SPECTRUM_FIELDS_V3
    elif version < 8:
        expected = _SPECTRUM_FIELDS_V5
    else:
        expected = _SPECTRUM_FIELDS_V8
    value = fields(document, expected, context)
    for name in ("configuration_sha256", "equilibrium_sha256"):
        digest = value.get(name)
        if digest is not None:
            if not isinstance(digest, str) or re.fullmatch(r"[0-9a-f]{64}", digest) is None:
                raise ValueError(
                    f"{context}.{name} must be null or 64 lowercase hex digits"
                )
    parity = integer(value["parity_class"], f"{context}.parity_class", 0)
    if parity not in (0, 1, 2):
        raise ValueError(f"{context}.parity_class must be 0, 1 or 2")
    field_periods = integer(value["field_periods"], f"{context}.field_periods", 1)
    modes, degree, angular_resolution = _problem_metadata(value, context, version)
    gamma = real(value["adiabatic_index"], f"{context}.adiabatic_index", 0.0)
    density = real(value["density_kg_m3"], f"{context}.density_kg_m3")
    floor = real(value["zero_floor"], f"{context}.zero_floor")
    if gamma <= 0.0 or density <= 0.0 or floor <= 0.0:
        raise ValueError(
            f"{context} adiabatic_index, density_kg_m3, and zero_floor "
            "must be positive"
        )
    counts, has_vector, vector = _component_vector(value, context)
    chart_metric = value["has_chart_metric"]
    if not isinstance(chart_metric, bool):
        raise ValueError(f"{context}.has_chart_metric must be a boolean")
    _validate_conventions(value, context, version)
    components = _certificate_components(value, context)
    return SpectrumResult(
        parity_class=parity,
        field_periods=field_periods,
        modes=modes,
        degree=degree,
        angular_resolution=angular_resolution,
        adiabatic_index=gamma,
        density_kg_m3=density,
        zero_floor=floor,
        negative_count=integer(value["negative_count"], f"{context}.negative_count"),
        floor_count=integer(value["floor_count"], f"{context}.floor_count"),
        lowest_eigenvalue=real(
            value["lowest_eigenvalue"], f"{context}.lowest_eigenvalue"
        ),
        certificate=components["certificate"],
        eigenpair_residual=components["eigenpair_residual"],
        eigenpair_resolution=components["eigenpair_resolution"],
        inertia_interval=components["inertia_interval"],
        eigenvector=vector,
        normal_unknowns=counts[0],
        eta_unknowns=counts[1],
        mu_unknowns=counts[2],
        has_chart_metric=chart_metric,
        has_eigenvector=has_vector,
        coordinate_handedness=value["coordinate_handedness"],
        boundary_condition=value["boundary_condition"],
        solver_tolerances=(
            SolverTolerances.from_dict(value["solver_tolerances"])
            if version >= 2
            else SolverTolerances.historical_defaults()
        ),
        discretization_revision=discretization_revision(value, version, context),
        configuration_sha256=value.get("configuration_sha256"),
        equilibrium_sha256=value.get("equilibrium_sha256"),
    )


def stability_result_to_dict(result: StabilityResult) -> Dict[str, Any]:
    """Return a validated canonical versioned result document."""
    if not isinstance(result, StabilityResult):
        raise TypeError("result must be a gliss.StabilityResult")
    version = SCHEMA_VERSION
    if len({item.solver_tolerances for item in result.classes}) != 1:
        raise ValueError("result parity classes have inconsistent solver tolerances")
    document = {
        "schema": RESULT_SCHEMA,
        "schema_version": version,
        "classes": [_spectrum_to_dict(item) for item in result.classes],
    }
    stability_result_from_dict(document)
    return document


def stability_result_from_dict(document: Mapping[str, Any]) -> StabilityResult:
    """Validate and construct a versioned result document."""
    value = fields(document, {"schema", "schema_version", "classes"}, "result")
    version = schema(value, RESULT_SCHEMA, "result", SCHEMA_VERSIONS)
    classes = value["classes"]
    if not isinstance(classes, list) or len(classes) not in (1, 2):
        raise ValueError(
            "result.classes must contain parity classes 1 then 2, or the "
            "coupled class 0"
        )
    result = StabilityResult(
        tuple(
            _spectrum_from_dict(item, index, version)
            for index, item in enumerate(classes)
        )
    )
    parities = tuple(item.parity_class for item in result.classes)
    if parities not in ((1, 2), (0,)):
        raise ValueError(
            "result parity classes must be 1 then 2, or the coupled class 0"
        )
    if parities == (0,) and version < 5:
        raise ValueError("the coupled class 0 requires schema version 5")
    reference = result.classes[0]
    for item in result.classes[1:]:
        shared = (
            "field_periods",
            "modes",
            "degree",
            "angular_resolution",
            "adiabatic_index",
            "density_kg_m3",
            "zero_floor",
            "has_chart_metric",
            "solver_tolerances",
            "discretization_revision",
            "boundary_condition",
            "configuration_sha256",
            "equilibrium_sha256",
        )
        if any(getattr(item, name) != getattr(reference, name) for name in shared):
            raise ValueError("result parity classes have inconsistent problem metadata")
    return result


def write_stability_result(result: StabilityResult, path: PathLike) -> None:
    """Atomically write a canonical versioned result document."""
    write_json(path, stability_result_to_dict(result))


def read_stability_result(path: PathLike) -> StabilityResult:
    """Read and strictly validate a versioned result document."""
    return stability_result_from_dict(read_json(path))
