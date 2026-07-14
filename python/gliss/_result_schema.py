"""Versioned JSON interchange for fixed-boundary stability results."""

from typing import Any, Dict, Mapping

import numpy as np

from ._schema_support import (
    SCHEMA_VERSION,
    fields,
    integer,
    read_json,
    real,
    schema,
    write_json,
)
from .equilibrium import PathLike
from ._stability_input import validate_modes as _validate_modes
from .stability import (
    _QUADRATURE,
    SpectrumResult,
    StabilityResult,
)
from .solver import SolverTolerances

RESULT_SCHEMA = "gliss.stability.result"
_SPECTRUM_FIELDS = {
    "parity_class",
    "field_periods",
    "modes",
    "radial_quadrature",
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
_SPECTRUM_FIELDS_V2 = _SPECTRUM_FIELDS | {"solver_tolerances"}


def _spectrum_to_dict(result: SpectrumResult) -> Dict[str, Any]:
    document = {
        "parity_class": result.parity_class,
        "field_periods": result.field_periods,
        "modes": [list(mode) for mode in result.modes],
        "radial_quadrature": result.radial_quadrature,
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
    if result.solver_tolerances != SolverTolerances.historical_defaults():
        document["solver_tolerances"] = result.solver_tolerances.to_dict()
    return document


def _problem_metadata(value: Mapping[str, Any], context: str) -> tuple:
    try:
        modes = _validate_modes(value["modes"])
    except (TypeError, ValueError) as error:
        raise ValueError(f"{context}.modes: {error}") from error
    quadrature = value["radial_quadrature"]
    if quadrature not in _QUADRATURE:
        raise ValueError(f"{context}.radial_quadrature must be 'midpoint'")
    angular = value["angular_resolution"]
    if not isinstance(angular, list) or len(angular) != 2:
        raise ValueError(f"{context}.angular_resolution must contain two integers")
    resolution = (
        integer(angular[0], f"{context}.angular_resolution[0]", 1),
        integer(angular[1], f"{context}.angular_resolution[1]", 1),
    )
    return modes, quadrature, resolution


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


def _validate_conventions(value: Mapping[str, Any], context: str) -> None:
    constants = {
        "eigenvalue_unit": "s^-2",
        "boundary_condition": "fixed",
        "normalization": "x.T @ M @ x = 1",
        "coordinate_handedness": "left-handed",
        "fourier_convention": "2*pi*(m*theta - n*zeta/N_T)",
    }
    for name, expected in constants.items():
        if value[name] != expected:
            raise ValueError(f"{context}.{name} must be {expected!r}")


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
    expected = _SPECTRUM_FIELDS_V2 if version == 2 else _SPECTRUM_FIELDS
    value = fields(document, expected, context)
    parity = integer(value["parity_class"], f"{context}.parity_class", 1)
    if parity not in (1, 2):
        raise ValueError(f"{context}.parity_class must be 1 or 2")
    field_periods = integer(value["field_periods"], f"{context}.field_periods", 1)
    modes, quadrature, angular_resolution = _problem_metadata(value, context)
    gamma = real(value["adiabatic_index"], f"{context}.adiabatic_index", 0.0)
    density = real(value["density_kg_m3"], f"{context}.density_kg_m3")
    floor = real(value["zero_floor"], f"{context}.zero_floor")
    if density <= 0.0 or floor <= 0.0:
        raise ValueError(f"{context} density_kg_m3 and zero_floor must be positive")
    counts, has_vector, vector = _component_vector(value, context)
    chart_metric = value["has_chart_metric"]
    if not isinstance(chart_metric, bool):
        raise ValueError(f"{context}.has_chart_metric must be a boolean")
    _validate_conventions(value, context)
    components = _certificate_components(value, context)
    return SpectrumResult(
        parity_class=parity,
        field_periods=field_periods,
        modes=modes,
        radial_quadrature=quadrature,
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
        solver_tolerances=(
            SolverTolerances.from_dict(value["solver_tolerances"])
            if version == 2
            else SolverTolerances.historical_defaults()
        ),
    )


def stability_result_to_dict(result: StabilityResult) -> Dict[str, Any]:
    """Return a validated canonical versioned result document."""
    if not isinstance(result, StabilityResult):
        raise TypeError("result must be a gliss.StabilityResult")
    version = 2 if any(
        item.solver_tolerances != SolverTolerances.historical_defaults()
        for item in result.classes
    ) else SCHEMA_VERSION
    if version == 2 and len({item.solver_tolerances for item in result.classes}) != 1:
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
    version = schema(value, RESULT_SCHEMA, "result", (1, 2))
    classes = value["classes"]
    if not isinstance(classes, list) or len(classes) != 2:
        raise ValueError("result.classes must contain parity classes 1 then 2")
    result = StabilityResult(
        (
            _spectrum_from_dict(classes[0], 0, version),
            _spectrum_from_dict(classes[1], 1, version),
        )
    )
    if tuple(item.parity_class for item in result.classes) != (1, 2):
        raise ValueError("result parity classes must be 1 then 2")
    reference = result.classes[0]
    for item in result.classes[1:]:
        shared = (
            "field_periods",
            "modes",
            "radial_quadrature",
            "angular_resolution",
            "adiabatic_index",
            "density_kg_m3",
            "zero_floor",
            "has_chart_metric",
            "solver_tolerances",
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
