"""Strict deterministic JSON helpers for GLISS interchange schemas."""

import json
import math
import numbers
import os
import tempfile
from pathlib import Path
from typing import Any, Dict, Mapping, Optional

from .equilibrium import PathLike

# Version 6 adds the free-boundary configuration (boundary_condition "free"
# with a vacuum model) and allows "free" in results. Version 7 adds the
# configuration's radial_cells (null: one cell per equilibrium surface).
# Version 8 records the canonical configuration digest in every native
# spectrum, including the radial mesh and exact vacuum/wall model.
SCHEMA_VERSION = 8
SCHEMA_VERSIONS = (1, 2, 3, 4, 5, 6, 7, 8)

# Revision of the assembled fixed-boundary operator. Documents record it so
# that a configuration is never replayed on a different discretization:
# 0 is the midpoint radial quadrature of schema versions 1 and 2, 1 the Gauss
# FEEC operator of versions 3 and 4, 2 the axis-conforming FEEC space of #16
# (tangential axis weight and the regular third unknown), and 3 the space of
# #35 that ties the leading |m|=1 coefficients of xi^s and eta and integrates
# the axis element in sqrt(s).
DISCRETIZATION_REVISION = 3


def discretization_revision(value: Mapping[str, Any], version: int, context: str) -> int:
    """Return the operator revision recorded by a versioned document."""
    if version in (1, 2):
        return 0
    if version in (3, 4):
        return 1
    revision = integer(
        value["discretization_revision"], f"{context}.discretization_revision"
    )
    if revision > DISCRETIZATION_REVISION:
        raise ValueError(
            f"{context}.discretization_revision {revision} is newer than this "
            f"GLISS ({DISCRETIZATION_REVISION})"
        )
    return revision


def document_path(path: PathLike, operation: str) -> Path:
    try:
        value = os.fspath(path)
    except TypeError as error:
        raise TypeError(
            f"{operation} path must be a string or path-like object"
        ) from error
    if not isinstance(value, str):
        raise TypeError(f"{operation} path must resolve to a string")
    if "\0" in value:
        raise ValueError(f"{operation} path contains a null byte")
    return Path(value)


def write_json(path: PathLike, document: Mapping[str, Any]) -> None:
    destination = document_path(path, "output")
    if not destination.parent.is_dir():
        raise FileNotFoundError(
            f"output directory does not exist: {destination.parent}"
        )
    try:
        text = (
            json.dumps(
                document, allow_nan=False, indent=2, sort_keys=True, ensure_ascii=False
            )
            + "\n"
        )
    except (TypeError, ValueError) as error:
        raise ValueError(f"document is not finite JSON data: {error}") from error
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=destination.parent,
            prefix=f".{destination.name}.",
            suffix=".tmp",
            delete=False,
        ) as output:
            temporary = Path(output.name)
            output.write(text)
            output.flush()
            os.fsync(output.fileno())
        os.replace(temporary, destination)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def read_json(path: PathLike) -> Dict[str, Any]:
    source = document_path(path, "input")
    if not source.exists():
        raise FileNotFoundError(f"input file does not exist: {source}")
    if not source.is_file():
        raise ValueError(f"input path is not a file: {source}")
    try:
        text = source.read_text(encoding="utf-8")
    except UnicodeDecodeError as error:
        raise ValueError(f"{source}: document is not valid UTF-8") from error
    try:
        document = json.loads(text, object_pairs_hook=_unique_object)
    except json.JSONDecodeError as error:
        raise ValueError(
            f"{source}: invalid or truncated JSON at line {error.lineno}, "
            f"column {error.colno}"
        ) from error
    except ValueError as error:
        raise ValueError(f"{source}: {error}") from error
    if not isinstance(document, dict):
        raise ValueError(f"{source}: top-level document must be an object")
    return document


def _unique_object(pairs: list) -> Dict[str, Any]:
    result = {}
    for name, value in pairs:
        if name in result:
            raise ValueError(f"duplicate field {name!r}")
        result[name] = value
    return result


def fields(value: Any, expected: set, context: str) -> Mapping[str, Any]:
    if not isinstance(value, dict):
        raise ValueError(f"{context} must be an object")
    unknown = sorted(set(value) - expected)
    if unknown:
        raise ValueError(f"{context} has unknown field {unknown[0]!r}")
    missing = sorted(expected - set(value))
    if missing:
        raise ValueError(f"{context} is missing field {missing[0]!r}")
    return value


def schema(
    value: Mapping[str, Any], expected: str, context: str, versions=(1,)
) -> int:
    for name in ("schema", "schema_version"):
        if name not in value:
            raise ValueError(f"{context} is missing field {name!r}")
    if value["schema"] != expected:
        raise ValueError(f"{context}.schema must be {expected!r}")
    version = value["schema_version"]
    if isinstance(version, bool) or not isinstance(version, int) or version not in versions:
        expected_versions = " or ".join(str(item) for item in versions)
        raise ValueError(
            f"{context}.schema_version is {version!r}; expected {expected_versions}"
        )
    return version


def integer(value: Any, name: str, minimum: int = 0) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{name} must be an integer")
    if value < minimum:
        raise ValueError(f"{name} must be at least {minimum}")
    return value


def real(value: Any, name: str, minimum: Optional[float] = None) -> float:
    if isinstance(value, bool) or not isinstance(value, numbers.Real):
        raise ValueError(f"{name} must be a real number")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    if minimum is not None and result < minimum:
        raise ValueError(f"{name} must be at least {minimum}")
    return result


def string(value: Any, name: str) -> str:
    if not isinstance(value, str) or not value:
        raise ValueError(f"{name} must be a nonempty string")
    return value
