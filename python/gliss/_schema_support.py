"""Strict deterministic JSON helpers for GLISS interchange schemas."""

import json
import math
import numbers
import os
import tempfile
from pathlib import Path
from typing import Any, Dict, Mapping, Optional

from .equilibrium import PathLike

SCHEMA_VERSION = 1


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


def schema(value: Mapping[str, Any], expected: str, context: str) -> None:
    if value["schema"] != expected:
        raise ValueError(f"{context}.schema must be {expected!r}")
    version = value["schema_version"]
    if isinstance(version, bool) or not isinstance(version, int) or version != 1:
        raise ValueError(f"{context}.schema_version is {version!r}; expected 1")


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
