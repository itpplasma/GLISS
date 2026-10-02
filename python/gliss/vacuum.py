"""Vacuum and conducting-wall model of the physical free-boundary problem."""

import ctypes
import math
import numbers
from dataclasses import dataclass
from typing import Any, Dict, Mapping, Optional, Tuple, Union

import numpy as np

WALL_NONE = 0
WALL_CONFORMAL = 1
WALL_SURFACE = 2


class _VacuumModel(ctypes.Structure):
    _fields_ = [
        ("struct_size", ctypes.c_size_t),
        ("edge_nu", ctypes.c_int32),
        ("edge_nv", ctypes.c_int32),
        ("wall_kind", ctypes.c_int32),
        ("wall_distance", ctypes.c_double),
        ("wall_nu", ctypes.c_int32),
        ("wall_nv", ctypes.c_int32),
        ("wall_xyz", ctypes.POINTER(ctypes.c_double)),
    ]


def _node_count(value: Any, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, numbers.Integral):
        raise TypeError(f"{name} must be an integer")
    value = int(value)
    if value < 3 or value > np.iinfo(np.int32).max:
        raise ValueError(f"{name} must be between 3 and {np.iinfo(np.int32).max}")
    return value


@dataclass(frozen=True, eq=False)
class VacuumModel:
    """Vacuum outside the plasma edge, optionally bounded by an ideal wall.

    ``edge_resolution`` is the (poloidal, toroidal) node count of the
    full-torus edge mesh on which the vacuum boundary-integral equation is
    discretized; each count must exceed twice the largest ``|m|`` and ``|n|`` of the
    mode table. ``wall`` is ``None`` (no wall, vacuum to infinity), a
    positive float (a conformal wall that distance in metres along the
    outward edge normal), or a Cartesian node array of shape (3, nu, nv) in
    metres, poloidal index first, that encloses the plasma.
    """

    edge_resolution: Tuple[int, int] = (32, 32)
    wall: Union[None, float, np.ndarray] = None

    def __post_init__(self) -> None:
        resolution = self.edge_resolution
        if not isinstance(resolution, (tuple, list)) or len(resolution) != 2:
            raise ValueError("edge_resolution must contain two node counts")
        object.__setattr__(
            self,
            "edge_resolution",
            (
                _node_count(resolution[0], "edge_resolution[0]"),
                _node_count(resolution[1], "edge_resolution[1]"),
            ),
        )
        wall = self.wall
        if wall is None:
            return
        if isinstance(wall, numbers.Real) and not isinstance(wall, bool):
            distance = float(wall)
            if not math.isfinite(distance) or distance <= 0.0:
                raise ValueError("a conformal wall distance must be positive")
            object.__setattr__(self, "wall", distance)
            return
        array = np.array(wall, dtype=np.float64, copy=True, order="F")
        if array.ndim != 3 or array.shape[0] != 3:
            raise ValueError("a wall surface must have shape (3, nu, nv)")
        _node_count(array.shape[1], "wall nu")
        _node_count(array.shape[2], "wall nv")
        if not np.all(np.isfinite(array)):
            raise ValueError("wall surface nodes must be finite")
        array.setflags(write=False)
        object.__setattr__(self, "wall", array)

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, VacuumModel):
            return NotImplemented
        if self.edge_resolution != other.edge_resolution:
            return False
        if isinstance(self.wall, np.ndarray) or isinstance(other.wall, np.ndarray):
            return (
                isinstance(self.wall, np.ndarray)
                and isinstance(other.wall, np.ndarray)
                and np.array_equal(self.wall, other.wall)
            )
        return self.wall == other.wall

    @property
    def wall_kind(self) -> str:
        """``"none"``, ``"conformal"`` or ``"surface"``."""
        if self.wall is None:
            return "none"
        if isinstance(self.wall, float):
            return "conformal"
        return "surface"

    def _native(self) -> Tuple[_VacuumModel, Optional[np.ndarray]]:
        model = _VacuumModel(
            struct_size=ctypes.sizeof(_VacuumModel),
            edge_nu=self.edge_resolution[0],
            edge_nv=self.edge_resolution[1],
            wall_kind=WALL_NONE,
            wall_distance=0.0,
        )
        keep = None
        if isinstance(self.wall, float):
            model.wall_kind = WALL_CONFORMAL
            model.wall_distance = self.wall
        elif isinstance(self.wall, np.ndarray):
            keep = np.asfortranarray(self.wall)
            model.wall_kind = WALL_SURFACE
            model.wall_nu = keep.shape[1]
            model.wall_nv = keep.shape[2]
            model.wall_xyz = keep.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        return model, keep

    def to_dict(self) -> Dict[str, Any]:
        """Return the JSON form recorded in configuration documents."""
        document: Dict[str, Any] = {
            "edge_resolution": list(self.edge_resolution),
            "wall_kind": self.wall_kind,
        }
        if isinstance(self.wall, float):
            document["wall_distance_m"] = self.wall
        elif isinstance(self.wall, np.ndarray):
            document["wall_nodes_m"] = self.wall.tolist()
        return document

    @classmethod
    def from_dict(cls, document: Mapping[str, Any]) -> "VacuumModel":
        """Validate the JSON form of a vacuum model."""
        if not isinstance(document, dict):
            raise ValueError("vacuum must be an object")
        kind = document.get("wall_kind")
        expected = {"edge_resolution", "wall_kind"}
        if kind == "conformal":
            expected.add("wall_distance_m")
        elif kind == "surface":
            expected.add("wall_nodes_m")
        elif kind != "none":
            raise ValueError(
                "vacuum.wall_kind must be 'none', 'conformal' or 'surface'"
            )
        if set(document) != expected:
            raise ValueError(
                f"vacuum fields must be {sorted(expected)}; got {sorted(document)}"
            )
        resolution = document["edge_resolution"]
        if not isinstance(resolution, list) or len(resolution) != 2:
            raise ValueError("vacuum.edge_resolution must contain two integers")
        wall: Any = None
        if kind == "conformal":
            wall = document["wall_distance_m"]
            if isinstance(wall, bool) or not isinstance(wall, numbers.Real):
                raise ValueError("vacuum.wall_distance_m must be a number")
        elif kind == "surface":
            wall = document["wall_nodes_m"]
            if not isinstance(wall, list):
                raise ValueError("vacuum.wall_nodes_m must be a nested array")
        try:
            return cls(tuple(resolution), wall)
        except (TypeError, ValueError) as error:
            raise ValueError(f"vacuum: {error}") from error
