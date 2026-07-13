"""Input validation shared by fixed-boundary Python interfaces."""

import math
import numbers
import operator
from typing import Any, Sequence, Tuple

import numpy as np

_INT32_MIN = -(2**31)
_INT32_MAX = 2**31 - 1


def real_parameter(value: Any, name: str, allow_zero: bool) -> float:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, numbers.Real):
        raise TypeError(f"{name} must be a real number")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    if allow_zero:
        if result < 0.0:
            raise ValueError(f"{name} must be nonnegative")
    elif result <= 0.0:
        raise ValueError(f"{name} must be positive")
    return result


def mode_integer(value: Any, name: str) -> int:
    if isinstance(value, (bool, np.bool_)):
        raise TypeError(f"{name} must be an integer")
    try:
        result = operator.index(value)
    except TypeError as error:
        raise TypeError(f"{name} must be an integer") from error
    if not _INT32_MIN <= result <= _INT32_MAX:
        raise ValueError(f"{name} must fit a signed 32-bit integer")
    return result


def validate_modes(
    modes: Sequence[Tuple[int, int]],
) -> Tuple[Tuple[int, int], ...]:
    try:
        entries = tuple(modes)
    except TypeError as error:
        raise TypeError("modes must be a sequence of (m, n) pairs") from error
    if not entries:
        raise ValueError("modes must not be empty")
    result = []
    seen = set()
    for index, entry in enumerate(entries):
        if not isinstance(entry, (tuple, list)) or len(entry) != 2:
            raise TypeError(f"modes[{index}] must be an (m, n) pair")
        m_value = mode_integer(entry[0], f"modes[{index}] poloidal mode")
        n_value = mode_integer(entry[1], f"modes[{index}] toroidal mode")
        if m_value < 0:
            raise ValueError(f"modes[{index}] poloidal mode must be nonnegative")
        if m_value == 0 and n_value < 0:
            raise ValueError(f"modes[{index}] axis mode requires nonnegative n")
        pair = (m_value, n_value)
        if pair in seen:
            raise ValueError(f"duplicate mode {pair!r}")
        seen.add(pair)
        result.append(pair)
    return tuple(result)
