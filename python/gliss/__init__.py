"""Python interface to the Global Linear Ideal Stability Solver.

Equilibria have explicit, context-managed native lifetimes. One-shot Mercier
functions remain available for scripts that evaluate an export once.
"""

import atexit
import ctypes
import ctypes.util
import os
import sys
from contextlib import ExitStack
from importlib.resources import as_file, files

__version__ = "0.0.1"

_ABI_VERSION = 1
_LIBRARY_NAME = "gliss_c"
_RESOURCE_STACK = ExitStack()
atexit.register(_RESOURCE_STACK.close)


def _bundled_library():
    suffix = {"win32": ".dll", "darwin": ".dylib"}.get(sys.platform, ".so")
    prefix = "" if sys.platform == "win32" else "lib"
    resource = files("gliss") / f"{prefix}{_LIBRARY_NAME}{suffix}"
    if not resource.is_file():
        return None
    return _RESOURCE_STACK.enter_context(as_file(resource))


def _open_library():
    configured = os.environ.get("GLISS_LIB")
    path = configured or _bundled_library() or ctypes.util.find_library(_LIBRARY_NAME)
    if not path:
        raise OSError(
            "GLISS shared library not found; install a platform wheel or set "
            "GLISS_LIB to a locally built libgliss_c shared library"
        )
    try:
        return ctypes.CDLL(os.fspath(path))
    except OSError as error:
        raise OSError(
            f"failed to load GLISS shared library {path!s}: {error}"
        ) from error


def _load_library():
    library = _open_library()
    try:
        abi_version = library.gliss_abi_version
    except AttributeError as error:
        raise OSError("GLISS shared library does not expose its ABI version") from error
    abi_version.argtypes = ()
    abi_version.restype = ctypes.c_int
    actual = abi_version()
    if actual != _ABI_VERSION:
        raise OSError(
            f"GLISS shared library ABI version {actual} is incompatible; "
            f"this Python package requires {_ABI_VERSION}"
        )
    return library


def _require_symbols(library, symbols, feature: str) -> None:
    missing = [symbol for symbol in symbols if not hasattr(library, symbol)]
    if missing:
        names = ", ".join(missing)
        raise OSError(
            f"GLISS {feature} requires native symbols {names}; install a "
            "matching Python package and shared library"
        )


def version() -> str:
    """Return the version of the loaded GLISS shared library."""
    library = _load_library()
    library.gliss_version.argtypes = (ctypes.c_char_p, ctypes.c_int)
    library.gliss_version.restype = None
    buffer = ctypes.create_string_buffer(32)
    library.gliss_version(buffer, len(buffer))
    return buffer.value.decode("ascii")


def get_include() -> str:
    """Return the directory containing the wheel-installed C header."""
    resource = files("gliss") / "include" / "gliss.h"
    if not resource.is_file():
        raise FileNotFoundError(
            "bundled gliss.h not found; source builds provide it in include/"
        )
    header = _RESOURCE_STACK.enter_context(as_file(resource))
    return os.fspath(header.parent)


from .equilibrium import (  # noqa: E402
    Equilibrium,
    GlissAllocationError,
    GlissArgumentError,
    GlissCapacityError,
    GlissComputationError,
    GlissError,
    GlissIOError,
    GlissInternalError,
)
from .mercier import mercier_objective, mercier_profile  # noqa: E402
from .stability import (  # noqa: E402
    SpectrumResult,
    StabilityProblem,
    StabilityResult,
)

__all__ = [
    "__version__",
    "Equilibrium",
    "GlissAllocationError",
    "GlissArgumentError",
    "GlissCapacityError",
    "GlissComputationError",
    "GlissError",
    "GlissIOError",
    "GlissInternalError",
    "get_include",
    "mercier_objective",
    "mercier_profile",
    "SpectrumResult",
    "StabilityProblem",
    "StabilityResult",
    "version",
]
