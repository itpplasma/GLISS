"""Python skeleton for GLISS, the Global Linear Ideal Stability Solver.

``version()`` and the Mercier stability bindings (``mercier_profile``,
``mercier_objective``) proving ctypes linkage to the compiled Fortran C API
(``src/gliss_capi.f90``, CMake target ``gliss_c``) are implemented. Build the
shared library with CMake before calling any of them.
"""
import ctypes
import ctypes.util
import os

__version__ = "0.1.0"

_LIBRARY_NAME = "gliss_c"


def _load_library():
    path = os.environ.get("GLISS_LIB") or ctypes.util.find_library(_LIBRARY_NAME)
    if not path:
        raise OSError(
            "libgliss_c shared library not found. Build it with CMake "
            "(target gliss_c) and set GLISS_LIB to the built path, e.g. "
            "build/libgliss_c.so, or add its directory to LD_LIBRARY_PATH."
        )
    try:
        return ctypes.CDLL(path)
    except OSError as error:
        raise OSError(f"failed to load {path}: {error}") from error


def version():
    """Return the compiled library version via the gliss_version C symbol."""
    library = _load_library()
    library.gliss_version.argtypes = (ctypes.c_char_p, ctypes.c_int)
    library.gliss_version.restype = None
    buffer = ctypes.create_string_buffer(32)
    library.gliss_version(buffer, len(buffer))
    return buffer.value.decode("ascii")


from .mercier import mercier_objective, mercier_profile  # noqa: E402

__all__ = ["version", "mercier_profile", "mercier_objective"]
