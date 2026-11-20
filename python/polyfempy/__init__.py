"""polyfempy Python package.

This repository contains:
- A pure-Python solve/runtime layer (`polyfempy.runtime.*`)
- A compiled C++ extension submodule built by CMake (nanobind backend)

Why this shape:
- The Python import path should remain stable: `import polyfempy as pf`

Important:
- The compiled extension is expected to be available as `polyfempy.polyfempy`
  (i.e., a submodule inside this package). When present, we re-export its public
  symbols (e.g., `Solver`, `version`, etc.) at the
  package top-level for convenience.
"""

# High-level (pure Python) solve/runtime surface
from .runtime import *  # noqa: F403

# Low-level C++ extension (optional at import time; required for real compute)
_cpp_backend_available = False
_cpp_backend_error = None
try:
    import importlib
    import importlib.util

    # This is the compiled extension module produced by CMake:
    # full import name:  polyfempy.polyfempy
    _spec = importlib.util.find_spec(__name__ + ".polyfempy")
    if _spec is None:
        raise ModuleNotFoundError(f"No module named '{__name__}.polyfempy' (C++ backend not built)")

    _core = importlib.import_module(__name__ + ".polyfempy")  # type: ignore

    # Re-export all public names from the C++ extension at package top-level.
    for _name in dir(_core):
        if not _name.startswith("_"):
            globals()[_name] = getattr(_core, _name)

    _cpp_backend_available = True
except Exception as e:
    # Keep import working even if the extension is not built/installed yet.
    _cpp_backend_error = e
    _core = None


def cpp_backend_available() -> bool:
    """Return True if the compiled C++ backend submodule is importable."""
    return _cpp_backend_available


def cpp_backend_error():
    """Return the import error (if any) for the compiled C++ backend."""
    return _cpp_backend_error
