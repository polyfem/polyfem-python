"""Handwritten solve/runtime layer for PolyFEM.

The recommended user-facing surface is intentionally small:

- ``solve`` for running simulations
- ``Result`` for structured solver output

Generated configuration helpers are the preferred direction for new authoring
APIs. This package owns runtime preparation, backend dispatch, and result
objects.
"""

from ._runtime import (
    configure_windows_runtime as _configure_windows_runtime,
    should_auto_configure_windows as _should_auto_configure_windows,
)

if _should_auto_configure_windows():
    _configure_windows_runtime()

from .result import Result as Result
from .solve import solve as solve

CORE_RUNTIME = [
    "solve",
    "Result",
]

__all__ = ["solve", "Result"]
