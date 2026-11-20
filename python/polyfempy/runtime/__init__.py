"""Forward simulation runtime."""

from .result import Result as Result
from .solve import solve as solve

CORE_RUNTIME = ["solve", "Result"]
__all__ = CORE_RUNTIME
