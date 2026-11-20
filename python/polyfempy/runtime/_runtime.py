"""Optional Windows console and OpenMP compatibility settings.

Importing the package leaves the host process unchanged by default. Set
POLYFEMPY_WINDOWS_AUTOCONFIG=1 before import to enable the legacy settings,
or call configure_windows_runtime() explicitly. The older opt-out variable
POLYFEMPY_SKIP_WINDOWS_AUTOCONFIG takes precedence when both are set.
"""

from __future__ import annotations

import io
import os
import sys
from typing import Dict, List

__all__ = ["configure_windows_runtime", "should_auto_configure_windows"]


_TRUTHY = {"1", "true", "yes", "on", "y", "t"}


def should_auto_configure_windows() -> bool:
    """Return whether the caller explicitly enabled the Windows settings."""
    skip = os.environ.get("POLYFEMPY_SKIP_WINDOWS_AUTOCONFIG", "").strip().lower()
    enable = os.environ.get("POLYFEMPY_WINDOWS_AUTOCONFIG", "").strip().lower()
    return enable in _TRUTHY and skip not in _TRUTHY


def configure_windows_runtime(*, force: bool = False) -> Dict[str, List[str]]:
    """Apply Windows-only runtime tweaks to the current Python process.

    The tweaks are:

    - Rewrap ``sys.stdout`` / ``sys.stderr`` as UTF-8 ``TextIOWrapper`` so
      Unicode math symbols do not become mojibake under the default CP-1252
      console. Skipped if the stream is already UTF-8.
    - Run ``chcp 65001`` so the Windows console code page is UTF-8 too.
    - Set ``KMP_DUPLICATE_LIB_OK=TRUE`` to let a process that links multiple
      OpenMP runtimes (e.g. PyTorch + polyfempy) continue instead of aborting.
      Only set when the variable is not already in the environment so we
      never override a deliberate user choice.

    Safe to call multiple times; operations are idempotent unless
    ``force=True`` is passed (which re-wraps ``sys.stdout`` / ``sys.stderr``
    even if they are already UTF-8).

    On non-Windows platforms this is a no-op.

    Returns:
        A dict with two keys:

        - ``applied``: list of action names that this call actually performed
          (e.g. ``["stdout_utf8", "chcp_65001", "kmp_duplicate_lib_ok"]``).
        - ``skipped``: list of reasons for every action we chose not to
          perform (e.g. ``["not_windows"]``, ``["kmp_already_set"]``, or
          ``["stdout_already_utf8"]``).
    """
    applied: List[str] = []
    skipped: List[str] = []

    if sys.platform != "win32":
        skipped.append("not_windows")
        return {"applied": applied, "skipped": skipped}

    # --- UTF-8 stdout / stderr wrapping ---
    for name in ("stdout", "stderr"):
        stream = getattr(sys, name, None)
        already_utf8 = (
            isinstance(stream, io.TextIOWrapper)
            and (getattr(stream, "encoding", "") or "").lower().replace("-", "") == "utf8"
        )
        if already_utf8 and not force:
            skipped.append(f"{name}_already_utf8")
            continue
        buf = getattr(stream, "buffer", None)
        if buf is None:
            skipped.append(f"{name}_no_buffer")
            continue
        try:
            setattr(sys, name, io.TextIOWrapper(buf, encoding="utf-8", errors="replace"))
            applied.append(f"{name}_utf8")
        except Exception as exc:  # pragma: no cover - defensive
            skipped.append(f"{name}_utf8_error:{type(exc).__name__}")

    # --- Windows console code page ---
    try:
        os.system("chcp 65001 >nul 2>&1")
        applied.append("chcp_65001")
    except Exception as exc:  # pragma: no cover - defensive
        skipped.append(f"chcp_error:{type(exc).__name__}")

    # --- OpenMP duplicate-runtime toleration ---
    if "KMP_DUPLICATE_LIB_OK" in os.environ and not force:
        skipped.append("kmp_already_set")
    else:
        os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
        applied.append("kmp_duplicate_lib_ok")

    return {"applied": applied, "skipped": skipped}
