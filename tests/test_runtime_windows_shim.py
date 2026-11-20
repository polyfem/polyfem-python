"""Check that Windows compatibility settings require explicit opt-in."""

from __future__ import annotations

import importlib
import io
import os
import sys
import unittest
import unittest.mock
from pathlib import Path

_REPO = Path(__file__).resolve().parents[1]
if str(_REPO / "python") not in sys.path:
    sys.path.insert(0, str(_REPO / "python"))

from polyfempy.runtime._runtime import (  # noqa: E402
    configure_windows_runtime,
    should_auto_configure_windows,
)


class PublicSurfaceTests(unittest.TestCase):
    def test_function_is_not_reexported_from_polyfempy_runtime(self):
        """Advanced callers use ``polyfempy.runtime._runtime`` explicitly."""
        runtime_mod = importlib.import_module("polyfempy.runtime")
        self.assertFalse(hasattr(runtime_mod, "configure_windows_runtime"))
        self.assertEqual(runtime_mod.__all__, ["solve", "Result"])

    def test_returns_applied_and_skipped_keys(self):
        result = configure_windows_runtime()
        self.assertIn("applied", result)
        self.assertIn("skipped", result)
        self.assertIsInstance(result["applied"], list)
        self.assertIsInstance(result["skipped"], list)

    def test_dead_private_stream_wrapper_is_not_present(self):
        runtime_mod = importlib.import_module("polyfempy.runtime._runtime")
        self.assertFalse(hasattr(runtime_mod, "_wrap_stream_as_utf8"))


class NonWindowsNoOpTests(unittest.TestCase):
    def test_linux_is_a_noop_and_reports_why(self):
        if sys.platform == "win32":
            self.skipTest("This test verifies the non-Windows no-op branch.")
        result = configure_windows_runtime()
        self.assertEqual(result["applied"], [])
        self.assertIn("not_windows", result["skipped"])

    def test_linux_does_not_touch_env_or_stdout(self):
        """On non-Windows, the helper must leave the process environment and
        ``sys.stdout`` / ``sys.stderr`` exactly as they were."""
        if sys.platform == "win32":
            self.skipTest("This test verifies the non-Windows no-op branch.")
        before_stdout = sys.stdout
        before_stderr = sys.stderr
        # Use a placeholder value that's highly unlikely to already be set.
        os.environ.pop("POLYFEMPY_T7_SENTINEL", None)

        configure_windows_runtime()

        self.assertIs(sys.stdout, before_stdout)
        self.assertIs(sys.stderr, before_stderr)
        self.assertNotIn("POLYFEMPY_T7_SENTINEL", os.environ)


class OptOutEnvVarTests(unittest.TestCase):
    def _set_env(self, value):
        if value is None:
            os.environ.pop("POLYFEMPY_SKIP_WINDOWS_AUTOCONFIG", None)
        else:
            os.environ["POLYFEMPY_SKIP_WINDOWS_AUTOCONFIG"] = value

    def setUp(self):
        self._saved = os.environ.get("POLYFEMPY_SKIP_WINDOWS_AUTOCONFIG")
        self._enabled = os.environ.pop("POLYFEMPY_WINDOWS_AUTOCONFIG", None)

    def tearDown(self):
        self._set_env(self._saved)
        if self._enabled is None:
            os.environ.pop("POLYFEMPY_WINDOWS_AUTOCONFIG", None)
        else:
            os.environ["POLYFEMPY_WINDOWS_AUTOCONFIG"] = self._enabled

    def test_default_does_not_auto_configure(self):
        self._set_env(None)
        self.assertFalse(should_auto_configure_windows())

    def test_empty_string_is_treated_as_unset(self):
        self._set_env("")
        self.assertFalse(should_auto_configure_windows())

    def test_truthy_values_disable_auto_configure(self):
        os.environ["POLYFEMPY_WINDOWS_AUTOCONFIG"] = "1"
        for v in ("1", "true", "TRUE", "Yes", "on", "  yes  "):
            with self.subTest(value=v):
                self._set_env(v)
                self.assertFalse(should_auto_configure_windows())

    def test_falsy_values_keep_auto_configure_enabled(self):
        os.environ["POLYFEMPY_WINDOWS_AUTOCONFIG"] = "1"
        for v in ("0", "false", "no", "off", "nope", "anything-else"):
            with self.subTest(value=v):
                self._set_env(v)
                self.assertTrue(should_auto_configure_windows())

    def test_opt_in_accepts_truthy_values(self):
        self._set_env(None)
        for value in ("1", "true", "TRUE", "Yes", "on", "  yes  "):
            with self.subTest(value=value):
                os.environ["POLYFEMPY_WINDOWS_AUTOCONFIG"] = value
                self.assertTrue(should_auto_configure_windows())


class WindowsBranchTests(unittest.TestCase):
    """The Windows-only branches. We patch ``sys.platform`` so these run on
    any host."""

    def setUp(self):
        self._saved_kmp = os.environ.get("KMP_DUPLICATE_LIB_OK")

    def tearDown(self):
        if self._saved_kmp is None:
            os.environ.pop("KMP_DUPLICATE_LIB_OK", None)
        else:
            os.environ["KMP_DUPLICATE_LIB_OK"] = self._saved_kmp

    def test_kmp_is_set_when_not_present(self):
        os.environ.pop("KMP_DUPLICATE_LIB_OK", None)
        with unittest.mock.patch.object(sys, "platform", "win32"), \
                unittest.mock.patch.object(sys, "stdout", io.StringIO()), \
                unittest.mock.patch.object(sys, "stderr", io.StringIO()), \
                unittest.mock.patch("os.system", return_value=0):
            result = configure_windows_runtime()
        self.assertEqual(os.environ.get("KMP_DUPLICATE_LIB_OK"), "TRUE")
        self.assertIn("kmp_duplicate_lib_ok", result["applied"])

    def test_kmp_not_overwritten_when_already_set(self):
        os.environ["KMP_DUPLICATE_LIB_OK"] = "FALSE"  # user's deliberate choice
        with unittest.mock.patch.object(sys, "platform", "win32"), \
                unittest.mock.patch.object(sys, "stdout", io.StringIO()), \
                unittest.mock.patch.object(sys, "stderr", io.StringIO()), \
                unittest.mock.patch("os.system", return_value=0):
            result = configure_windows_runtime()
        self.assertEqual(os.environ.get("KMP_DUPLICATE_LIB_OK"), "FALSE")
        self.assertIn("kmp_already_set", result["skipped"])

    def test_kmp_overwritten_when_force_is_true(self):
        os.environ["KMP_DUPLICATE_LIB_OK"] = "FALSE"
        with unittest.mock.patch.object(sys, "platform", "win32"), \
                unittest.mock.patch.object(sys, "stdout", io.StringIO()), \
                unittest.mock.patch.object(sys, "stderr", io.StringIO()), \
                unittest.mock.patch("os.system", return_value=0):
            result = configure_windows_runtime(force=True)
        self.assertEqual(os.environ.get("KMP_DUPLICATE_LIB_OK"), "TRUE")
        self.assertIn("kmp_duplicate_lib_ok", result["applied"])

    def test_stdout_already_utf8_is_skipped(self):
        """If ``sys.stdout`` is already a UTF-8 TextIOWrapper, the helper
        must not re-wrap it (would double-buffer and lose writes)."""
        utf8_stream = io.TextIOWrapper(
            io.BytesIO(), encoding="utf-8", errors="replace"
        )
        with unittest.mock.patch.object(sys, "platform", "win32"), \
                unittest.mock.patch.object(sys, "stdout", utf8_stream), \
                unittest.mock.patch.object(sys, "stderr", utf8_stream), \
                unittest.mock.patch("os.system", return_value=0):
            result = configure_windows_runtime()
            self.assertIs(sys.stdout, utf8_stream)
        self.assertIn("stdout_already_utf8", result["skipped"])
        self.assertIn("stderr_already_utf8", result["skipped"])

    def test_idempotency_second_call_applies_nothing_new(self):
        """Calling the function twice in a row must not keep re-wrapping
        streams or re-setting the KMP flag."""
        os.environ.pop("KMP_DUPLICATE_LIB_OK", None)

        class _FakeStream:
            def __init__(self):
                self.buffer = io.BytesIO()

        fake_out = _FakeStream()
        fake_err = _FakeStream()
        with unittest.mock.patch.object(sys, "platform", "win32"), \
                unittest.mock.patch.object(sys, "stdout", fake_out), \
                unittest.mock.patch.object(sys, "stderr", fake_err), \
                unittest.mock.patch("os.system", return_value=0):
            first = configure_windows_runtime()
            second = configure_windows_runtime()

        # First call applies the full set of actions.
        self.assertIn("kmp_duplicate_lib_ok", first["applied"])
        self.assertIn("stdout_utf8", first["applied"])
        # The second call keeps the configured streams and environment value.
        self.assertIn("kmp_already_set", second["skipped"])
        self.assertIn("stdout_already_utf8", second["skipped"])
        self.assertIn("stderr_already_utf8", second["skipped"])


if __name__ == "__main__":
    unittest.main()
