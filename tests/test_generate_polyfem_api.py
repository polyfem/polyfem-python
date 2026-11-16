from __future__ import annotations

import contextlib
import importlib.util
import io
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = ROOT / "tools" / "generate_polyfem_api.py"
GENERATOR_SCRIPT = ROOT / "python-from-jse" / "tools" / "generate_with_overrides.py"
POLYFEM_SOURCE_DIR = ROOT / "polyfem"
POLYFEM_SCHEMA_FILE = POLYFEM_SOURCE_DIR / "json-specs" / "input-spec.json"
GENERATED_DIR = ROOT / "python" / "polyfempy" / "generated_api"


def import_workflow(module_name: str = "generate_polyfem_api_for_test"):
    spec = importlib.util.spec_from_file_location(module_name, SCRIPT_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class GeneratePolyfemApiWorkflowTests(unittest.TestCase):
    def test_default_workflow_runs_generator_from_repo_root(self):
        workflow = import_workflow()

        with (
            mock.patch.object(workflow, "resolve_include_spec_dirs", return_value=[]),
            mock.patch.object(workflow, "missing_required_paths", return_value=[]),
            mock.patch.object(workflow.subprocess, "run") as run_mock,
        ):
            run_mock.return_value = mock.Mock(returncode=0)

            with contextlib.redirect_stdout(io.StringIO()):
                result = workflow.main([])

        self.assertEqual(0, result)
        self.assertEqual(
            [
                [
                    sys.executable,
                    str(GENERATOR_SCRIPT),
                    "--schema-file",
                    str(POLYFEM_SCHEMA_FILE),
                    "--output-file",
                    str(GENERATED_DIR / "generated_class.py"),
                    "--api-output-file",
                    str(GENERATED_DIR / "generated_api.py"),
                    "--manifest-dir",
                    str(GENERATED_DIR),
                    "--relationships",
                    str(ROOT / "generator-config" / "id_relationships.json"),
                    "--api-aliases",
                    str(ROOT / "generator-config" / "api_aliases.json"),
                    "--model-entry",
                    "polyfem.model",
                ],
            ],
            [call.args[0] for call in run_mock.call_args_list],
        )
        self.assertEqual(
            [ROOT],
            [call.kwargs["cwd"] for call in run_mock.call_args_list],
        )

    def test_polyfem_source_dir_override_changes_schema_location(self):
        workflow = import_workflow()
        custom_source_dir = ROOT / "custom-polyfem-source"
        custom_schema_file = custom_source_dir / "json-specs" / "input-spec.json"

        with (
            mock.patch.object(workflow, "resolve_include_spec_dirs", return_value=[]),
            mock.patch.object(workflow, "missing_required_paths", return_value=[]),
            mock.patch.object(workflow.subprocess, "run") as run_mock,
        ):
            run_mock.return_value = mock.Mock(returncode=0)

            with contextlib.redirect_stdout(io.StringIO()):
                result = workflow.main([
                    "--polyfem-source-dir",
                    str(custom_source_dir),
                ])

        self.assertEqual(0, result)
        self.assertIn(
            str(custom_schema_file),
            run_mock.call_args_list[0].args[0],
        )

    def test_include_spec_dir_is_forwarded_to_generator(self):
        workflow = import_workflow()
        include_dir = ROOT / "custom-linked-specs"

        with (
            mock.patch.object(
                workflow,
                "resolve_include_spec_dirs",
                return_value=[include_dir],
            ),
            mock.patch.object(workflow, "missing_required_paths", return_value=[]),
            mock.patch.object(workflow.subprocess, "run") as run_mock,
        ):
            run_mock.return_value = mock.Mock(returncode=0)

            with contextlib.redirect_stdout(io.StringIO()):
                result = workflow.main([
                    "--include-spec-dir",
                    str(include_dir),
                ])

        self.assertEqual(0, result)
        self.assertIn(
            str(include_dir),
            run_mock.call_args_list[0].args[0],
        )

    def test_resolved_schema_file_skips_polysolve_cache_resolution(self):
        workflow = import_workflow()
        resolved_schema_file = ROOT / "build" / "polyfem-resolved-spec.json"

        with (
            mock.patch.object(workflow, "resolve_include_spec_dirs") as resolve_mock,
            mock.patch.object(workflow, "missing_required_paths", return_value=[]) as missing_mock,
            mock.patch.object(workflow.subprocess, "run") as run_mock,
        ):
            run_mock.return_value = mock.Mock(returncode=0)

            with contextlib.redirect_stdout(io.StringIO()):
                result = workflow.main([
                    "--resolved-schema-file",
                    str(resolved_schema_file),
                ])

        self.assertEqual(0, result)
        resolve_mock.assert_not_called()
        missing_mock.assert_called_once_with(
            resolved_schema_file,
            [],
            require_linked_specs=False,
        )
        generator_command = run_mock.call_args_list[0].args[0]
        self.assertIn(str(resolved_schema_file), generator_command)
        self.assertNotIn(str(POLYFEM_SCHEMA_FILE), generator_command)
        self.assertNotIn("--include-spec-dir", generator_command)

    def test_resolved_schema_file_keeps_explicit_include_dirs_without_cache_resolution(self):
        workflow = import_workflow()
        resolved_schema_file = ROOT / "build" / "polyfem-resolved-spec.json"
        include_dir = ROOT / "custom-linked-specs"

        with (
            mock.patch.object(workflow, "resolve_include_spec_dirs") as resolve_mock,
            mock.patch.object(workflow, "missing_required_paths", return_value=[]),
            mock.patch.object(workflow.subprocess, "run") as run_mock,
        ):
            run_mock.return_value = mock.Mock(returncode=0)

            with contextlib.redirect_stdout(io.StringIO()):
                result = workflow.main([
                    "--resolved-schema-file",
                    str(resolved_schema_file),
                    "--include-spec-dir",
                    str(include_dir),
                ])

        self.assertEqual(0, result)
        resolve_mock.assert_not_called()
        generator_command = run_mock.call_args_list[0].args[0]
        self.assertIn(str(resolved_schema_file), generator_command)
        self.assertIn("--include-spec-dir", generator_command)
        self.assertIn(str(include_dir), generator_command)

    def test_missing_resolved_schema_file_points_to_build_dump(self):
        workflow = import_workflow()
        resolved_schema_file = ROOT / "build" / "polyfem-resolved-spec.json"

        def fake_exists(path: Path) -> bool:
            return path != resolved_schema_file

        with (
            mock.patch.object(workflow.Path, "exists", fake_exists),
            mock.patch.object(workflow.subprocess, "run") as run_mock,
        ):
            stderr = io.StringIO()
            with contextlib.redirect_stderr(stderr):
                result = workflow.main([
                    "--resolved-schema-file",
                    str(resolved_schema_file),
                ])

        self.assertEqual(1, result)
        self.assertIn("--resolved-schema-file", stderr.getvalue())
        self.assertIn("embedded/resolved spec", stderr.getvalue())
        run_mock.assert_not_called()

    def test_missing_linked_solver_specs_are_cached_from_polysolve_pin(self):
        workflow = import_workflow()

        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)
            polyfem_source_dir = tmp_dir / "polyfem"
            json_spec_dir = polyfem_source_dir / "json-specs"
            cmake_dir = polyfem_source_dir / "cmake" / "recipes"
            cache_dir = tmp_dir / "build" / "polysolve-json-specs"
            json_spec_dir.mkdir(parents=True)
            cmake_dir.mkdir(parents=True)
            schema_file = json_spec_dir / "input-spec.json"
            schema_file.write_text("{}", encoding="utf-8")
            cache_dir.mkdir(parents=True)
            for name in workflow.LINKED_SOLVER_SPEC_FILES:
                (cache_dir / name).write_text("old unversioned spec", encoding="utf-8")

            def fake_urlretrieve(url, target):
                Path(target).write_text(url, encoding="utf-8")

            with (
                mock.patch.object(workflow, "POLYSOLVE_SPEC_CACHE_DIR", cache_dir),
                mock.patch.object(workflow.urllib.request, "urlretrieve", side_effect=fake_urlretrieve) as download,
            ):
                for revision, downloads in [("abc123", 2), ("def456", 4), ("abc123", 4)]:
                    (cmake_dir / "polysolve.cmake").write_text(
                        f'CPMAddPackage("gh:polyfem/polysolve#{revision}")', encoding="utf-8"
                    )
                    dirs = workflow.resolve_include_spec_dirs(schema_file, [], polyfem_source_dir)
                    self.assertEqual([cache_dir / revision], dirs)
                    self.assertEqual(downloads, download.call_count)
                    for name in workflow.LINKED_SOLVER_SPEC_FILES:
                        self.assertEqual(
                            f"https://raw.githubusercontent.com/polyfem/polysolve/{revision}/{name}",
                            (dirs[0] / name).read_text(encoding="utf-8"),
                        )

    def test_check_workflow_runs_backend_free_checks_after_generation(self):
        workflow = import_workflow()

        with (
            mock.patch.object(workflow, "resolve_include_spec_dirs", return_value=[]),
            mock.patch.object(workflow, "missing_required_paths", return_value=[]),
            mock.patch.object(workflow.subprocess, "run") as run_mock,
        ):
            run_mock.return_value = mock.Mock(returncode=0)

            with contextlib.redirect_stdout(io.StringIO()):
                result = workflow.main(["--check"])

        self.assertEqual(0, result)
        self.assertEqual(
            [
                [
                    sys.executable,
                    str(GENERATOR_SCRIPT),
                    "--schema-file",
                    str(POLYFEM_SCHEMA_FILE),
                    "--output-file",
                    str(GENERATED_DIR / "generated_class.py"),
                    "--api-output-file",
                    str(GENERATED_DIR / "generated_api.py"),
                    "--manifest-dir",
                    str(GENERATED_DIR),
                    "--relationships",
                    str(ROOT / "generator-config" / "id_relationships.json"),
                    "--api-aliases",
                    str(ROOT / "generator-config" / "api_aliases.json"),
                    "--model-entry",
                    "polyfem.model",
                ],
                [
                    sys.executable,
                    "-m",
                    "pytest",
                    "tests/test_examples_public_surface.py",
                    "tests/test_generated_api_example.py",
                ],
            ],
            [call.args[0] for call in run_mock.call_args_list],
        )
        self.assertEqual(
            [ROOT, ROOT],
            [call.kwargs["cwd"] for call in run_mock.call_args_list],
        )
        self.assertNotIn("regenerate_and_test.py", str(run_mock.call_args_list))

    def test_workflow_reports_missing_required_directory(self):
        workflow = import_workflow()

        def fake_exists(path: Path) -> bool:
            return path != ROOT / "python-from-jse"

        with (
            mock.patch.object(workflow, "resolve_include_spec_dirs", return_value=[]),
            mock.patch.object(workflow.Path, "exists", fake_exists),
            mock.patch.object(workflow.subprocess, "run") as run_mock,
        ):
            stderr = io.StringIO()
            with contextlib.redirect_stderr(stderr):
                result = workflow.main([])

        self.assertEqual(1, result)
        self.assertIn("Missing required path", stderr.getvalue())
        self.assertIn("python-from-jse", stderr.getvalue())
        self.assertIn("submodule", stderr.getvalue())
        run_mock.assert_not_called()

    def test_missing_polyfem_schema_points_to_submodule_setup(self):
        workflow = import_workflow()

        def fake_exists(path: Path) -> bool:
            return path != POLYFEM_SCHEMA_FILE

        with (
            mock.patch.object(workflow, "resolve_include_spec_dirs", return_value=[]),
            mock.patch.object(workflow.Path, "exists", fake_exists),
            mock.patch.object(workflow.subprocess, "run") as run_mock,
        ):
            stderr = io.StringIO()
            with contextlib.redirect_stderr(stderr):
                result = workflow.main([])

        self.assertEqual(1, result)
        self.assertIn("polyfem", stderr.getvalue())
        self.assertIn("submodule", stderr.getvalue())
        run_mock.assert_not_called()

    def test_missing_linked_solver_spec_points_to_include_spec_dir(self):
        workflow = import_workflow()

        def fake_exists(path: Path) -> bool:
            return path not in {
                POLYFEM_SCHEMA_FILE.parent / "linear-solver-spec.json",
                POLYFEM_SCHEMA_FILE.parent / "nonlinear-solver-spec.json",
            }

        with (
            mock.patch.object(workflow, "resolve_include_spec_dirs", return_value=[]),
            mock.patch.object(workflow.Path, "exists", fake_exists),
            mock.patch.object(workflow.subprocess, "run") as run_mock,
        ):
            stderr = io.StringIO()
            with contextlib.redirect_stderr(stderr):
                result = workflow.main([])

        self.assertEqual(1, result)
        self.assertIn("linear-solver-spec.json", stderr.getvalue())
        self.assertIn("--include-spec-dir", stderr.getvalue())
        run_mock.assert_not_called()


if __name__ == "__main__":
    unittest.main()
