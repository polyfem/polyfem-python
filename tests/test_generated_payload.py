"""Tests for the generated-class payload boundary.

These tests stay in pure Python. They protect the first generated-only layer:
``Root.as_dict()``-style payloads should be prepared for the backend without
constructing a second config object.
"""

from __future__ import annotations

import copy
import sys
import tempfile
import unittest
from pathlib import Path

_REPO = Path(__file__).resolve().parents[1]
if str(_REPO / "python") not in sys.path:
    sys.path.insert(0, str(_REPO / "python"))

from polyfempy.runtime._solve_contract import (  # noqa: E402
    generated_payload_from_config,
    prepare_canonical_solve_input,
    prepare_generated_backend_payload,
)


def _canonical_path_text(path: Path) -> str:
    return str(path.resolve())


class GeneratedPayloadTests(unittest.TestCase):
    def test_generated_payload_helpers_live_in_solve_contract(self):
        self.assertFalse((_REPO / "python" / "polyfempy" / "api" / "_generated_payload.py").exists())

    def test_generated_payload_from_config_requires_dict(self):
        class BadGeneratedConfig:
            def as_dict(self):
                return ["not", "a", "dict"]

        with self.assertRaisesRegex(TypeError, "as_dict\\(\\) must return dict"):
            generated_payload_from_config(BadGeneratedConfig())

    def test_prepare_payload_does_not_mutate_input_and_preserves_backend_fields(self):
        payload = {
            "geometry": [{"mesh": "beam.msh", "enabled": True, "unused": None}],
            "materials": [{"type": "LinearElasticity", "E": 20.0, "nu": 0.3}],
            "solver": {"new_backend_option": False},
            "output": {"directory": "", "stats": False},
        }
        snapshot = copy.deepcopy(payload)

        prepared = prepare_generated_backend_payload(payload)

        self.assertEqual(payload, snapshot)
        self.assertNotIn("unused", prepared["geometry"][0])
        self.assertFalse(prepared["solver"]["new_backend_option"])
        self.assertEqual(prepared["output"]["directory"], "")
        self.assertFalse(prepared["output"]["stats"])

    def test_prepare_payload_leaves_output_schema_validation_to_backend(self):
        payload = {
            "geometry": [{"mesh": "beam.msh"}],
            "materials": [{"type": "LinearElasticity", "E": 20.0, "nu": 0.3}],
            "output": {
                "directory": "out",
                "future_backend_option": {"enabled": True},
            },
        }

        prepared = prepare_generated_backend_payload(payload)

        self.assertEqual(prepared["output"], payload["output"])

    def test_prepare_payload_restores_single_global_material_object_shape(self):
        payload = {
            "materials": [
                {"type": "NeoHookean", "E": 10000.0, "nu": 0.4, "rho": 1000.0}
            ]
        }

        prepared = prepare_generated_backend_payload(payload)

        self.assertEqual(
            {"type": "NeoHookean", "E": 10000.0, "nu": 0.4, "rho": 1000.0},
            prepared["materials"],
        )

    def test_prepare_payload_keeps_id_materials_as_list(self):
        payload = {
            "materials": [
                {"id": 1, "type": "LinearElasticity", "E": 20.0, "nu": 0.3}
            ]
        }

        prepared = prepare_generated_backend_payload(payload)

        self.assertEqual(payload["materials"], prepared["materials"])

    def test_prepare_payload_resolves_relative_geometry_mesh_from_root_path(self):
        with tempfile.TemporaryDirectory() as tmp:
            case_dir = Path(tmp) / "case"
            case_dir.mkdir()
            mesh_file = case_dir / "beam.msh"
            mesh_file.write_text("placeholder", encoding="utf-8")

            payload = {
                "root_path": str(case_dir / "config.json"),
                "geometry": [{"mesh": "beam.msh"}],
                "materials": [{"type": "LinearElasticity", "E": 20.0, "nu": 0.3}],
            }

            prepared = prepare_generated_backend_payload(payload)

        self.assertEqual(
            prepared["geometry"][0]["mesh"],
            _canonical_path_text(mesh_file),
        )

    def test_prepare_payload_resolves_relative_selection_files_from_root_path(self):
        with tempfile.TemporaryDirectory() as tmp:
            case_dir = Path(tmp) / "case"
            case_dir.mkdir()
            volume_file = case_dir / "volume.txt"
            surface_file = case_dir / "surface.txt"
            volume_file.write_text("1", encoding="utf-8")
            surface_file.write_text("2", encoding="utf-8")

            payload = {
                "root_path": str(case_dir / "config.json"),
                "geometry": [
                    {
                        "mesh": "beam.msh",
                        "volume_selection": {"file": "volume.txt"},
                        "surface_selection": {"file": "surface.txt"},
                    }
                ],
            }

            prepared = prepare_generated_backend_payload(payload)

        self.assertEqual(
            prepared["geometry"][0]["volume_selection"]["file"],
            _canonical_path_text(volume_file),
        )
        self.assertEqual(
            prepared["geometry"][0]["surface_selection"]["file"],
            _canonical_path_text(surface_file),
        )

    def test_prepare_payload_resolves_relative_collision_mesh_files_from_root_path(self):
        with tempfile.TemporaryDirectory() as tmp:
            case_dir = Path(tmp) / "case"
            case_dir.mkdir()
            mesh_file = case_dir / "proxy.ply"
            map_file = case_dir / "coarse-to-proxy.hdf5"
            mesh_file.write_text("placeholder", encoding="utf-8")
            map_file.write_text("placeholder", encoding="utf-8")

            payload = {
                "root_path": str(case_dir / "config.json"),
                "contact": {
                    "collision_mesh": {
                        "mesh": "proxy.ply",
                        "linear_map": "coarse-to-proxy.hdf5",
                    }
                },
            }

            prepared = prepare_generated_backend_payload(payload)

        self.assertEqual(
            prepared["contact"]["collision_mesh"]["mesh"],
            _canonical_path_text(mesh_file),
        )
        self.assertEqual(
            prepared["contact"]["collision_mesh"]["linear_map"],
            _canonical_path_text(map_file),
        )

    def test_prepare_payload_restores_uniform_scale_scalar_for_backend(self):
        payload = {
            "geometry": [
                {
                    "mesh": "ball.obj",
                    "transformation": {
                        "translation": [],
                        "rotation": [],
                        "rotation_mode": "xyz",
                        "scale": [0.04],
                    },
                }
            ],
        }

        prepared = prepare_generated_backend_payload(payload)

        self.assertEqual(0.04, prepared["geometry"][0]["transformation"]["scale"])
        self.assertNotIn("translation", prepared["geometry"][0]["transformation"])
        self.assertNotIn("rotation", prepared["geometry"][0]["transformation"])
        self.assertNotIn("rotation_mode", prepared["geometry"][0]["transformation"])

    def test_prepare_payload_drops_generated_solver_advanced_defaults(self):
        payload = {
            "solver": {
                "advanced": {
                    "cache_size": 900000,
                    "lump_mass_matrix": True,
                    "lagged_regularization_weight": 0.0,
                    "lagged_regularization_iterations": 1,
                    "check_inversion": "Discrete",
                    "jacobian_threshold": 0.0,
                    "characteristic_length": -1.0,
                    "characteristic_force_density": 10000.0,
                }
            }
        }

        prepared = prepare_generated_backend_payload(payload)

        self.assertEqual({"lump_mass_matrix": True}, prepared["solver"]["advanced"])

    def test_prepare_payload_keeps_backend_nonlinear_solver_names(self):
        payload = {
            "solver": {
                "nonlinear": {
                    "solver": "Newton",
                    "x_delta_tol": 1e-12,
                    "grad_norm_tol": 1e-5,
                    "rel_grad_norm_tol": 1e-10,
                    "newton_decrement_tol": 0.0,
                    "rel_x_delta_tol": 0.0,
                    "first_grad_norm_tol": 1e-12,
                    "norm_type": "L2",
                    "max_iterations": 500,
                    "allow_out_of_iterations": False,
                    "iterations_per_strategy": 1,
                    "line_search": {
                        "method": "RobustArmijo",
                        "use_grad_norm_tol": 1e-6,
                        "min_step_size": 1e-10,
                        "max_step_size_iter": 30,
                        "min_step_size_final": 1e-20,
                        "max_step_size_iter_final": 100,
                        "default_init_step_size": 1.0,
                        "step_ratio": 0.5,
                    },
                    "advanced": {
                        "f_delta_tol": 1e-5,
                        "f_delta_step_tol": 100,
                        "derivative_along_delta_x_tol": 0,
                        "apply_gradient_fd": "None",
                        "gradient_fd_eps": 1e-7,
                    },
                }
            }
        }

        prepared = prepare_generated_backend_payload(payload)
        nonlinear = prepared["solver"]["nonlinear"]

        self.assertEqual(1e-12, nonlinear["x_delta_tol"])
        self.assertEqual(1e-5, nonlinear["grad_norm_tol"])
        self.assertNotIn("x_delta", nonlinear)
        self.assertNotIn("grad_norm", nonlinear)
        self.assertNotIn("rel_grad_norm_tol", nonlinear)
        self.assertNotIn("newton_decrement_tol", nonlinear)
        self.assertNotIn("rel_x_delta_tol", nonlinear)
        self.assertNotIn("first_grad_norm_tol", nonlinear)
        self.assertNotIn("norm_type", nonlinear)
        self.assertNotIn("max_iterations", nonlinear)
        self.assertEqual({"method": "RobustArmijo"}, nonlinear["line_search"])
        self.assertEqual({"f_delta_tol": 1e-5}, nonlinear["advanced"])

    def test_prepare_payload_accepts_legacy_f_delta_without_mutating_input(self):
        for augmented in (False, True):
            with self.subTest(augmented=augmented):
                solver = {"nonlinear": {"advanced": {"f_delta": 1e-5}}}
                payload = {"solver": {"augmented_lagrangian": solver} if augmented else solver}
                snapshot = copy.deepcopy(payload)

                prepared = prepare_generated_backend_payload(payload)["solver"]
                if augmented:
                    prepared = prepared["augmented_lagrangian"]

                self.assertEqual({"f_delta_tol": 1e-5}, prepared["nonlinear"]["advanced"])
                self.assertEqual(snapshot, payload)

    def test_prepare_payload_prefers_current_f_delta_tol_including_zero(self):
        for value in (0.0, 2e-5):
            with self.subTest(value=value):
                payload = {"solver": {"nonlinear": {"advanced": {
                    "f_delta": 1e-5, "f_delta_tol": value, "track_info": True,
                }}}}

                advanced = prepare_generated_backend_payload(payload)["solver"]["nonlinear"]["advanced"]

                expected = {"track_info": True, "f_delta_tol": value}
                if value == 0:
                    expected.pop("f_delta_tol")
                self.assertEqual(expected, advanced)

    def test_prepare_payload_drops_zero_tolerance_default_for_both_names(self):
        for key in ("f_delta", "f_delta_tol"):
            with self.subTest(key=key):
                payload = {"solver": {"nonlinear": {"advanced": {key: 0.0}}}}

                prepared = prepare_generated_backend_payload(payload)

                self.assertNotIn("advanced", prepared["solver"]["nonlinear"])

    def test_model_solver_preserves_f_delta_tol_for_the_backend(self):
        import polyfem

        model = polyfem.model()
        model.mesh(mesh="triangle.obj")
        model.solver(nonlinear=polyfem.nonlinear(
            advanced=polyfem.nonlinear_advanced(f_delta_tol=1e-5),
        ))
        cfg = model.config()
        snapshot = cfg.as_dict()

        canonical = prepare_canonical_solve_input(vertices=None, cells=None, cfg=cfg)

        self.assertEqual(
            {"f_delta_tol": 1e-5},
            canonical.backend_settings["solver"]["nonlinear"]["advanced"],
        )
        self.assertEqual(snapshot, cfg.as_dict())

    def test_dictionary_tolerance_names_for_file_and_array_meshes(self):
        cases = (
            ({"f_delta": 1e-5}, 1e-5),
            ({"f_delta_tol": 1e-5}, 1e-5),
            ({"f_delta": 1e-5, "f_delta_tol": 2e-5}, 2e-5),
            ({"f_delta": 1e-5, "f_delta_tol": 0.0}, 0.0),
        )
        for advanced, expected in cases:
            for augmented in (False, True):
                for array_mesh in (False, True):
                    with self.subTest(advanced=advanced, augmented=augmented, array_mesh=array_mesh):
                        solver = {"nonlinear": {"advanced": dict(advanced, track_info=True)}}
                        cfg: dict[str, object] = {
                            "solver": {"augmented_lagrangian": solver} if augmented else solver,
                        }
                        vertices, cells = None, None
                        if array_mesh:
                            vertices, cells = [[0, 0], [1, 0], [0, 1]], [[0, 1, 2]]
                        else:
                            cfg["geometry"] = [{"mesh": "triangle.obj"}]
                        snapshot = copy.deepcopy(cfg)

                        settings = prepare_canonical_solve_input(
                            vertices=vertices, cells=cells, cfg=cfg,
                        ).backend_settings["solver"]
                        if augmented:
                            settings = settings["augmented_lagrangian"]

                        self.assertEqual(
                            {"f_delta_tol": expected, "track_info": True},
                            settings["nonlinear"]["advanced"],
                        )
                        self.assertEqual(snapshot, cfg)


if __name__ == "__main__":
    unittest.main()
