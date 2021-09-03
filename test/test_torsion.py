import unittest

import polyfempy as pf
# from .utils import plot

import os
import platform


class TorsionTest(unittest.TestCase):
    def test_run(self):
        root_folder = os.path.join(
            "..", "3rdparty.nosync" if platform.system() == 'Darwin' else "3rdparty", "data")
        dir_path = os.path.dirname(os.path.realpath(__file__))
        mesh_path = os.path.join(dir_path, root_folder, "square_beam_h.HYBRID")

        settings = pf.Settings(
            discr_order=1, pde=pf.PDEs.NonLinearElasticity, nl_solver_rhs_steps=5)

        settings.set_material_params("E", 200)
        settings.set_material_params("nu", 0.35)

        problem = pf.Torsion(
            fixed_boundary=5, turning_boundary=6, axis_coordiante=2, n_turns=0.5)
        settings.problem = problem

        solver = pf.Solver()
        solver.settings(settings)
        solver.load_mesh_from_path(
            mesh_path, normalize_mesh=False, vismesh_rel_area=0.001)

        solver.solve()

        pts, tets, el_id, disp = solver.get_sampled_solution()
        vertices = pts + disp
        mises, _ = solver.get_sampled_mises_avg()

        # plot(vertices, tets, mises)


if __name__ == '__main__':
    unittest.main()
