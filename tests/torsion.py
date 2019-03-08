import unittest

import polyfempy as pf
# from .utils import plot

import os


class TorsionTest(unittest.TestCase):
    def test_run(self):
        dir_path = os.path.dirname(os.path.realpath(__file__))
        mesh_path = os.path.join(dir_path, "../3rdparty/data/square_beam_h.HYBRID")

        settings = pf.Settings()
        settings.discr_order = 1
        settings.normalize_mesh = False



        settings.set_material_params("E", 200)
        settings.set_material_params("nu", 0.35)

        settings.nl_solver_rhs_steps = 5
        settings.tensor_formulation = pf.TensorFormulations.NeoHookean

        problem = pf.Torsion()
        problem.fixed_boundary = 5
        problem.turning_boundary = 6
        problem.axis_coordiante = 2
        problem.n_turns = 0.5


        settings.vismesh_rel_area = 0.001

        settings.set_problem(problem)


        solver = pf.Solver()


        solver.settings(str(settings))
        solver.load_mesh_from_path(mesh_path)

        solver.solve()

        [pts, tets, disp] = solver.get_sampled_solution()
        vertices = pts + disp
        mises, _ = solver.get_sampled_mises_avg()

        # plot(vertices, tets, mises)


if __name__ == '__main__':
    unittest.main()
