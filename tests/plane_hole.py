import unittest

import polyfempy as pf
# from .utils import plot
import os


class BendingTest(unittest.TestCase):
    def test_run(self):
        dir_path = os.path.dirname(os.path.realpath(__file__))
        mesh_path = os.path.join(dir_path, "../3rdparty/data/plane_hole.obj")

        settings = pf.Settings()
        settings.discr_order = 1
        settings.normalize_mesh = True
        # settings.vismesh_rel_area = 0.1


        settings.set_material_params("E", 210000)
        settings.set_material_params("nu", 0.3)


        settings.tensor_formulation = pf.TensorFormulations.LinearElasticity

        problem = pf.GenericTensor()
        problem.add_dirichlet_value(1, [0, 0], [True, False])
        problem.add_dirichlet_value(4, [0, 0], [False, True])

        problem.add_neumann_value(3, [100, 0])

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
