import unittest

import polyfempy as pf
# from .utils import plot
import os
import platform


class BendingTest(unittest.TestCase):
    def test_run(self):
        #     self.run_one(1)
        #     self.run_one(2)

        # def run_one(self, discr_order):
        discr_order = 2
        root_folder = os.path.join("..", "data", "data")

        dir_path = os.path.dirname(os.path.realpath(__file__))
        mesh_path = os.path.join(dir_path, root_folder, "plane_hole.obj")

        settings = pf.Settings(discr_order=discr_order,
                               pde=pf.PDEs.LinearElasticity)

        settings.set_material_params("E", 210000)
        settings.set_material_params("nu", 0.3)

        problem = pf.Problem()
        problem.set_x_symmetric(1)
        problem.set_y_symmetric(4)

        problem.add_neumann_value(3, [100, 0])
        # problem.add_neumann_value(3, lambda x,y,z: [100, 0, 0])

        settings.problem = problem

        solver = pf.Solver()

        solver.settings(settings)
        solver.load_mesh_from_path(mesh_path, normalize_mesh=True)

        solver.solve()

        pts, tets, el_id, bid, disp = solver.get_sampled_solution()
        vertices = pts + disp
        mises, _ = solver.get_sampled_mises_avg()

        solver.compute_errors()

        # plot(vertices, tets, mises)


if __name__ == '__main__':
    unittest.main()
