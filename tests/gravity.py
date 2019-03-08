import unittest

import numpy as np
import polyfempy as pf
# from .utils import plot


class Gravity(unittest.TestCase):
    def test_run(self):
        n_pts = 3

        extend = np.linspace(0, 1, n_pts)
        x, y = np.meshgrid(extend, extend, sparse=False, indexing='xy')
        pts = np.column_stack((x.ravel(), y.ravel()))


        # Create connectivity
        faces = np.ndarray([(n_pts-1)**2, 4],dtype=np.int32)

        index = 0
        for i in range(n_pts-1):
            for j in range(n_pts-1):
                faces[index, :] = np.array([
                    j + i * n_pts,
                    j+1 + i * n_pts,
                    j+1 + (i+1) * n_pts,
                    j + (i+1) * n_pts
                ])
                index = index + 1

        # create settings
        settings = pf.Settings()
        settings.discr_order = 1

        settings.set_material_params("E", 2100)
        settings.set_material_params("nu", 0.3)

        # We use a linear material model
        settings.tensor_formulation = pf.TensorFormulations.LinearElasticity

        problem = pf.Gravity()
        settings.set_problem(problem)

        solver = pf.Solver()
        solver.settings(str(settings))

        # This time we are using pts and faces instead of loading from the disk
        solver.set_mesh(pts, faces)

        solver.solve()

        #Animation frame, last one
        frame = -1

        # Get the solution
        pts = solver.get_sampled_points_frames()
        tets = solver.get_sampled_connectivity_frames()
        disp = solver.get_sampled_solution_frames()


        # diplace the mesh
        vertices = pts[frame] + disp[frame]

        # and get the stresses
        mises = solver.get_sampled_mises_avg_frames()

        # finally plot
        # plot(vertices, tets[frame], mises[frame])


if __name__ == '__main__':
    unittest.main()
