import unittest

import numpy as np
import polyfempy as pf
# from .utils import plot


class PythonicTest(unittest.TestCase):
    def test_run(self):
        n_pts = 3

        extend = np.linspace(0, 1, n_pts)
        x, y = np.meshgrid(extend, extend, sparse=False, indexing='xy')
        pts = np.column_stack((x.ravel(), y.ravel()))

        # Create connectivity
        faces = np.ndarray([(n_pts-1)**2, 4], dtype=np.int32)

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

        def sideset(p):
            if p[0] <= 0.001:
                return 4
            return 1

        solver = pf.solve(
            vertices=pts,
            cells=faces,
            sidesets_func=sideset,
            diriclet_bc=[{"id": 4, "value": [0, 0]}],
            materials=[{"id": 0, "E": 2100, "nu": 0.3}],
            rhs=[0, 0.1],
            pde=pf.PDEs.LinearElasticity,
            discr_order=1
        )

        log = solver.get_log()
        print(log["time_solving"])

        # Get the solution
        pts, tris, el_id, bid, fun = solver.get_sampled_solution()


if __name__ == '__main__':
    unittest.main()
