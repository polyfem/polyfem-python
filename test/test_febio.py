import unittest

import platform
import polyfempy as pf
# from .utils import plot
import os


class FeBioTest(unittest.TestCase):
    def test_run(self):
        root_folder = os.path.join(
            "..", "3rdparty.nosync" if platform.system() == 'Darwin' else "3rdparty", "data")

        dir_path = os.path.dirname(os.path.realpath(__file__))
        febio_file = os.path.join(dir_path, root_folder, "test.feb")

        opts = {
            "solver_type": "Eigen::SparseLU",
            "solver_params": {
                "gradNorm": 1e-1,
                "nl_iterations": 10
            }
        }

        pf.solve_febio(
            febio_file,
            output_path="test.vtu",
            log_level=1,
            **opts)


if __name__ == '__main__':
    unittest.main()
