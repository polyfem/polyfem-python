import unittest

import polyfempy as pf
import numpy as np
import platform

# from .utils import plot

import os


class InflationTest(unittest.TestCase):
	def test_run(self):
		root_folder = os.path.join("..", "3rdparty.nosync" if platform.system() == 'Darwin' else "3rdparty", "data")

		solver = pf.Solver()

		#some setup
		dir_path = os.path.dirname(os.path.realpath(__file__))
		mesh_path = os.path.join(dir_path, root_folder, "circle2.msh")

		settings = pf.Settings(discr_order=2, pde=pf.PDEs.Laplacian)

		problem = pf.Problem(rhs=0)
		problem.add_dirichlet_value("all", 10)
		settings.set_problem(problem)


		solver.settings(settings)
		solver.load_mesh_from_path(mesh_path, normalize_mesh=True, vismesh_rel_area=0.00001)

		solver.solve()
		sol = solver.get_solution()
		##########################################################################

		# now we got the solution of the first laplacian, we use it as rhs for the second one
		# setup zero bc and use sol as rhs
		problem = pf.Problem(rhs=0)
		problem.add_dirichlet_value("all", 0)
		settings.problem = problem

		#reload the parameters and mesh
		solver.settings(settings)
		solver.load_mesh_from_path(mesh_path, normalize_mesh=True, vismesh_rel_area=0.00001)

		#set the rhs as prev sol
		solver.set_rhs(sol)

		solver.solve()

		#get the solution on a densly sampled mesh
		vertices, tris, sol = solver.get_sampled_solution()



		#the dense mesh is disconnected, so we need to connecit it back
		_, res = np.unique(vertices, axis=0, return_inverse=True)
		vertices, resi = np.unique(vertices, axis=0, return_index=True)

		faces = np.ndarray([len(tris), 3], dtype=np.int64)
		vv = np.ndarray([len(vertices), 3], dtype=np.float64)

		for i in range(len(tris)):
			faces[i] = np.array([res[tris[i][0]], res[tris[i][1]], res[tris[i][2]]])

		for i in range(len(vv)):
			vv[i] = np.array([vertices[i][0], vertices[i][1], sol[resi[i]][0]])

		# plot(vv, faces, None)

		# #save obj
		# with open(output, "w") as file:
		# 	# use sol as z
		# 	for i in range(len(vertices)):
		# 		file.write("v {} {} {}\n".format())

		# 	for i in range(len(tris)):
		# 		i0 = res[tris[i][0]]
		# 		i1 = res[tris[i][1]]
		# 		i2 = res[tris[i][2]]
		# 		file.write("f {} {} {}\n".format(i0+1, i1+1, i2+1))


if __name__ == '__main__':
	unittest.main()
