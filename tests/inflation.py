import unittest

import polyfempy as pf
import numpy as np

import os


class InflationTest(unittest.TestCase):
	def test_run(self):
		solver = pf.Solver()

		#some setup
		dir_path = os.path.dirname(os.path.realpath(__file__))
		mesh_path = os.path.join(dir_path, "../3rdparty/data/circle2.msh")
		print(mesh_path)
		output = "inflation.obj"

		settings = pf.Settings()
		settings.discr_order = 2
		settings.normalize_mesh = True
		settings.vismesh_rel_area = 0.00001
		settings.scalar_formulation = pf.ScalarFormulations.Laplacian

		problem = pf.GenericScalar()
		problem.add_dirichlet_value("all", 10)
		problem.rhs = 0
		settings.set_problem(problem)



		solver.settings(settings.serialize())
		solver.load_mesh(mesh_path)

		solver.solve()
		sol = solver.get_solution()
		##########################################################################

		# now we got the solution of the first laplacian, we use it as rhs for the second one
		# setup zero bc and use sol as rhs
		problem = pf.GenericScalar()
		problem.add_dirichlet_value("all", 0)
		problem.rhs = 0
		settings.set_problem(problem)

		#reload the parameters and mesh
		solver.settings(settings.serialize())
		solver.load_mesh(mesh_path)

		#set the rhs as prev sol
		solver.set_rhs(sol)

		solver.solve()

		#get the solution on a densly sampled mesh
		[vertices, tris, sol] = solver.get_sampled_solution()



		#the dense mesh is disconnected, so we need to connecit it back
		_, res = np.unique(vertices, axis=0, return_inverse=True)
		vertices, resi = np.unique(vertices, axis=0, return_index=True)


		#save obj
		with open(output, "w") as file:
			# use sol as z
			for i in range(len(vertices)):
				file.write("v {} {} {}\n".format(vertices[i][0], vertices[i][1], sol[resi[i]][0]))

			for i in range(len(tris)):
				i0 = res[tris[i][0]]
				i1 = res[tris[i][1]]
				i2 = res[tris[i][2]]
				file.write("f {} {} {}\n".format(i0+1, i1+1, i2+1))


if __name__ == '__main__':
	unittest.main()
