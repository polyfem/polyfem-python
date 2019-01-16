import polyfempy
import json
import numpy as np

with open("test.json", 'r') as f:
	problem_params = json.load(f)


#some setup
mesh_path = "/Users/teseo/Desktop/Higher_Order_Meshes/apps/inflate/meshes/rocket.msh"
bc_tag_path = "/Users/teseo/Desktop/Higher_Order_Meshes/apps/inflate/meshes/rocket.txt"
output = "test.obj"

problem_params["discr_order"] = 2

#density of the output mesh
problem_params["vismesh_rel_area"] = 0.01

##########################################################################
solver = polyfempy.Solver()
solver.set_log_level(0)
##########################################################################


#first solve problem with 0 rhs and some bc on the edges
problem_params["problem_params"]["dirichlet_boundary"][0]["value"] = 0.01
problem_params["problem_params"]["rhs"] = 0

solver.load_parameters(json.dumps(problem_params))
solver.load_mesh(mesh_path, bc_tag_path)

solver.solve()

sol = solver.get_solution()
##########################################################################

#now we got the solution of the first laplacian, we use it as rhs for the second one

#setup zero bc and use sol as rhs
problem_params["problem_params"]["dirichlet_boundary"][0]["value"] = 0

#reload the parameters and mesh
solver.load_parameters(json.dumps(problem_params))
solver.load_mesh(mesh_path, bc_tag_path)

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
	cnt = 0

	# use sol as z
	for i in range(len(vertices)):
		file.write("v {} {} {}\n".format(vertices[i][0], vertices[i][1], sol[resi[i]][0]))

	for i in range(len(tris)):
		i0 = res[tris[i][0]]
		i1 = res[tris[i][1]]
		i2 = res[tris[i][2]]
		file.write("f {} {} {}\n".format(i0+1, i1+1, i2+1))
