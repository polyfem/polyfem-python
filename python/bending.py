import polyfempy
import Settings
import Problems



mesh_path = "/../3rdparty/data/square_beam.mesh"
tag_path = "../3rdparty/data/square_beam.txt"

settings = Settings.Settings()
settings.discr_order = 1
settings.normalize_mesh = False



settings.set_material_params("E", 200000)
settings.set_material_params("nu", 0.35)


settings.tensor_formulation = polyfempy.TensorFormulations.LinearElasticity

problem = Problems.GenericTensor()
problem.add_dirichlet_value(2, [0, 0, 0])
problem.add_neumann_value(1,[0, -100, 0])


settings.vismesh_rel_area = 0.00001

settings.set_problem(problem)


solver = polyfempy.Solver()

solver.settings(settings.serialize())
solver.load_mesh(mesh_path, tag_path)

solver.solve()

solver.export_vtu("bending.vtu")
