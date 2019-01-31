import polyfempy as pf





mesh_path = "../3rdparty/data/square_beam_h.HYBRID"

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


settings.vismesh_rel_area = 0.00001

settings.set_problem(problem)


solver = pf.Solver()


solver.settings(settings.serialize())
solver.load_mesh(mesh_path)

solver.solve()

solver.export_vtu("torsion.vtu")
