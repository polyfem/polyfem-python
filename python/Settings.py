import json


class Settings:
	def __init__(self):
		self.discr_order = 1
		self.pressure_discr_order = 1

		self.scalar_formulation = "Laplacian"
		self.tensor_formulation = "LinearElasticity"

		self.normalize_mesh = True

		self.n_refs = 0
		self.vismesh_rel_area = 0.00001



		self.problem = "Franke"


		self.tend = 1
		self.time_steps = 10


		self.quadrature_order = 4

		self.params = {}
		self.export = {}
		self.advanced_options = {}

		self.problem_params = {}


	def set_problem(self, problem):
		self.problem = problem.name()
		self.problem_params = problem.params()


	def set_material_params(self, name, value):
		self.params[name] = value


	def set_vtu_export_path(self, path, bounda_only=False):
		self.export["vis_mesh"] = path
		self.export["vis_boundary_only"] = bounda_only


	def set_wireframe_export_path(self, path):
		self.export["wire_mesh"] = path


	def set_isolines_export_path(self, path):
		self.export["iso_mesh"] = path


	def set_solution_export_path(self, path):
		self.export["solution"] = path


	def set_advanced_option(self, key, value):
		self.advanced_options[key] = value


	def serialize(self):
		return json.dumps(self, default=lambda o: o.__dict__, sort_keys=True, indent=4)
