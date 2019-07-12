import json
import polyfempy

class Settings:
	"""Class that encodes the settings of the solver, it models the input json file"""

	def __init__(self):
		self.discr_order = 1
		self.pressure_discr_order = 1
		self.__is_scalar = True

		self.scalar_formulation = "Laplacian"
		self.tensor_formulation = "LinearElasticity"

		self.normalize_mesh = True

		self.n_refs = 0
		self.vismesh_rel_area = 0.00001
		self.nl_solver_rhs_steps = 1



		self.problem = "Franke"


		self.tend = 1
		self.time_steps = 10


		self.quadrature_order = 4

		self.params = {}
		self.export = {}
		self.advanced_options = {}

		self.problem_params = {}



	def set_problem(self, problem):
		"""Sets the problem, use any of the problems in Problems or the Problem"""

		if isinstance(problem, polyfempy.Problem):
			self.problem = "GenericScalar" if self.__is_scalar else "GenericTensor"
			if problem.rhs is None:
				problem.rhs = 0 if self.__is_scalar else [0, 0, 0]
		else:
			self.problem = problem.name()
		self.problem_params = problem.params()


	def set_pde(self, pde):
		"""Sets the PDE to solve, use any of the polyfempy.PDEs"""

		if pde == "NonLinearElasticity":
			pde = "NeoHookean"

		self.__is_scalar = not polyfempy.polyfempy.is_tensor(pde)

		if self.__is_scalar:
			self.scalar_formulation = pde
		else:
			self.tensor_formulation = pde

	def set_material_params(self, name, value):
		"""set the material parameters, for instance set_material_params("E", 200) sets the Young's modulus E to 200. See https://polyfem.github.io/documentation/#formulations for full list"""

		self.params[name] = value


	def set_vtu_export_path(self, path, bounda_only=False):
		"""Sets the path to export a vtu file with the results, use bounda_only to export only one layer of the mesh in 3d"""

		self.export["vis_mesh"] = path
		self.export["vis_boundary_only"] = bounda_only


	def set_wireframe_export_path(self, path):
		"""Sets the path to export a wireframe of the mesh"""

		self.export["wire_mesh"] = path


	def set_isolines_export_path(self, path):
		"""Sets the path to export the isolines of the solution"""

		self.export["iso_mesh"] = path


	def set_solution_export_path(self, path):
		"""Sets the path to save the solution"""

		self.export["solution"] = path


	def set_advanced_option(self, key, value):
		"""Used to set any advanced option not present in this class, for instance set_advanced_option("use_spline",True), see https://polyfem.github.io/documentation/ for full list"""

		self.advanced_options[key] = value

	def __str__(self):
		"""stringygied json description of this class, used to run the solver"""

		tmp = dict(
			(key, value)
			for (key, value) in self.__dict__.items())
		tmp.pop('advanced_options', None)
		tmp.update(self.advanced_options)

		return json.dumps(tmp, sort_keys=True, indent=4)

	def serialize(self):
		"""stringyfied json description of this class, used to run the solver"""

		return str(self)
