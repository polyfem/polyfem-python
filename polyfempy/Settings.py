import json
from logging import raiseExceptions
import polyfempy
import sys

class Settings:
    """Class that encodes the settings of the solver, it models the input json file"""

    def __init__(self, discr_order=1, pressure_discr_order=1, pde=polyfempy.PDEs.Laplacian, contact_problem=False, BDF_order=1, nl_solver_rhs_steps=-1, tend=1, time_steps=10, dhat=0.03):
        self.discr_order = discr_order
        self.pressure_discr_order = pressure_discr_order
        self._is_scalar = True
        self.dhat = dhat

        self.BDF_order = 1

        self.scalar_formulation = "Laplacian"
        self.tensor_formulation = "LinearElasticity"
        self.has_collision = contact_problem
        self.BDF_order = BDF_order

        self.nl_solver_rhs_steps = nl_solver_rhs_steps

        self._problem = "Franke"
        self._python_problem = None

        self.tend = tend
        self.time_steps = time_steps

        self.quadrature_order = 4

        self.params = {}
        self.body_params = []
        self.export = {}
        self.advanced_options = {}

        self.problem_params = {}

        self.pde = pde

        self.selection = None
        raise RuntimeError("Old Version Deprecated. Use version <0.5.2 on conda for the old interface")

    def get_problem(self):
        """Get the problem"""
        return self._problem

    def set_problem(self, problem):
        """Sets the problem, use any of the problems in Problems or the Problem"""
        if isinstance(problem, str):
            self._problem = problem
            return

        if isinstance(problem, polyfempy.Problem):
            self._problem = "GenericScalar" if self._is_scalar else "GenericTensor"
            if problem.rhs is None:
                problem.rhs = 0 if self._is_scalar else [0, 0, 0]
        else:
            self._problem = problem.name()
        self.problem_params = problem.params()
        self._python_problem = problem

    def get_pde(self, pde):
        """Get the PDE"""
        if self._is_scalar:
            return self.scalar_formulation
        else:
            self.tensor_formulation

    def set_pde(self, pde):
        """Sets the PDE to solve, use any of the polyfempy.PDEs"""

        if pde == "NonLinearElasticity":
            pde = "NeoHookean"

        self._is_scalar = not polyfempy.polyfempy.is_tensor(pde)

        if self._is_scalar:
            self.scalar_formulation = pde
        else:
            self.tensor_formulation = pde

    def set_material_params(self, name, value):
        """set the material parameters, for instance set_material_params("E", 200) sets the Young's modulus E to 200. See https://polyfem.github.io/documentation/#formulations for full list"""

        self.params[name] = value

    def set_body_params(self, id, **kwargs):
        """set the material parameters, for a body id. For instance set_body_params(1, E=200, nu=0.3, rho=1000) sets the Young's modulus E to 200, nu=0.3 and density=100 body body 1. See https://polyfem.github.io/documentation/#formulations for full list"""

        tmp = {}
        tmp["id"] = id
        for key, value in kwargs.items():
            tmp[key] = value

        self.body_params.append(tmp)

    def set_vtu_export_path(self, path, boundary_only=False):
        """Sets the path to export a vtu file with the results, use boundary_only to export only one layer of the mesh in 3d"""

        self.export["vis_mesh"] = path
        self.export["vis_boundary_only"] = boundary_only

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
        tmp.pop('_problem', None)
        tmp.pop('_Settings_is_scalar', None)
        tmp.pop('_python_problem', None)
        tmp.pop('selection', None)

        if self.selection is not None:
            tmp["body_ids"] = self.selection.body_ids
            tmp["boundary_sidesets"] = self.selection.boundary_sidesets

        if len(self.body_params) <= 0:
            tmp.pop("body_params", None)

        tmp["problem"] = self.problem
        tmp.update(self.advanced_options)

        json_str = json.dumps(tmp, sort_keys=True, indent=4)

        return json_str

    def serialize(self):
        """stringyfied json description of this class, used to run the solver"""

        return str(self)

    problem = property(get_problem, set_problem)
    pde = property(get_pde, set_pde)
