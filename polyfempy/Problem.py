import json
import types


class Problem:
    """Generic problem problem, scalar or tensor depending on the pde. Warning, this problem needs to be used with the `set_pde` function in settings"""

    def __init__(self, rhs=None, exact=None, is_time_dependent=False):
        self.rhs = rhs
        self.is_time_dependent = is_time_dependent
        self.exact = exact
        self.dirichlet_boundary = []
        self.neumann_boundary = []
        self.pressure_boundary = []

        self.initial_solution = None
        self.initial_velocity = None
        self.initial_acceleration = None

        # self.dirichlet_boundary_lambda = []
        # self.neumann_boundary_lambda = []
        # self.pressure_boundary_lambda = []

    def set_dirichlet_value(self, id, value, is_dirichlet_dim=None, linear_ramp_to=None):
        """set the Dirichlet value value for the sideset id. Note the value must be a scalar, vector in 2D, or 3D depending on the problem. is_dirichlet_dim is a vector of boolean specifying which dimentions are fixed, only for vector-based problems."""
        self.add_dirichlet_value(id, value, is_dirichlet_dim, linear_ramp_to)

    def set_neumann_value(self, id, value, linear_ramp_to=None):
        """set the Neumann value value for the sideset id. Note the value must be a scalar, vector in 2D, or 3D depending on the problem"""
        self.add_neumann_value(id, value, linear_ramp_to)

    def set_pressure_value(self, id, value, linear_ramp_to=None):
        """set the Pressure value value for the sideset id. Note the value must be a scalar"""
        self.add_pressure_value(id, value, linear_ramp_to)

    def add_dirichlet_value(self, id, value, is_dirichlet_dim=None, linear_ramp_to=None):
        """set the Dirichlet value value for the sideset id. Note the value must be a scalar, vector in 2D, or 3D depending on the problem. is_dirichlet_dim is a vector of boolean specifying which dimentions are fixed, only for vector-based problems."""

        tmp = {}
        tmp["id"] = id
        tmp["value"] = value
        if is_dirichlet_dim is not None:
            assert(len(is_dirichlet_dim) == 3 or len(is_dirichlet_dim) == 2)
            assert(len(value) == len(is_dirichlet_dim))
            tmp["dimension"] = is_dirichlet_dim

        if linear_ramp_to is not None:
            tmp["linear_ramp"] = {"to": linear_ramp_to}

        if isinstance(value, types.LambdaType) or isinstance(value, types.FunctionType):
            pass
            # self.dirichlet_boundary_lambda.append(tmp)
        else:
            self.dirichlet_boundary.append(tmp)

    def add_neumann_value(self, id, value, linear_ramp_to=None):
        """set the Neumann value value for the sideset id. Note the value must be a scalar, vector in 2D, or 3D depending on the problem"""

        tmp = {}
        tmp["id"] = id
        tmp["value"] = value

        if linear_ramp_to is not None:
            tmp["linear_ramp"] = {"to": linear_ramp_to}

        if isinstance(value, types.LambdaType) or isinstance(value, types.FunctionType):
            pass
            # self.neumann_boundary_lambda.append(tmp)
        else:
            self.neumann_boundary.append(tmp)

    def add_pressure_value(self, id, value, linear_ramp_to=None):
        """set the Pressure value value for the sideset id. Note the value must be a scalar"""

        tmp = {}
        tmp["id"] = id
        tmp["value"] = value

        if linear_ramp_to is not None:
            tmp["linear_ramp"] = {"to": linear_ramp_to}

        if isinstance(value, types.LambdaType) or isinstance(value, types.FunctionType):
            pass
            # self.pressure_boundary_lambda.append(tmp)
        else:
            self.pressure_boundary.append(tmp)

    def set_initial_solution(self, value):
        """set the initial solution for time dependent problems"""
        self.initial_solution = value

    def set_initial_velocity(self, value):
        """set the initial velocity for time dependent problems"""
        self.initial_velocity = value

    def set_initial_acceleration(self, value):
        """set the initial acceleration for time dependent problems"""
        self.initial_acceleration = value

    def set_velocity(self, id, value, is_dim_fixed=None, linear_ramp_to=None):
        """set the velocity value for the sideset id. Note the value must be a vector in 2D or 3D depending on the problem"""
        self.add_dirichlet_value(id, value, is_dim_fixed, linear_ramp_to)

    def set_displacement(self, id, value, is_dim_fixed=None, linear_ramp_to=None):
        """set the displacement value for the sideset id. Note the value must be a vector in 2D or 3D depending on the problem"""
        self.add_dirichlet_value(id, value, is_dim_fixed, linear_ramp_to)

    def set_force(self, id, value, linear_ramp_to=None):
        """set the force value for the sideset id. Note the value must be a vector in 2D or 3D depending on the problem"""
        self.add_neumann_value(id, value, linear_ramp_to)

    def set_x_symmetric(self, id):
        """x coorinate is fixed, y is allowed to move (Neumann)"""
        self.add_dirichlet_value(id, [0, 0], [True, False])

    def set_y_symmetric(self, id):
        """y coorinate is fixed, x is allowed to move (Neumann)"""
        self.add_dirichlet_value(id, [0, 0], [False, True])

    def set_xy_symmetric(self, id):
        """xy coorinates are fixed, z is allowed to move (Neumann)"""
        self.add_dirichlet_value(id, [0, 0, 0], [True, True, False])

    def set_xz_symmetric(self, id):
        """xz coorinates are fixed, y is allowed to move (Neumann)"""
        self.add_dirichlet_value(id, [0, 0, 0], [True, False, True])

    def set_yz_symmetric(self, id):
        """yz coorinates are fixed, x is allowed to move (Neumann)"""
        self.add_dirichlet_value(id, [0, 0, 0], [False, True, True])

    def params(self):
        """return a dictionary representation of the problem"""
        tmp = dict(
            (key, value)
            for (key, value) in self.__dict__.items())

        if self.initial_solution is None:
            tmp.pop('initial_solution', None)

        if self.initial_velocity is None:
            tmp.pop('initial_velocity', None)

        if self.initial_acceleration is None:
            tmp.pop('initial_acceleration', None)

        # tmp.pop('dirichlet_boundary_lambda', None)
        # tmp.pop('neumann_boundary_lambda', None)
        # tmp.pop('pressure_boundary_lambda', None)

        return tmp
