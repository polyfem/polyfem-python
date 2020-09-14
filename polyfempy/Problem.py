import json
import types

class Problem:
    """Generic problem problem, scalar or tensor depending on the pde. Warning, this problem needs to be used with the `set_pde` function in settings"""

    def __init__(self, rhs=None, exact=None):
        self.rhs = rhs
        self.exact = exact
        self.dirichlet_boundary = []
        self.neumann_boundary = []
        self.pressure_boundary = []

        self.dirichlet_boundary_lambda = []
        self.neumann_boundary_lambda = []
        self.pressure_boundary_lambda = []

    def set_dirichlet_value(self, id, value, is_dirichlet_dim=None):
        """set the Dirichlet value value for the sideset id. Note the value must be a scalar, vector in 2D, or 3D depending on the problem. is_dirichlet_dim is a vector of boolean specifying which dimentions are fixed, only for vector-based problems."""
        self.add_dirichlet_value(id, value, is_dirichlet_dim)

    def set_neumann_value(self, id, value):
        """set the Neumann value value for the sideset id. Note the value must be a scalar, vector in 2D, or 3D depending on the problem"""
        self.add_neumann_value(id, value)

    def set_pressure_value(self, id, value):
        """set the Pressure value value for the sideset id. Note the value must be a scalar"""
        self.add_pressure_value(id, value)

    def add_dirichlet_value(self, id, value, is_dirichlet_dim=None):
        """set the Dirichlet value value for the sideset id. Note the value must be a scalar, vector in 2D, or 3D depending on the problem. is_dirichlet_dim is a vector of boolean specifying which dimentions are fixed, only for vector-based problems."""

        tmp = {}
        tmp["id"] = id
        tmp["value"] = value
        if is_dirichlet_dim is not None:
            assert(len(is_dirichlet_dim) == 3 or len(is_dirichlet_dim) == 2)
            assert(len(value) == len(is_dirichlet_dim))
            tmp["dimension"] = is_dirichlet_dim

        if isinstance(value, types.LambdaType) or isinstance(value, types.FunctionType):
            pass
            # self.dirichlet_boundary_lambda.append(tmp)
        else:
            self.dirichlet_boundary.append(tmp)

    def add_neumann_value(self, id, value):
        """set the Neumann value value for the sideset id. Note the value must be a scalar, vector in 2D, or 3D depending on the problem"""

        tmp = {}
        tmp["id"] = id
        tmp["value"] = value

        if isinstance(value, types.LambdaType) or isinstance(value, types.FunctionType):
            pass
            # self.neumann_boundary_lambda.append(tmp)
        else:
            self.neumann_boundary.append(tmp)

    def add_pressure_value(self, id, value):
        """set the Pressure value value for the sideset id. Note the value must be a scalar"""

        tmp = {}
        tmp["id"] = id
        tmp["value"] = value

        if isinstance(value, types.LambdaType) or isinstance(value, types.FunctionType):
            pass
            # self.pressure_boundary_lambda.append(tmp)
        else:
            self.pressure_boundary.append(tmp)


    def set_velocity(self, id, value, is_dim_fixed=None):
        """set the velocity value for the sideset id. Note the value must be a vector in 2D or 3D depending on the problem"""
        self.add_dirichlet_value(id, value, is_dim_fixed)

    def set_displacement(self, id, value, is_dim_fixed=None):
        """set the displacement value for the sideset id. Note the value must be a vector in 2D or 3D depending on the problem"""
        self.add_dirichlet_value(id, value, is_dim_fixed)

    def set_force(self, id, value):
        """set the force value for the sideset id. Note the value must be a vector in 2D or 3D depending on the problem"""
        self.add_neumann_value(id, value)

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
        """removing the lambda list from the result"""
        tmp = dict(
            (key, value)
            for (key, value) in self.__dict__.items())
        tmp.pop('dirichlet_boundary_lambda', None)
        tmp.pop('neumann_boundary_lambda', None)
        tmp.pop('pressure_boundary_lambda', None)

        return tmp
