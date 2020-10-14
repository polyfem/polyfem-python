

#################################################
############Scalar###############################
class Franke:
    """Franke problem with exact solution https://polyfem.github.io/documentation/#franke"""

    def name(self):
        return "Franke"

    def params(self):
        return {}


class GenericScalar:
    """Generic scalar problem https://polyfem.github.io/documentation/#genericscalar"""

    def __init__(self):
        self.rhs = 0
        self.exact = None
        self.dirichlet_boundary = []
        self.neumann_boundary = []

    def add_dirichlet_value(self, id, value):
        """add the Dirichlet value value for the sideset id"""

        tmp = {}
        tmp["id"] = id
        tmp["value"] = value
        self.dirichlet_boundary.append(tmp)

    def add_neumann_value(self, id, value):
        """add the Neumann value value for the sideset id"""

        tmp = {}
        tmp["id"] = id
        tmp["value"] = value
        self.neumann_boundary.append(tmp)

    def name(self):
        return "GenericScalar"

    def params(self):
        return self.__dict__


#################################################
############Elasticity###########################
class Gravity:
    """time dependent gravity problem https://polyfem.github.io/documentation/#gravity"""

    def __init__(self, force=0.1):
        self.force = force

    def name(self):
        return "Gravity"

    def params(self):
        return self.__dict__


class Torsion:
    """3D torsion problem, specify which sideset to fix (fixed_boundary) and which one turns turning_boundary https://polyfem.github.io/documentation/#torsionelastic"""

    def __init__(self, axis_coordiante=None, axis_coordinate=2, n_turns=0.5, fixed_boundary=5, turning_boundary=6):
        if axis_coordiante:
            self.axis_coordiante = axis_coordiante
        else:
            self.axis_coordiante = axis_coordinate
        self.n_turns = n_turns
        self.fixed_boundary = fixed_boundary
        self.turning_boundary = turning_boundary

    def name(self):
        return "TorsionElastic"

    def params(self):
        return self.__dict__


class GenericTensor:
    """Generic tensor problem https://polyfem.github.io/documentation/#generictensor"""

    def __init__(self):
        self.rhs = [0, 0, 0]
        self.exact = None
        self.dirichlet_boundary = []
        self.neumann_boundary = []

    def set_velocity(self, id, value, is_dim_fixed=None):
        """set the velocity value for the sideset id. Note the value must be a vector in 2D or 3D depending on the problem"""
        self.add_dirichlet_value(id, value, is_dim_fixed)

    def set_displacement(self, id, value, is_dim_fixed=None):
        """set the displacement value for the sideset id. Note the value must be a vector in 2D or 3D depending on the problem"""
        self.add_dirichlet_value(id, value, is_dim_fixed)

    def set_force(self, id, value):
        """set the force value for the sideset id. Note the value must be a vector in 2D or 3D depending on the problem"""
        self.add_neumann_value(id, value)

    def add_dirichlet_value(self, id, value, is_dirichlet_dim=None):
        """add the Dirichlet value value for the sideset id. Note the value must be a vector in 2D or 3D depending on the problem. is_dirichlet_dim is a vector of boolean specifying which dimentions are fixed."""
        assert(len(value) == 3 or len(value) == 2)
        tmp = {}
        tmp["id"] = id
        tmp["value"] = value
        if is_dirichlet_dim is not None:
            assert(len(is_dirichlet_dim) == 3 or len(is_dirichlet_dim) == 2)
            tmp["dimension"] = is_dirichlet_dim

        self.dirichlet_boundary.append(tmp)

    def add_neumann_value(self, id, value):
        """add the Neumann value value for the sideset id. Note the value must be a vector in 2D or 3D depending on the problem"""

        tmp = {}
        tmp["id"] = id
        tmp["value"] = value
        self.neumann_boundary.append(tmp)

    def name(self):
        return "GenericTensor"

    def params(self):
        return self.__dict__


#################################################
############Stokes###############################
class Flow:
    """Inflow/outflow problem for fluids. You can specify the sideset for the moving fluxes, the axial direction of the flow, and the list of obstacle sidesets. https://polyfem.github.io/documentation/#flow"""

    def __init__(self, inflow=1, outflow=3, inflow_amout=0.25, outflow_amout=0.25, direction=0, obstacle=[7]):
        self.inflow = inflow
        self.outflow = outflow
        self.inflow_amout = inflow_amout
        self.outflow_amout = outflow_amout
        self.direction = direction
        self.obstacle = obstacle

    def name(self):
        return "Flow"

    def params(self):
        return self.__dict__


class DrivenCavity:
    """Classical driven cavity problem in fluid simulation"""

    def name(self):
        return "DrivenCavity"

    def params(self):
        return {}


class FlowWithObstacle:
    """FLuid Obstacle problem"""

    def __init__(self, U=1.5, time_dependent=True):
        self.U = U
        self.time_dependent = time_dependent

    def name(self):
        return "FlowWithObstacle"

    def params(self):
        return self.__dict__
