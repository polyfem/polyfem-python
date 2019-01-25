

#################################################
############Scalar###############################
class Franke:
	def name(self):
		return "Franke"

	def params(self):
		return {}


class GenericScalar:
	def __init__(self):
		self.rhs = 0
		self.dirichlet_boundary = []
		self.neumann_boundary = []

	def add_dirichlet_value(self, id, value):
		tmp = {}
		tmp["id"] = id
		tmp["value"] = value
		self.dirichlet_boundary.append(tmp)

	def add_neumann_value(self, id, value):
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
	def name(self):
		return "Gravity"

	def params(self):
		return {}


class Torsion:
	def __init__(self):
		self.axis_coordiante = 2
		self.n_turns = 0.5
		self.fixed_boundary = 5
		self.turning_boundary = 6


	def name(self):
		return "TorsionElastic"


	def params(self):
		return self.__dict__


class GenericTensor:
	def __init__(self):
		self.rhs = [0, 0, 0]
		self.dirichlet_boundary = []
		self.neumann_boundary = []

	def add_dirichlet_value(self, id, value, is_dirichlet_dim):
		assert(len(value) == 3)
		tmp = {}
		tmp["id"] = id
		tmp["value"] = value
		if is_dirichlet_dim is not None:
			assert(len(is_dirichlet_dim) == 3)
			tmp["dimension"] = is_dirichlet_dim

		self.dirichlet_boundary.append(tmp)

	def add_neumann_value(self, id, value):
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
	def __init__(self):
		self.inflow = 1
		self.outflow = 3
		self.inflow_amout = 0.25
		self.outflow_amout = 0.25
		self.obstacle = [7]

	def name(self):
		return "Flow"

	def params(self):
		return self.__dict__


class DrivenCavity:
	def name(self):
		return "DrivenCavity"

	def params(self):
		return {}
