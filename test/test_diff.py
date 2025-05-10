import polyfempy as pf
import json
import numpy as np
import torch

torch.set_default_dtype(torch.float64)

# Differentiable simulator that computes shape derivatives
class Simulate(torch.autograd.Function):

    @staticmethod
    def forward(ctx, solvers, vertices):
        solutions = []
        for solver in solvers:
            solver.mesh().set_vertices(vertices)
            solver.set_cache_level(pf.CacheLevel.Derivatives) # enable backward derivatives
            solver.solve()
            sol = torch.tensor(solver.get_solutions())
            solutions.append(sol)
        
        ctx.save_for_backward(vertices)
        ctx.solvers = solvers
        return tuple(solutions)

    @staticmethod
    @torch.autograd.function.once_differentiable
    def backward(ctx, *grad_output):
        vertices, = ctx.saved_tensors
        grad = torch.zeros_like(vertices)
        for i, solver in enumerate(ctx.solvers):
            solver.solve_adjoint(grad_output[i])
            grad += torch.tensor(pf.shape_derivative(solver))
        return None, grad


root = "../data/differentiable/input"
with open(root + "/initial-contact.json", "r") as f:
    config = json.load(f)

config["root_path"] = root + "/initial-contact.json"

# Simulation 1

solver1 = pf.Solver()
solver1.set_settings(json.dumps(config), False)
solver1.set_log_level(2)
solver1.load_mesh_from_settings()
# solver1.solve()

mesh = solver1.mesh()
v = mesh.vertices()
vertices = torch.tensor(solver1.mesh().vertices(), requires_grad=True)

# Simulation 2

config["initial_conditions"]["velocity"][0]["value"] = [3, 0]
solver2 = pf.Solver()
solver2.set_settings(json.dumps(config), False)
solver2.set_log_level(2)
solver2.load_mesh_from_settings()
# solver2.solve()

# Verify gradient

def loss(vertices):
    solutions1, solutions2 = Simulate.apply([solver1, solver2], vertices)
    return torch.linalg.norm(solutions1[:, -1]) * torch.linalg.norm(solutions2[:, -1])

torch.set_printoptions(12)

param = vertices.clone().detach().requires_grad_(True)
theta = torch.randn_like(param)
l = loss(param)
l.backward()
grad = param.grad
t = 1e-6
with torch.no_grad():
    analytic = torch.tensordot(grad, theta)
    f1 = loss(param + theta * t)
    f2 = loss(param - theta * t)
    fd = (f1 - f2) / (2 * t)
    print(f'grad {analytic}, fd {fd} {(f1 - l) / t} {(l - f2) / t}, relative err {abs(analytic - fd) / abs(analytic):.3e}')
    print(f'f(x+dx)={f1}, f(x)={l.detach()}, f(x-dx)={f2}')
    assert(abs(analytic - fd) <= 1e-4 * abs(analytic))
