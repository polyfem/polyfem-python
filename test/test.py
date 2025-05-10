# %%
import polyfempy as pf
import json
import numpy as np

# %%
root = "../data/differentiable/input"
with open(root + "/initial-contact.json", "r") as f:
    config = json.load(f)

config["contact"]["use_convergent_formulation"] = True
config["root_path"] = root + "/initial-contact.json"

solver = pf.Solver()
solver.set_settings(json.dumps(config), True)
solver.set_log_level(2)
solver.load_mesh_from_settings()

# %%
mesh = solver.mesh()

print(mesh.n_vertices())
print(mesh.n_elements())
print(mesh.n_cell_vertices(1))
print(mesh.element_vertex(3, 0))
print(mesh.boundary_element_vertex(3, 0))
assert(mesh.is_boundary_vertex(1))

min, max = mesh.bounding_box()

# %%
config = solver.settings()
t0 = config["time"]["t0"]
dt = config["time"]["dt"]

# inits stuff
solver.build_basis()
solver.assemble()
sol = solver.init_timestepping(t0, dt)

for i in range(1, 5):
    
    # substepping
    for t in range(1):
        sol = solver.step_in_time(sol, t0, dt, t+1)

    t0 += dt
    solver.export_vtu("step_" + str(i) + ".vtu", sol, np.zeros((0, 0)), t0, dt)

prob = solver.nl_problem()

h = prob.hessian(sol)
reduced_sol = prob.full_to_reduced(sol)
full_sol = prob.reduced_to_full(reduced_sol)

assert(np.linalg.norm(full_sol - sol.flatten()) < 1e-12)

cache = solver.get_solution_cache()

print(cache.solution(1).shape)
print(cache.velocity(2).shape)
print(cache.acceleration(3).shape)
print(cache.hessian(4).shape)
