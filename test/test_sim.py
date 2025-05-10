import polyfempy as pf
import json
import numpy as np

pf.polyfem_command(json="../data/contact/examples/2D/unit-tests/5-squares.json", log_level=2, max_threads=16)

root = "../data/contact/examples/2D/unit-tests"
with open(root + "/5-squares.json", "r") as f:
    config = json.load(f)

config["root_path"] = root + "/5-squares.json"

solver = pf.Solver()
solver.set_settings(json.dumps(config), False)
solver.set_log_level(2)
solver.load_mesh_from_settings()

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
    solver.export_vtu(sol, np.zeros((0, 0)), t0, dt, "step_" + str(i) + ".vtu")
