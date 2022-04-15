import polyfempy as pf
import json
import numpy as np


root = "/Users/teseo/Downloads"
with open(root + "/pushbox.json", "r") as f:
    config = json.load(f)

config["root_path"] = root + "/pushbox.json"

dt = config["dt"]
t0 = config["t0"]


solver = pf.Solver()
solver.set_settings(json.dumps(config))
solver.load_mesh_from_settings()

# inits stuff
solver.init_timestepping(t0, dt)

for i in range(1, 6):
    # substepping
    for t in range(1):
        solver.step_in_time(t0, dt, t+1)
    solver.update_dirichlet_boundary(1, [
        "0.5*t",
        "0",
        "0"
    ])
    pixles = solver.render(width=84,
                           height=84,
                           camera_position=np.array([[0., 0., 0.]]),
                           camera_fov=90,
                           ambient_light=np.array([[0., 0., 0.]]),
                           diffuse_color=np.array([[1., 0., 0.]]),
                           specular_color=np.array([[1., 0., 0.]]))
    break

    print(pixles)
    t0 += dt


solver.export_vtu("xxx.vtu")
