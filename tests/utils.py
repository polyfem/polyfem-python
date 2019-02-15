import plotly.offline as plotly
import plotly.graph_objs as go

import numpy as np


def plot(vertices, connectivity, function):
    x = vertices[:,0]
    y = vertices[:,1]

    if vertices.shape[1] == 3:
        z = vertices[:,2]
    else:
        z = np.zeros(x.shape, dtype=x.dtype)

    if connectivity.shape[1] == 3:
        f = connectivity
    else:
        f = np.ndarray([len(connectivity)*4, 3], dtype=np.int64)
        for i in range(len(connectivity)):
            f[i*4+0] = np.array([connectivity[i][1], connectivity[i][0], connectivity[i][2]])
            f[i*4+1] = np.array([connectivity[i][0], connectivity[i][1], connectivity[i][3]])
            f[i*4+2] = np.array([connectivity[i][1], connectivity[i][2], connectivity[i][3]])
            f[i*4+3] = np.array([connectivity[i][2], connectivity[i][0], connectivity[i][3]])

    mesh = go.Mesh3d(x=x, y=y, z=z,
                     i=f[:,0], j=f[:,1], k=f[:,2],
                     intensity=function, flatshading=function is not None,
                     colorscale='Viridis',
                     contour=dict(show=True, color='#fff'),
                     hoverinfo="all")
    layout = go.Layout(scene=go.layout.Scene(aspectmode='data'))
    fig = go.Figure(data=[mesh], layout=layout)

    plotly.plot(fig)
