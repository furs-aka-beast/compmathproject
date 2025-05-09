import numpy as np
import plotly.graph_objects as go
from solvers.advection_equation import solve_advection_equation

length = 1.0
nx = 100
nt = 200
dt = 0.005
a = 1.0

x = np.linspace(0, length, nx)
initial = np.where((x > 0.3) & (x < 0.5), 1.0, 0.0)

u = solve_advection_equation(a, length, nx, nt, dt, initial, bc=0.0)

fig = go.Figure(data=[go.Surface(z=u.T, x=x, y=np.arange(nt) * dt)])
fig.update_layout(title="Upwind Scheme: Square Pulse")
fig.write_html("tests/advection_plot.html")