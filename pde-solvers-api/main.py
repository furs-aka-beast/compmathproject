from fastapi import FastAPI
from fastapi.responses import JSONResponse
import numpy as np
from schemas import AdvectionEquationInput
from solvers.advection_equation import solve_advection_equation
import plotly.graph_objects as go
import base64

app = FastAPI(title="PDE Solver API")

@app.post("/solve/heat-equation/")
async def solve_advection(params: AdvectionEquationInput): # asynchronous call to computing server
    # Convert input to NumPy
    u0 = np.array(params.initial_condition)
    
    # Solve PDE
    solution = solve_advection_equation(
        a=params.a,
        length=params.length,
        nx=params.nx,
        nt=params.nt,
        dt=params.dt,
        initial_condition=u0,
        bc=params.boundary_condition
    )

    fig = go.Figure(data=[go.Surface(z=solution.T)])
    fig.update_layout(title="Advection Equation (Upwind Scheme)")
    plot_html = fig.to_html(full_html=False)

    return {
        "solution": solution.tolist(),
        "visualization": plot_html
    }