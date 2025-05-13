from fastapi import FastAPI
from fastapi.responses import JSONResponse
import numpy as np
from schemas import AdvectionEquationInput
from solvers.advection_equation import solve_advection_equation, solve_advection_lax_wendroff, solve_advection_tvd
import plotly.graph_objects as go

app = FastAPI(title="Advection Solver API")

@app.post("/solve/advection-equation/")
async def solve_advection_eq(params: AdvectionEquationInput):
    u0 = np.array(params.initial_condition)
    
    # Выбор метода
    if params.method == "tvd":
        solution = solve_advection_tvd(
            c=params.c, length=params.length, nx=params.nx, nt=params.nt,
            dt=params.dt, initial_condition=u0, bc=params.boundary_conditions
        )
    elif params.method == "lax_wendroff":
        solution = solve_advection_lax_wendroff(
            c=params.c, length=params.length, nx=params.nx, nt=params.nt,
            dt=params.dt, initial_condition=u0, bc=params.boundary_conditions
        )
    else:
        solution = solve_advection_equation(
            c=params.c, length=params.length, nx=params.nx, nt=params.nt,
            dt=params.dt, initial_condition=u0, bc=params.boundary_conditions
        )

    # Визуализация
    fig = go.Figure(data=[go.Surface(z=solution.T)])
    fig.update_layout(title=f"Advection Equation Solution ({params.method.upper()})")
    plot_html = fig.to_html(full_html=False)

    return {
        "solution": solution.tolist(),
        "visualization": plot_html
    }