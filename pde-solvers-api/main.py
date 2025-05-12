from fastapi import FastAPI
from fastapi.responses import JSONResponse
import numpy as np
from schemas import AdvectionEquationInput
from solvers.advection_equation import solve_advection_equation
import plotly.graph_objects as go

app = FastAPI(title="Advection Solver API")

@app.post("/solve/advection-equation/")
async def solve_advection_eq(params: AdvectionEquationInput):
    # Конвертируем входные данные в NumPy
    u0 = np.array(params.initial_condition)
    
    # Решаем уравнение адвекции
    solution = solve_advection_equation(
        c=params.c,
        length=params.length,
        nx=params.nx,
        nt=params.nt,
        dt=params.dt,
        initial_condition=u0,
        bc=params.boundary_conditions
    )

    # Генерируем визуализацию с помощью Plotly
    fig = go.Figure(data=[go.Surface(z=solution.T)])
    fig.update_layout(title="Advection Equation Solution")
    plot_html = fig.to_html(full_html=False)

    # Возвращаем решение и визуализацию
    return {
        "solution": solution.tolist(),
        "visualization": plot_html
    }