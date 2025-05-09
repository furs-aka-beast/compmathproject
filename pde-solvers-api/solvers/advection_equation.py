import numpy as np

def solve_advection_equation(a: float, length: float, nx: int, nt: int, dt: float,
initial_condition: np.ndarray, bc: float) -> np.ndarray:
    dx = length / (nx - 1)
    u = np.zeros((nx, nt))
    u[:, 0] = initial_condition
    cfl = a * dt / dx
    if cfl > 1.0:
        raise ValueError("CFL condition violated: a * dt / dx must be ≤ 1")

    for t in range(1, nt):
        u[0, t] = bc  # граничное условие Дирихле
        for x in range(1, nx):
            u[x, t] = u[x, t-1] - cfl * (u[x, t-1] - u[x-1, t-1])
    return u