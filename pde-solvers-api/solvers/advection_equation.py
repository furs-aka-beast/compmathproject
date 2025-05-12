import numpy as np

def solve_advection_equation(c: float, length: float, nx: int, nt: int, dt: float,
                            initial_condition: np.ndarray, bc: dict) -> np.ndarray:
    """
    Решает одномерное уравнение адвекции методом "вверх по потоку".
    Args:
        c: Скорость переноса
        length: Длина области
        nx: Количество узлов по пространству
        nt: Количество временных шагов
        dt: Шаг по времени
        initial_condition: Начальное условие u(x,0)
        bc: Граничные условия, например, {"left": 0.0}
    Returns:
        Матрица решения u(x,t) в виде NumPy массива.
    """
    dx = length / (nx - 1)
    u = np.zeros((nx, nt))
    u[:, 0] = initial_condition  # Задаем начальное условие

    # Проверка условия Куранта
    courant = c * dt / dx
    if courant > 1:
        raise ValueError(f"Нарушено условие Куранта: c * dt / dx = {courant} > 1")

    # Метод "вверх по потоку" (предполагаем c > 0)
    for t in range(1, nt):
        u[0, t] = bc.get("left", 0.0)  # Граничное условие слева
        for x in range(1, nx):
            u[x, t] = u[x, t-1] - courant * (u[x, t-1] - u[x-1, t-1])
        # Если нужны периодические граничные условия:
        # u[0, t] = u[-1, t-1]

    return u

def solve_advection_lax_wendroff(c: float, length: float, nx: int, nt: int, dt: float,
                                 initial_condition: np.ndarray, bc: dict) -> np.ndarray:
    dx = length / (nx - 1)
    u = np.zeros((nx, nt))
    u[:, 0] = initial_condition
    courant = c * dt / dx

    for t in range(1, nt):
        u[0, t] = bc.get("left", 0.0)
        for x in range(1, nx-1):
            u[x, t] = u[x, t-1] - 0.5 * courant * (u[x+1, t-1] - u[x-1, t-1]) + \
                      0.5 * courant**2 * (u[x+1, t-1] - 2*u[x, t-1] + u[x-1, t-1])
        u[-1, t] = bc.get("right", u[-1, t-1])  # Правая граница
    return u