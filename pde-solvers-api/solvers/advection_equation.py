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
        bc: Граничные условия, например, {"left": 0.0} или {} для периодических
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
        u[0, t] = u[-1, t-1] # Периодические граничные условия:

    return u

def solve_advection_lax_wendroff(c: float, length: float, nx: int, nt: int, dt: float,
                                 initial_condition: np.ndarray, bc: dict) -> np.ndarray:
    """
    Решает одномерное уравнение адвекции методом Лакса-Вендроффа.
    Args:
        c: Скорость переноса
        length: Длина области
        nx: Количество узлов по пространству
        nt: Количество временных шагов
        dt: Шаг по времени
        initial_condition: Начальное условие u(x,0)
        bc: Граничные условия, например, {"left": 0.0} или {} для периодических
    Returns:
        Матрица решения u(x,t) в виде NumPy массива.
    """
    dx = length / (nx - 1)
    u = np.zeros((nx, nt))
    u[:, 0] = initial_condition
    courant = c * dt / dx
    if abs(courant) > 1:
        raise ValueError(f"Нарушено условие Куранта: |c * dt / dx| = {abs(courant)} > 1")

    for t in range(1, nt):
        if not bc: # Периодические граничные условия
            u[0, t] = u[-1, t-1]
            u[-1, t] = u[0, t-1]
        else:
            u[0, t] = bc.get("left", 0.0)
            u[-1, t] = bc.get("right", 0.0)

        # Схема Лакса-Вендроффа
        for x in range(1, nx-1):
            u[x, t] = u[x, t-1] - 0.5 * courant * (u[x+1, t-1] - u[x-1, t-1]) + \
                      0.5 * courant**2 * (u[x+1, t-1] - 2*u[x, t-1] + u[x-1, t-1])
    
    return u

def superbee_limiter(r):
    """Ограничитель Superbee для скаляра или массива"""
    if isinstance(r, (np.ndarray, list)):
        # Для массивов используем NumPy
        return np.maximum(0, np.minimum(2 * r, 1), np.minimum(r, 2))
    else:
        # Для скаляров используем Python min/max
        return max(0, min(2 * r, 1), min(r, 2))

def solve_advection_tvd(c: float, length: float, nx: int, nt: int, dt: float,
                        initial_condition: np.ndarray, bc: dict) -> np.ndarray:
    """
    Решает одномерное уравнение адвекции методом TVD с ограничителем Superbee.
    Args:
        c: Скорость переноса
        length: Длина области
        nx: Количество узлов по пространству
        nt: Количество временных шагов
        dt: Шаг по времени
        initial_condition: Начальное условие u(x,0)
        bc: Граничные условия, например, {"left": 0.0} или {} для периодических
    Returns:
        Матрица решения u(x,t) в виде NumPy массива.
    """
    dx = length / (nx - 1)
    u = np.zeros((nx, nt))
    u[:, 0] = initial_condition
    courant = c * dt / dx
    if abs(courant) > 1:
        raise ValueError(f"Нарушено условие Куранта: |c * dt / dx| = {abs(courant)} > 1")

    for t in range(1, nt):
        # Периодические граничные условия
        if not bc:
            u[0, t] = u[-1, t-1]
            u[-1, t] = u[0, t-1]
        else:
            u[0, t] = bc.get("left", 0.0)
            u[-1, t] = bc.get("right", 0.0)

        # Вычисление потоков
        F = np.zeros(nx)
        if c > 0:
            for i in range(1, nx-1):
                r = (u[i, t-1] - u[i-1, t-1]) / (u[i+1, t-1] - u[i, t-1] + 1e-10)  # Избегаем деления на 0
                phi = superbee_limiter(r)
                F[i] = c * u[i, t-1] + 0.5 * c * (1 - courant) * phi * (u[i+1, t-1] - u[i, t-1])
            # Обновление решения
            for i in range(1, nx-1):
                u[i, t] = u[i, t-1] - (dt / dx) * (F[i] - F[i-1])
        else:
            for i in range(1, nx-1):
                r = (u[i+1, t-1] - u[i, t-1]) / (u[i, t-1] - u[i-1, t-1] + 1e-10)
                phi = superbee_limiter(r)
                F[i] = c * u[i, t-1] - 0.5 * c * (1 + courant) * phi * (u[i, t-1] - u[i-1, t-1])
            for i in range(1, nx-1):
                u[i, t] = u[i, t-1] - (dt / dx) * (F[i+1] - F[i])

    return u