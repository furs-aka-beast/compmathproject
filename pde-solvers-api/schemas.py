from pydantic import BaseModel
from typing import List

class AdvectionEquationInput(BaseModel):
    c: float               # Скорость переноса (латинская 'c')
    length: float          # Длина области
    nx: int                # Количество узлов по пространству
    nt: int                # Количество временных шагов
    dt: float              # Шаг по времени
    initial_condition: List[float]  # Начальное условие u(x,0)
    boundary_conditions: dict       # Граничные условия, например, {"left": 0.0}