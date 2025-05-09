from pydantic import BaseModel
from typing import List, Literal

class AdvectionEquationInput(BaseModel):
    a: float                         # скорость переноса
    length: float                    # длина области
    nx: int                          # число узлов по x
    nt: int                          # число шагов по времени
    dt: float                        # шаг по времени
    initial_condition: List[float]   # u(x, 0)
    boundary_condition: float        # значение u на левом краю (x=0)