# compmathproject
Репозиторий для солвера одномерного уравнения адвекции.

## Описание
Проект решает **одномерное уравнение адвекции**:

$$ 
\frac{\partial u}{\partial t} + c \frac{\partial u}{\partial x} = 0, $$

где $  u(x, t)$ — функция, (c) — скорость переноса, $ x \in [0, L] , ~ t \in [0, T]  $.

Задача с условиями Дирихле:
- Начальное: $ u(x, 0) = u_0(x) ,$
- Граничные: $ u(0, t) = u_L(t),~  u(L, t) = u_R(t).$

Поддерживаются периодические условия $ u(0, t) = u(L, t) $.

### Особенности
- Методы: Upwind, Lax-Wendroff, TVD (Superbee).
- API на **FastAPI**, запуск через **Uvicorn**.
- Визуализация с **[Plotly](https://plotly.com/python/animations/)** (3D-поверхности).
- Тесты в папке `tests`.

## Установка
1. Склонируйте:

   git clone https://github.com/furs-aka-beast/compmathproject
   
   cd compmathproject
2. Запустите сервер:

    uvicorn main:app --reload
3. Запустите тесты:

    python tests/advection_tests.py

## Выполнили
-Лавыгин Кирилл (213)

-Пермяков Максим (205)

Зоны ответственности:
Лавыгин Кирилл: FastAPI, методы Upwind, Lax-Wendroff.

Пермяков Максим: TVD-солвер.
## Выводы:
Реализованы методы для уравнения адвекции:
- Upwind: Прост, но сглаживает решение.
- Lax-Wendroff: Прост, но даёт осцилляции.
- TVD (Superbee): Cбалансирован, минимизирует диффузию и осцилляции.

TVD — оптимальный выбор для разрывных и гладких решений.