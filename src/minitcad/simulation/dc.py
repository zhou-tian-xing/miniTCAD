"""单偏压求解编排（第 9.2 节）：`sweep` 与 `run` 的唯一调用点。

只按 `solver.method` 选入口，不写“Newton 失败就退回 Gummel”的兜底——那样结果里的
`solver` 标签与实际执行的方法不符。求解失败时 `SolveResult.state is None`，本层不抛
异常，由 `sweep` 决定停止还是回滚；只有参数非法才抛 `ValueError`。
"""

from __future__ import annotations

from time import perf_counter

from minitcad.core.device import NormalizedDevice, SolverConfig
from minitcad.core.state import SolveResult, State
from minitcad.solvers import solve_gummel, solve_newton

_SOLVERS = {"gummel": solve_gummel, "newton": solve_newton}


def solve_single_bias(
    device: NormalizedDevice,
    bias_V: float,
    initial_state: State,
    solver: SolverConfig,
) -> tuple[SolveResult, float]:
    """解一个偏压点，返回 `(结果, 墙钟时间[s])`。

    `solver.method` 已在配置层校验；这里的 `ValueError` 对应“绕过配置直接构造
    `SolverConfig`”，属编程错误。
    """
    try:
        solve = _SOLVERS[solver.method]
    except KeyError as exc:
        raise ValueError(f"unknown solver method: {solver.method!r}") from exc
    started = perf_counter()
    result = solve(
        device,
        bias_V,
        initial_state,
        tol_residual=solver.tol_residual,
        tol_update=solver.tol_update,
        max_iter=solver.max_iter,
        alpha_min=solver.alpha_min,
        n_floor=solver.n_floor,
    )
    return result, perf_counter() - started
