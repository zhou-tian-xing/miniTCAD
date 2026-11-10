"""单偏压求解编排（第 9.2 节）。

`solve_single_bias` 是 `sweep` 与 `run` 的唯一调用点：把求解器设置、初值校验和计时
收在一处，保证 GUI/CLI 走同一条路径（§10.2）。

**方法选择**：只按 `solver.method` 选入口，两种方法的其余参数逐个透传。这里刻意不做
“Newton 失败就退回 Gummel”的兜底——那样结果里的 `solver` 标签会与实际执行的方法不符，
Gate 2 的对照实验也会把两条路径混在一起。失败就是失败，由 `sweep` 决定回滚。

失败语义：求解失败时 `SolveResult.state is None`（求解器保证），本层不抛异常，
由 `sweep` 决定停止还是回滚；只有**参数非法**才抛 `ValueError`。
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

    `solver.method` 必须是 `SOLVER_METHODS` 之一；配置层已经校验过，这里的失败是
    “有人绕过配置直接构造 `SolverConfig`”，属于编程错误，用 `ValueError` 明确暴露。
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
