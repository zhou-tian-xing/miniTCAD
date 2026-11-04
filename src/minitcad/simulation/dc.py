"""单偏压求解编排（第 9.2 节）。

`solve_single_bias` 是 `sweep` 与 `run` 的唯一调用点：把求解器设置、初值校验和计时
收在一处，保证 GUI/CLI 走同一条路径（§10.2）。

失败语义：求解失败时 `SolveResult.state is None`（A 的实现保证），本层不抛异常，
由 `sweep` 决定停止还是回滚；只有**参数非法**才抛 `ValueError`。
"""

from __future__ import annotations

from time import perf_counter

from minitcad.core.device import NormalizedDevice, SolverConfig
from minitcad.core.state import SolveResult, State
from minitcad.solvers import solve_gummel


def solve_single_bias(
    device: NormalizedDevice,
    bias_V: float,
    initial_state: State,
    solver: SolverConfig,
) -> tuple[SolveResult, float]:
    """解一个偏压点，返回 `(结果, 墙钟时间[s])`。

    `solver.method == "newton"` 时直接拒绝：A 的 Newton 在 W5 落地，
    悄悄退回 Gummel 会让结果里的方法标签与实际不符。
    """
    if solver.method == "newton":
        raise NotImplementedError("solve_newton lands in W5 (A); use solver.method='gummel'")
    started = perf_counter()
    result = solve_gummel(
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
