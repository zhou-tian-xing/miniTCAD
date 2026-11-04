"""偏压扫描：W4 只有固定步进；自适应与回滚是 W5-B2（第 9.2 节）。

固定步进的失败语义是**立即停止**：下一个偏压仍用最后一个已收敛解热启动，
不猜测中间状态。这样 `RunResult.successful` 里永远是可以交给 writer 的真实解，
未尝试的偏压不会出现在 `attempts` 中（§11.3 不允许把失败掩盖成成功）。
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

from minitcad.core.device import NormalizedDevice, SolverConfig
from minitcad.core.state import State
from minitcad.simulation.dc import solve_single_bias
from minitcad.simulation.run_result import STATUS_FAIL, STATUS_PASS, Attempt

Progress = Callable[[float, str], None]
ShouldStop = Callable[[], bool]


@dataclass(frozen=True, slots=True)
class SweepOutcome:
    """扫描结果：`attempts` 按偏压升序，`cancelled` 表示被 `should_stop` 中断。"""

    attempts: tuple[Attempt, ...]
    cancelled: bool = False

    @property
    def last_converged_state(self) -> State | None:
        for attempt in reversed(self.attempts):
            if attempt.result is not None and attempt.result.converged:
                return attempt.result.state
        return None


def run_fixed_sweep(
    device: NormalizedDevice,
    biases: tuple[float, ...],
    initial_state: State,
    solver: SolverConfig,
    *,
    progress: Progress | None = None,
    should_stop: ShouldStop | None = None,
) -> SweepOutcome:
    """按给定偏压表逐点求解；第一个失败点后停止。

    `progress(bias_V, status)` 只在**尝试之后**回调，`status` 取 `PASS`/`FAIL`，
    调用方不能把它当成完成百分比（§10.2）。
    """
    stop = should_stop or (lambda: False)
    attempts: list[Attempt] = []
    state = initial_state

    for bias_V in biases:
        if stop():
            return SweepOutcome(attempts=tuple(attempts), cancelled=True)
        result, runtime = solve_single_bias(device, bias_V, state, solver)
        status = STATUS_PASS if result.converged else STATUS_FAIL
        attempts.append(
            Attempt(bias_V=float(bias_V), status=status, runtime_s=float(runtime), result=result)
        )
        if progress is not None:
            progress(float(bias_V), status)
        if not result.converged:
            break
        state = result.state

    return SweepOutcome(attempts=tuple(attempts), cancelled=False)


def equilibrium_start_state(device: NormalizedDevice) -> State:
    """0 V 平衡初值（A 的标量 Poisson 解），给扫描第一个偏压做热启动。

    该标量方程假设 `np = 1`，只在 0 V 成立；扫描起点不是 0 V 时仍用它做初值，
    由第一个外层迭代把状态带到目标偏压，这比用中性猜测更接近解。
    """
    from minitcad.simulation.equilibrium import solve_equilibrium_scalar

    targets = device.contact_targets(0.0)
    solution = solve_equilibrium_scalar(
        device.mesh.x_nodes_nd, device.cbar_nd, targets[0], targets[3]
    )
    if not solution.converged or solution.state is None:
        raise ArithmeticError(
            f"EQUILIBRIUM_FAILED: scalar Poisson did not converge ({solution.reason})"
        )
    return solution.state
