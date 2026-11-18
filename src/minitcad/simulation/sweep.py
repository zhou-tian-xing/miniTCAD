"""偏压扫描：固定步进与自适应继续法 + 回滚（§9.2）。

两种模式的失败语义不同：

* 固定步进（`run_fixed_sweep`）：失败即停，下一个偏压仍用最后一个已收敛解热启动。
* 自适应（`run_adaptive_sweep`）：失败就减小步长，从最后一个已接受解重新逼近；
  失败的那一次仍写入 `attempts` 并标 `FAIL`。

“重试同一个目标偏压”只能理解为：总体目标 `stop_V` 不变，把失败的一跳拆成更小的跳。
求解器是确定性的，用同样的阶跃从同一个状态重试必然再失败。

为什么需要它：Newton 在固定 0.05 V 步长下会在 0.10 V 失败，把步长降到 0.02 V 以内就能
一路收敛到 0.30 V（`reports/week05_acceptance.md` §4.1）。
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

from minitcad.core.device import NormalizedDevice, SolverConfig, SweepConfig
from minitcad.core.state import SolveResult, State
from minitcad.simulation.dc import solve_single_bias
from minitcad.simulation.run_result import STATUS_FAIL, STATUS_PASS, Attempt

Progress = Callable[[float, str], None]
ShouldStop = Callable[[], bool]
Solve = Callable[[NormalizedDevice, float, State, SolverConfig], tuple[SolveResult, float]]


@dataclass(frozen=True, slots=True)
class SweepOutcome:
    """扫描结果：`attempts` 按尝试顺序，`cancelled` 表示被 `should_stop` 中断。

    自适应回滚会让 `attempts` 里出现同一段偏压的多次尝试，所以它不是偏压表。
    """

    attempts: tuple[Attempt, ...]
    cancelled: bool = False

    @property
    def last_converged_state(self) -> State | None:
        for attempt in reversed(self.attempts):
            if attempt.result is not None and attempt.result.converged:
                return attempt.result.state
        return None

    @property
    def successful(self) -> tuple[Attempt, ...]:
        return tuple(attempt for attempt in self.attempts if attempt.status == STATUS_PASS)

    @property
    def reached_V(self) -> float | None:
        """最后一个收敛点的偏压；一个都没收敛时返回 `None`。"""
        passed = self.successful
        return passed[-1].bias_V if passed else None


def run_fixed_sweep(
    device: NormalizedDevice,
    biases: tuple[float, ...],
    initial_state: State,
    solver: SolverConfig,
    *,
    progress: Progress | None = None,
    should_stop: ShouldStop | None = None,
    solve: Solve = solve_single_bias,
) -> SweepOutcome:
    """按给定偏压表逐点求解；第一个失败点后停止。

    `progress(bias_V, status)` 在每次尝试之后回调，`status` 只有 `PASS`/`FAIL`。
    """
    stop = should_stop or (lambda: False)
    attempts: list[Attempt] = []
    state = initial_state

    for bias_V in biases:
        if stop():
            return SweepOutcome(attempts=tuple(attempts), cancelled=True)
        result, runtime = solve(device, bias_V, state, solver)
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


def run_adaptive_sweep(
    device: NormalizedDevice,
    sweep: SweepConfig,
    initial_state: State,
    solver: SolverConfig,
    *,
    progress: Progress | None = None,
    should_stop: ShouldStop | None = None,
    solve: Solve = solve_single_bias,
) -> SweepOutcome:
    """自适应步长 + 回滚地走到 `sweep.stop_V`。

    * 成功：`step = min(step * grow_factor, max_step_V)`，下一个目标 `= 已接受偏压 + step`；
    * 失败：记 `FAIL`，`step = max(step * shrink_factor, min_step_V)`，从已接受解重试；
    * 步长已经贴着 `min_step_V` 仍失败：返回，由 `run()` 判 `REJECTED`。

    乘性而非加性：`eta` 随偏压线性增长、通量随 `exp` 增长，收敛半径是指数级的。

    `solve` 可注入，供 T13 伪造一次求解失败；生产路径上是 `solve_single_bias`。
    """
    stop = should_stop or (lambda: False)
    attempts: list[Attempt] = []
    state = initial_state
    step = sweep.step_V
    accepted_V = sweep.start_V

    # 起点也是一次尝试：通常是 0 V 平衡解，但要在 attempts 里留痕，否则 CSV 缺第一个偏压。
    result, runtime = solve(device, accepted_V, state, solver)
    status = STATUS_PASS if result.converged else STATUS_FAIL
    attempts.append(
        Attempt(bias_V=float(accepted_V), status=status, runtime_s=float(runtime), result=result)
    )
    if progress is not None:
        progress(float(accepted_V), status)
    if not result.converged:
        return SweepOutcome(attempts=tuple(attempts), cancelled=False)
    state = result.state

    while accepted_V < sweep.stop_V:
        if stop():
            return SweepOutcome(attempts=tuple(attempts), cancelled=True)

        # 取 12 位小数，与 `SweepConfig.biases()` 同口径，避免累积出 0.08750000000000001。
        target_V = round(min(accepted_V + step, sweep.stop_V), 12)
        result, runtime = solve(device, target_V, state, solver)
        if result.converged:
            status = STATUS_PASS
            note = None
        else:
            status = STATUS_FAIL
            note = (
                f"rollback: step {step:.6g} -> {max(step * sweep.shrink_factor, sweep.min_step_V):.6g}"
            )
        attempts.append(
            Attempt(
                bias_V=float(target_V),
                status=status,
                runtime_s=float(runtime),
                result=result,
                note=note,
            )
        )
        if progress is not None:
            progress(float(target_V), status)

        if result.converged:
            state = result.state
            accepted_V = target_V
            step = min(step * sweep.grow_factor, sweep.max_step_V)
            continue

        shrunk = max(step * sweep.shrink_factor, sweep.min_step_V)
        if shrunk >= step:
            # 贴着 min_step_V 了：再缩是同一个步长，重试只会重复失败。
            return SweepOutcome(attempts=tuple(attempts), cancelled=False)
        step = shrunk

    return SweepOutcome(attempts=tuple(attempts), cancelled=False)


def equilibrium_start_state(device: NormalizedDevice) -> State:
    """0 V 平衡初值（A 的标量 Poisson 解），给第一个偏压做热启动。

    该标量方程假设 `np = 1`，只在 0 V 成立；起点不是 0 V 时仍拿它当初值，
    由第一个外层迭代带到目标偏压，比中性猜测更接近解。
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
