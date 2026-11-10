"""偏压扫描：固定步进（W4）与自适应继续法 + 回滚（W5-B2，§9.2）。

两种模式的失败语义不同，这是本节最重要的一件事：

* **固定步进**（`run_fixed_sweep`）：失败即停。下一个偏压仍用最后一个已收敛解热启动，
  不猜测中间状态。`RunResult.successful` 里永远是可以交给 writer 的真实解，未尝试的偏压
  不会出现在 `attempts` 中（§11.3 不允许把失败掩盖成成功）。
* **自适应**（`run_adaptive_sweep`）：失败就把步长减半，从一个**更近的**偏压重新逼近。
  失败的那一次仍然写入 `attempts` 并标 `FAIL`——回滚不等于失败没发生。

**“重试同一个目标偏压”的正确读法**：求解器是确定性的，用同样的阶跃从同一个状态重试
只会得到同样的失败，所以“缩小步长后重试同一目标”必然被理解成：**总体目标（`stop_V`）
不变，把失败的这一跳拆成几个更小的跳**。因此回滚后下一次尝试的偏压是
`last_accepted + new_step`，而不是那个刚失败的偏压。

**为什么需要它**（不是可选的健壮性装饰）：Newton 在 0.05 V 的固定步长下会在 0.10 V
`LINE_SEARCH_FAILED`——从 0.05 V 的解直接跳到 0.10 V 会在左接触形成 `1e12` 量级的虚假
边界层通量，线搜索只能把 `alpha` 压到 `1e-3` 量级，60 次迭代也走不出这一层。把步长降到
0.02 V 以内，同一个 Newton 能一路收敛到 0.30 V（见 `reports/week05_acceptance.md`）。
步长自适应是把“求解器能力”和“继续法步长”这两件事解耦的标准做法。
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

    `attempts` 可能包含**同一个偏压的多次尝试**（自适应回滚），因此它不是偏压表；
    要取“最终每个偏压的结论”，用 `successful` 或按 `status` 过滤。
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

    `progress(bias_V, status)` 只在**尝试之后**回调，`status` 取 `PASS`/`FAIL`，
    调用方不能把它当成完成百分比（§10.2）。
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

    步长状态机：

    * 成功：`step = min(step * grow_factor, max_step_V)`，下一个偏压 `= 已接受偏压 + step`；
    * 失败：记 `FAIL` 并回滚到最后一个已接受解，`step = max(step * shrink_factor, min_step_V)`，
      下一次尝试`= 已接受偏压 + 新 step`（即把一个失败的大跳拆成若干小跳）；
    * 缩到 `min_step_V` 仍失败：无法再缩小，返回（`run()` 据此判 `REJECTED`）。

    **为什么用乘性而不是加性步长**：收敛半径随偏压变化是指数级的（`eta` 随 `V/V_T` 线性
    增长，通量随 `exp` 增长），加性步长要么在高偏压处太大、要么在低偏压处浪费一半的步数。

    `solve` 可注入，供 T13 的故障注入使用：注入点必须能伪造一次求解失败，而不必真的
    构造一个不收敛的偏压。生产路径上它就是 `solve_single_bias`。
    """
    stop = should_stop or (lambda: False)
    attempts: list[Attempt] = []
    state = initial_state
    step = sweep.step_V
    accepted_V = sweep.start_V

    # 起点本身也是一次尝试：它通常是 0 V 平衡解，一次迭代就通过，
    # 但必须在 `attempts` 里留痕，否则 CSV 的 sweep 行会缺第一个偏压（§12.1）。
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

        # 取 12 位小数与 `SweepConfig.biases()` 同一口径：否则反复 `accepted + step`
        # 会在 CSV/报告里留下 `0.08750000000000001` 这种由二进制加法累积出来的偏差。
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
            # 已经贴着 min_step_V：再缩也还是这一步，重试只会重复同一个失败。
            return SweepOutcome(attempts=tuple(attempts), cancelled=False)
        step = shrunk

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
