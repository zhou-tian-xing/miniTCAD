"""T13：自适应继续法与回滚（§9.2、W5-B2）。

被测对象是 `run_adaptive_sweep` 的步长状态机，故障用 `solve` 钩子注入——注入点伪造一次
“求解失败”，不必真的构造一个不收敛的偏压，因此这个测试不依赖求解器的具体收敛行为。

四条被钉住的语义：

1. 失败只影响步长，不影响已接受的解：下一次尝试从**最后一个已接受解**出发；
2. 失败尝试留在 `attempts` 里，并带 `rollback:` 说明（§11.3 不允许把失败掩盖成成功）；
3. 缩步后仍能到达目标偏压，整体判决可以是 ACCEPTED；
4. 缩到 `min_step_V` 还失败就放弃，不会无限重试。
"""

from __future__ import annotations

import dataclasses
import pathlib

import numpy as np
import pytest

from minitcad.api import run
from minitcad.core.device import build_device
from minitcad.core.state import SolveResult
from minitcad.io.config import load_config
from minitcad.simulation.dc import solve_single_bias
from minitcad.simulation.sweep import (
    equilibrium_start_state,
    run_adaptive_sweep,
    run_fixed_sweep,
)

EXAMPLES = pathlib.Path(__file__).resolve().parents[1] / "examples"
BASELINE = EXAMPLES / "pn_baseline_0_to_030.yaml"
NEWTON_BASELINE = EXAMPLES / "pn_newton_adaptive_0_to_030.yaml"


def _adaptive_config(**overrides):
    config = load_config(BASELINE)
    sweep = dataclasses.replace(config.sweep, adaptive=True, **overrides)
    return dataclasses.replace(config, sweep=sweep)


@pytest.fixture(scope="module")
def device_and_start():
    config = load_config(BASELINE)
    device = build_device(config)
    return device, equilibrium_start_state(device)


class _FailAbove:
    """在 `bias_V >= threshold` 上伪造 `times` 次失败；`times=None` 表示一直失败。"""

    def __init__(self, threshold: float, times: int | None = 1) -> None:
        self.threshold = threshold
        self.times = times
        self.calls: list[tuple[float, object]] = []

    def __call__(self, device, bias_V, state, solver):
        self.calls.append((float(bias_V), state))
        result, runtime = solve_single_bias(device, bias_V, state, solver)
        if bias_V >= self.threshold and (self.times is None or self.times > 0):
            if self.times is not None:
                self.times -= 1
            return (
                SolveResult(
                    bias_V=float(bias_V),
                    state=None,
                    converged=False,
                    iterations=result.iterations,
                    residual_norm=np.inf,
                    update_norm=np.inf,
                    reason="INJECTED_FAILURE",
                ),
                runtime,
            )
        return result, runtime


def test_rollback_shrinks_the_step_and_still_reaches_the_target(device_and_start):
    device, start = device_and_start
    config = _adaptive_config(stop_V=0.3, step_V=0.05, min_step_V=1e-3)
    injector = _FailAbove(threshold=0.3, times=1)

    outcome = run_adaptive_sweep(
        device, config.sweep, start, config.solver, solve=injector
    )

    assert outcome.cancelled is False
    assert outcome.reached_V == pytest.approx(0.3)
    failures = [attempt for attempt in outcome.attempts if attempt.status == "FAIL"]
    assert len(failures) == 1
    assert failures[0].bias_V == pytest.approx(0.3)
    assert failures[0].note is not None and failures[0].note.startswith("rollback:")

    # 回滚之后的那一次成功，其步长必须比失败时的那一步更小。
    index = outcome.attempts.index(failures[0])
    accepted = [attempt for attempt in outcome.attempts[:index] if attempt.status == "PASS"]
    failed_step = 0.3 - accepted[-1].bias_V
    retried = outcome.attempts[index + 1]
    assert retried.status == "PASS"
    assert retried.bias_V - accepted[-1].bias_V < failed_step


def test_retry_starts_from_the_last_accepted_state(device_and_start):
    """回滚只回退步长，不回退状态：失败后的下一次求解必须用最后一次收敛的解。"""
    device, start = device_and_start
    config = _adaptive_config()
    injector = _FailAbove(threshold=0.3, times=1)
    outcome = run_adaptive_sweep(
        device, config.sweep, start, config.solver, solve=injector
    )

    index = next(i for i, a in enumerate(outcome.attempts) if a.status == "FAIL")
    last_pass = outcome.attempts[index - 1]
    assert last_pass.status == "PASS"
    state_after_failure = injector.calls[index + 1][1]
    assert state_after_failure is last_pass.result.state


def test_exhausted_step_gives_up_instead_of_retrying_forever(device_and_start):
    device, start = device_and_start
    config = _adaptive_config(min_step_V=0.01, step_V=0.05)
    injector = _FailAbove(threshold=0.29, times=None)

    outcome = run_adaptive_sweep(
        device, config.sweep, start, config.solver, solve=injector
    )

    assert outcome.reached_V is not None and outcome.reached_V < 0.29
    assert all(attempt.bias_V != pytest.approx(0.3) for attempt in outcome.successful)
    # 缩到 min_step_V 就停：失败次数是 log_{1/shrink}(step/min_step) 量级，不是无限。
    failures = [attempt for attempt in outcome.attempts if attempt.status == "FAIL"]
    assert 1 <= len(failures) <= 10


def test_fixed_sweep_still_stops_at_the_first_failure(device_and_start):
    """W4 的固定步进语义不能被自适应改写：失败即停，不留后续偏压。"""
    device, start = device_and_start
    config = load_config(BASELINE)
    injector = _FailAbove(threshold=0.2, times=1)
    outcome = run_fixed_sweep(
        device, config.sweep.biases(), start, config.solver, solve=injector
    )
    assert outcome.attempts[-1].status == "FAIL"
    assert outcome.attempts[-1].bias_V == pytest.approx(0.2)
    assert outcome.reached_V == pytest.approx(0.15)


def test_run_accepts_a_config_that_needed_rollbacks():
    """判据是“到达 `stop_V`”，不是“一次都没失败”：回滚成功的运行仍是 ACCEPTED。"""
    config = load_config(NEWTON_BASELINE)
    result = run(config)
    assert result.accepted and result.failure_reason is None
    assert any(attempt.status == "FAIL" for attempt in result.attempts)
    assert result.successful[-1].bias_V == pytest.approx(config.sweep.stop_V)
