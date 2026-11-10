"""W5-B3：`run()` 的进度与取消契约（§10.2）。

取消是**协作式**的：`should_stop()` 只在下一次尝试偏压之前被检查，正在进行的求解不会被打断。
这个取舍是明确的——中途打断一个 Newton 迭代会留下一个既没收敛也没记录的状态，而
“已经算完的偏压点完整保留”对 GUI 更有用（用户可以继续画已求出的曲线）。

契约条目：

* 取消 → `decision=REJECTED`、`failure_reason="CANCELLED"`，**不是** ACCEPTED；
* 已接受的偏压点与三文件照样能写出来，且 `run_id` 一致；
* `progress` 只回报**已经尝试过**的偏压，不回报计划中的偏压；
* 自适应模式下取消同样成立（回滚与取消是两条独立的退出路径）。
"""

from __future__ import annotations

import dataclasses
import pathlib

import pytest

from minitcad.api import REQUIRED_FILES, run, write_required_outputs
from minitcad.io.config import load_config
from minitcad.io.results_csv import read_rows

EXAMPLES = pathlib.Path(__file__).resolve().parents[1] / "examples"
BASELINE = EXAMPLES / "pn_baseline_0_to_030.yaml"
NEWTON_BASELINE = EXAMPLES / "pn_newton_adaptive_0_to_030.yaml"


@pytest.fixture(scope="module")
def config():
    return load_config(BASELINE)


def _stop_after(limit: int):
    """第 `limit + 1` 次检查时返回 True，即“再做 `limit` 次尝试就取消”。"""
    state = {"checks": 0}

    def stop() -> bool:
        state["checks"] += 1
        return state["checks"] > limit

    return stop


def test_cancel_keeps_accepted_points_and_is_not_accepted(config):
    result = run(config, should_stop=_stop_after(3))
    assert result.cancelled
    assert result.decision == "REJECTED"
    assert result.failure_reason == "CANCELLED"
    assert not result.accepted
    assert 2 <= len(result.attempts) <= 3
    # 取消只可能发生在两次尝试之间：已记录的尝试必须都是真解。
    assert all(attempt.result is not None and attempt.result.converged for attempt in result.attempts)
    assert len(result.successful) == len(result.attempts)


def test_cancel_still_writes_the_three_files(config, tmp_path):
    """取消不是“什么都不能交”：已接受的偏压点必须能落成可审计的三文件。"""
    result = run(config, should_stop=_stop_after(2))
    paths = write_required_outputs(result, str(tmp_path))
    assert [pathlib.Path(p).name for p in paths] == list(REQUIRED_FILES)
    for path in paths:
        assert pathlib.Path(path).is_file()

    rows = read_rows(paths[0])
    sweep_rows = [row for row in rows if row["record_type"] == "sweep"]
    assert {row["run_id"] for row in rows} == {result.run_id}
    # 每个尝试偏压写多行（施加电压/迭代数/耗时/收敛标志/终端电流），按偏压去重再比。
    assert len({row["bias_V"] for row in sweep_rows}) == len(result.attempts)
    assert len({row["config_sha256"] for row in rows}) == 1
    for row in rows:
        assert row["status"] == "PASS"


def test_progress_never_reports_unattempted_biases(config):
    """`progress` 的偏价序列必须与 `attempts` 逐个对齐——GUI 的进度条靠这个。"""
    seen: list[tuple[float, str]] = []
    result = run(
        config,
        progress=lambda bias, status: seen.append((float(bias), status)),
        should_stop=_stop_after(3),
    )
    assert [bias for bias, _ in seen] == [attempt.bias_V for attempt in result.attempts]
    assert [status for _, status in seen] == [attempt.status for attempt in result.attempts]
    assert all(bias <= config.sweep.stop_V for bias, _ in seen)


def test_cancel_before_the_first_attempt_has_no_side_effects(config):
    result = run(config, should_stop=lambda: True)
    assert result.attempts == ()
    assert result.successful == ()
    assert result.cancelled and not result.accepted


def test_cancel_in_adaptive_mode():
    """自适应模式多了一条回滚路径，但取消语义必须完全一致。"""
    config = load_config(NEWTON_BASELINE)
    assert config.sweep.adaptive
    result = run(config, should_stop=_stop_after(4))
    assert result.cancelled and result.decision == "REJECTED"
    assert result.failure_reason == "CANCELLED"
    assert result.attempts
    # 自适应路径上失败是**回滚**而不是取消：取消只在两次尝试之间发生，
    # 因此 attempts 里可以有 FAIL，但不允许有“没跑完的”尝试。
    assert all(attempt.result is not None for attempt in result.attempts)
    assert result.attempts[-1].bias_V < config.sweep.stop_V


def test_should_stop_is_not_consulted_inside_the_solver(config):
    """每次尝试只查一次 `should_stop`：求解过程本身没有取消点。"""
    checks = {"count": 0}

    def stop() -> bool:
        checks["count"] += 1
        return False

    result = run(config, should_stop=stop)
    assert len(result.attempts) == len(config.sweep.biases())
    # 固定步进在每次尝试之前查一次，最后一个偏压之后不再查。
    assert checks["count"] == len(result.attempts)


def test_cancelled_run_result_still_reports_a_reason_string(config):
    result = run(config, should_stop=_stop_after(1))
    assert isinstance(result.failure_reason, str) and result.failure_reason == "CANCELLED"
    assert dataclasses.fields(result)  # RunResult 仍是 dataclass，GUI 可以反射字段
