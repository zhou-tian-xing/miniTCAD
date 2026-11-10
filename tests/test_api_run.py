"""W4-B3/B4：`api.run()` 的公开语义（§10.2）。

覆盖：正常接受、协作式取消、失败原因、`progress` 只报已尝试偏压、
以及 W5 起 `newton` 走同一条流水线。
"""

from __future__ import annotations

import dataclasses
import pathlib

import pytest

from minitcad.api import run, write_required_outputs
from minitcad.core.device import OutputConfig
from minitcad.io.config import config_fingerprint, load_config

_EXAMPLES = pathlib.Path(__file__).resolve().parents[1] / "examples"
BASELINE = _EXAMPLES / "pn_baseline_0_to_030.yaml"
NEWTON_BASELINE = _EXAMPLES / "pn_newton_adaptive_0_to_030.yaml"


@pytest.fixture(scope="module")
def config():
    return load_config(BASELINE)


def test_run_returns_accepted_result_with_traceable_ids(config):
    result = run(config)
    assert result.accepted and result.failure_reason is None
    assert result.run_id and result.config_sha256 == config_fingerprint(config)
    assert result.config is config
    assert result.device.mesh.node_count == 201


def test_progress_reports_only_attempted_biases(config):
    seen: list[tuple[float, str]] = []
    result = run(config, progress=lambda bias, status: seen.append((float(bias), status)))

    assert [bias for bias, _ in seen] == pytest.approx([a.bias_V for a in result.attempts])
    assert {status for _, status in seen} == {"PASS"}
    assert seen[0][0] == 0.0 and seen[-1][0] == pytest.approx(0.3)


def test_should_stop_cancels_before_the_first_attempt(config):
    result = run(config, should_stop=lambda: True)
    assert result.attempts == ()
    assert result.decision == "REJECTED"
    assert result.cancelled and result.failure_reason == "CANCELLED"
    assert not result.accepted


def test_should_stop_keeps_the_accepted_points(config):
    """取消后已接受的偏压点必须完整保留，供 GUI 画已求出的曲线（§10.2）。"""
    calls = {"count": 0}

    def stop() -> bool:
        calls["count"] += 1
        return calls["count"] > 3

    result = run(config, should_stop=stop)
    assert result.cancelled
    assert 1 <= len(result.attempts) <= 4
    assert all(attempt.result.converged for attempt in result.attempts)
    assert len(result.successful) == len(result.attempts)


def test_failure_is_reported_with_bias_and_reason(config):
    strict = dataclasses.replace(config, solver=dataclasses.replace(config.solver, max_iter=1))
    result = run(strict)
    assert not result.accepted
    assert result.failure_reason is not None
    assert "V" in result.failure_reason and ":" in result.failure_reason
    # 失败尝试必须留在 attempts 里，最后一个成功点之前的解仍然可读。
    assert result.attempts
    assert result.attempts[-1].status == "FAIL"


def test_newton_runs_through_the_same_pipeline():
    """W5：`method: newton` 已接入 `run()`；同一器件用自适应 + Newton 走到 0.30 V。

    固定 0.05 V 步长下 Newton 会在 0.10 V 失败（见 `reports/week05_acceptance.md`），
    因此这里用 `examples/pn_newton_adaptive_0_to_030.yaml`。判决只看到没到终点：
    中途回滚的失败尝试留在 `attempts` 里，但整体仍可以是 ACCEPTED。
    """
    newton = load_config(NEWTON_BASELINE)
    assert newton.solver.method == "newton"
    result = run(newton)
    assert result.config.solver.method == "newton"
    assert result.accepted, result.failure_reason
    assert result.attempts[-1].bias_V == pytest.approx(0.3)
    assert result.attempts[-1].status == "PASS"
    # 自适应路径上出现过回滚：失败记录必须留着，不能被“最后成功了”抹掉。
    assert any(attempt.status == "FAIL" for attempt in result.attempts)
    assert len(result.successful) < len(result.attempts)


def test_writer_requires_overwrite_for_existing_files(config, tmp_path):
    result = run(config)
    paths = write_required_outputs(result, str(tmp_path))
    assert all(pathlib.Path(path).is_file() for path in paths)
    with pytest.raises(FileExistsError):
        write_required_outputs(result, str(tmp_path))
    write_required_outputs(result, str(tmp_path), overwrite=True)


def test_si_output_only_changes_the_writer(config, tmp_path):
    """同一器件换成 SI 输出：CSV 表头/量名/单位变，物理量与 run_id 语义不变。"""
    si_config = dataclasses.replace(config, output=OutputConfig(directory=str(tmp_path), unit_system="si"))
    result = run(si_config)
    paths = write_required_outputs(result, str(tmp_path))
    header = pathlib.Path(paths[0]).read_text(encoding="utf-8-sig").splitlines()[0]
    assert "position_m" in header and "position_cm" not in header
    assert result.device.l0_cm == pytest.approx(run(config).device.l0_cm, rel=0, abs=0)
