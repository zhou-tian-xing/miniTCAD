"""W4-B2：标准基准（0 → 0.30 V）的数据快照回归。

快照 `examples/baseline_0_to_030_snapshot.json` 由同一次真实运行生成。
断言只钉住**物理结果**：每个偏压的终端电流密度与总体决定。
迭代次数不进断言——跟踪算法改进会合法地改变它，钉死会阻碍优化。
"""

from __future__ import annotations

import json
import pathlib

import pytest

from minitcad.api import run
from minitcad.io.config import config_fingerprint, load_config
from minitcad.simulation.diagnostics import fields_for_attempt

ROOT = pathlib.Path(__file__).resolve().parents[1]
CONFIG = ROOT / "examples" / "pn_baseline_0_to_030.yaml"
SNAPSHOT = ROOT / "examples" / "baseline_0_to_030_snapshot.json"

CURRENT_RTOL = 1e-6


@pytest.fixture(scope="module")
def snapshot() -> dict:
    return json.loads(SNAPSHOT.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def current_run(snapshot):
    config = load_config(CONFIG)
    assert config_fingerprint(config) == snapshot["config_sha256"], (
        "基准配置被改过：若要更新基准，请重新生成 snapshot 并在报告中记录原因"
    )
    return config, run(config)


def test_decision_matches_snapshot(snapshot, current_run):
    _, result = current_run
    assert result.decision == snapshot["decision"]


def test_every_snapshot_bias_converges_with_the_recorded_current(snapshot, current_run):
    config, result = current_run
    recorded = {item["bias_V"]: item for item in snapshot["biases"]}
    assert [attempt.bias_V for attempt in result.successful] == pytest.approx(list(recorded))

    for attempt in result.successful:
        expected = recorded[attempt.bias_V]["Jterminal_A_cm2"]
        measured = fields_for_attempt(config, result.device, attempt).terminal_current_density_A_cm2()
        if expected == 0.0:
            assert measured == 0.0
        else:
            assert measured == pytest.approx(expected, rel=CURRENT_RTOL), attempt.bias_V


def test_validation_metrics_stay_within_snapshot_order_of_magnitude(snapshot, current_run):
    from minitcad.io.report import validation_metrics

    _, result = current_run
    metrics = validation_metrics(result)
    for key, expected in snapshot["validation"].items():
        if expected == 0.0:
            assert metrics[key] == pytest.approx(0.0, abs=1e-12)
        else:
            assert metrics[key] == pytest.approx(expected, rel=1e-3, abs=1e-12), key
