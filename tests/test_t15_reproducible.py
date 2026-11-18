"""T15：可复现性（§14.3）。

同一配置、同一版本与依赖连跑两次，关键值必须一致；运行元数据、配置 SHA 与完整配置
必须能追溯。这里只用固定步进的基准，因为自适应扫压的**尝试序列**在失败边界附近对末位
敏感（`reports/week05_acceptance.md` §4.2），拿它做复现性断言会是不稳定的测试。
"""

from __future__ import annotations

import pathlib

import pytest

from minitcad.api import run, write_required_outputs
from minitcad.io.results_csv import read_rows
from minitcad.io.config import config_fingerprint, load_config

EXAMPLES = pathlib.Path(__file__).resolve().parents[1] / "examples"
BASELINE = EXAMPLES / "pn_baseline_0_to_030.yaml"


@pytest.fixture(scope="module")
def config():
    return load_config(BASELINE)


@pytest.fixture(scope="module")
def first(config):
    return run(config)


def test_config_fingerprint_is_stable(config):
    """同一份 YAML 解析两次得到同一指纹——指纹不能依赖对象地址或时间。"""
    assert config_fingerprint(load_config(BASELINE)) == config_fingerprint(config)
    assert run(config).config_sha256 == config_fingerprint(config)


def test_repeated_runs_agree_on_every_key_value(config, first):
    """两次运行的偏压点、迭代数、残差与终端电流必须逐位相同。"""
    second = run(config)
    assert second.config_sha256 == first.config_sha256
    assert second.decision == first.decision
    assert [attempt.bias_V for attempt in second.attempts] == [
        attempt.bias_V for attempt in first.attempts
    ]
    for left, right in zip(first.attempts, second.attempts):
        assert right.status == left.status
        assert right.result.iterations == left.result.iterations
        assert right.result.residual_norm == left.result.residual_norm
        assert right.result.update_norm == left.result.update_norm
        assert tuple(right.result.state.u) == tuple(left.result.state.u)
        assert tuple(right.result.state.n) == tuple(left.result.state.n)


def test_run_id_is_unique_but_both_runs_are_traceable(config, first):
    """`run_id` 每次不同（时间戳+随机），但两次都必须能和同一 `config_sha256` 对上。"""
    second = run(config)
    assert second.run_id != first.run_id
    assert second.config_sha256 == first.config_sha256


def test_written_files_carry_the_metadata(config, tmp_path):
    """§12.4：`run_id` 与配置哈希必须写进文件的每一行，报告还要有版本信息。"""
    result = run(config)
    csv_path, _, md_path = write_required_outputs(result, str(tmp_path))
    rows = read_rows(csv_path)
    assert {row["run_id"] for row in rows} == {result.run_id}
    assert {row["config_sha256"] for row in rows} == {config_fingerprint(config)}

    report = pathlib.Path(md_path).read_text(encoding="utf-8")
    assert result.run_id in report
    assert result.config_sha256 in report
    for marker in ("- platform:", "- python:", "- decision:"):
        assert marker in report


def test_csv_is_byte_identical_across_runs(config, tmp_path):
    """行内容与顺序也必须稳定：只有 `run_id`/耗时/`runtime_s` 允许不同。"""
    first_run = run(config)
    first_dir = tmp_path / "a"
    write_required_outputs(first_run, str(first_dir))
    second_run = run(config)
    second_dir = tmp_path / "b"
    write_required_outputs(second_run, str(second_dir))

    def normalise(path: pathlib.Path) -> list[tuple[str, ...]]:
        """抹掉两处本来就允许变化的量：`run_id`（时间戳+随机）与墙钟耗时。"""
        out = []
        for row in read_rows(path):
            values = dict(row)
            values["run_id"] = ""
            if values["quantity"] == "runtime_s":
                values["numeric_value"] = ""
            out.append(tuple(values.items()))
        return out

    assert normalise(first_dir / "comparison.csv") == normalise(second_dir / "comparison.csv")
