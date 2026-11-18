"""T14：三个固定输出文件的独立验收（§12.4）。

`scripts/accept_outputs.py` 不 import writer，所以这里既验“正常产物能通过”，也逐条伪造
损坏产物，确认每一种违规都能被抓到——只测前者等于只测了写入，没测验收。
"""

from __future__ import annotations

import csv
import pathlib
import shutil
import sys

import pytest

_ROOT = pathlib.Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from minitcad.api import run, write_required_outputs  # noqa: E402
from minitcad.io.config import load_config  # noqa: E402
from scripts.accept_outputs import accept  # noqa: E402

EXAMPLES = _ROOT / "examples"
BASELINE = EXAMPLES / "pn_baseline_0_to_030.yaml"
SI = EXAMPLES / "pn_left_p_right_n_si.yaml"
NEWTON = EXAMPLES / "pn_newton_adaptive_0_to_030.yaml"


def _write(config_path: pathlib.Path, directory: pathlib.Path) -> pathlib.Path:
    result = run(load_config(config_path))
    write_required_outputs(result, str(directory))
    return directory


@pytest.fixture(scope="module")
def accepted(tmp_path_factory):
    return _write(BASELINE, tmp_path_factory.mktemp("accepted"))


@pytest.fixture(scope="module")
def si_output(tmp_path_factory):
    return _write(SI, tmp_path_factory.mktemp("si"))


@pytest.fixture(scope="module")
def rejected(tmp_path_factory):
    import dataclasses

    config = load_config(BASELINE)
    config = dataclasses.replace(config, solver=dataclasses.replace(config.solver, method="newton"))
    directory = tmp_path_factory.mktemp("rejected")
    write_required_outputs(run(config), str(directory))
    return directory


def _codes(result) -> set[str]:
    return {finding.code for finding in result.failures}


def test_accepted_output_passes_every_check(accepted):
    result = accept(accepted)
    assert result.ok, [finding.code for finding in result.failures]
    assert len(result.findings) >= 15


def test_si_output_passes_with_the_si_schema(si_output):
    """§12.4：默认 cm 与显式 SI 各自按确定的 schema 验收。"""
    result = accept(si_output)
    assert result.ok, [finding.code for finding in result.failures]
    header = (si_output / "comparison.csv").read_text(encoding="utf-8-sig").splitlines()[0]
    assert "position_m" in header and "position_cm" not in header


def test_rejected_run_is_recognised_not_rejected_for_having_files(rejected):
    """有 FAIL 记录、最终判决 REJECTED 的运行，文件本身仍应通过验收。"""
    assert accept(rejected, expected_decision="REJECTED").ok
    codes = _codes(accept(rejected, expected_decision="ACCEPTED"))
    assert "REPORT_EXPECTED" in codes


def test_missing_file_is_caught(accepted, tmp_path):
    target = tmp_path / "out"
    shutil.copytree(accepted, target)
    (target / "comparison.png").unlink()
    assert "FILES" in _codes(accept(target))


def test_corrupt_png_is_caught(accepted, tmp_path):
    target = tmp_path / "out"
    shutil.copytree(accepted, target)
    (target / "comparison.png").write_bytes(b"not a png at all")
    assert "PNG" in _codes(accept(target))


def test_missing_report_section_is_caught(accepted, tmp_path):
    target = tmp_path / "out"
    shutil.copytree(accepted, target)
    path = target / "convergence.md"
    text = path.read_text(encoding="utf-8").replace("## 7. Validation", "## 7. Checks")
    path.write_text(text, encoding="utf-8")
    assert "REPORT_SECTIONS" in _codes(accept(target))


def _rewrite(directory: pathlib.Path, mutate) -> pathlib.Path:
    """读出 CSV、就地改若干行、按原表头写回。"""
    path = directory / "comparison.csv"
    with open(path, "r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        columns = list(reader.fieldnames or ())
        rows = list(reader)
    rows = mutate(rows)
    with open(path, "w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return directory


def test_inconsistent_theory_error_is_caught(accepted, tmp_path):
    """把数值改掉而误差列不动，验收器必须重算出不一致。"""
    target = tmp_path / "out"
    shutil.copytree(accepted, target)

    def mutate(rows):
        for row in rows:
            if row["quantity"] == "built_in_potential_V":
                row["numeric_value"] = str(float(row["numeric_value"]) + 0.1)
                break
        return rows

    _rewrite(target, mutate)
    assert "THEORY_CONSISTENCY" in _codes(accept(target))


def test_theory_values_on_a_none_row_are_caught(accepted, tmp_path):
    target = tmp_path / "out"
    shutil.copytree(accepted, target)

    def mutate(rows):
        for row in rows:
            if row["quantity"] == "runtime_s":
                row["theory_kind"] = "none"
                row["theory_value"] = "1.0"
                break
        return rows

    _rewrite(target, mutate)
    assert "THEORY_CONSISTENCY" in _codes(accept(target))


def test_dropped_node_rows_are_caught_as_data_incomplete(accepted, tmp_path):
    target = tmp_path / "out"
    shutil.copytree(accepted, target)

    def mutate(rows):
        dropped = 0
        kept = []
        for row in rows:
            if row["record_type"] == "node" and row["quantity"] == "phi_V" and dropped < 5:
                dropped += 1
                continue
            kept.append(row)
        return kept

    _rewrite(target, mutate)
    assert "DATA_COMPLETE" in _codes(accept(target))


def test_wrong_quantity_unit_is_caught(accepted, tmp_path):
    target = tmp_path / "out"
    shutil.copytree(accepted, target)

    def mutate(rows):
        for row in rows:
            if row["quantity"] == "phi_V":
                row["unit"] = "V/cm"
                break
        return rows

    _rewrite(target, mutate)
    assert "QUANTITY_UNIT" in _codes(accept(target))


def test_decision_that_contradicts_the_last_attempt_is_caught(accepted, tmp_path):
    target = tmp_path / "out"
    shutil.copytree(accepted, target)

    def mutate(rows):
        for row in rows:
            if row["record_type"] == "sweep" and row["quantity"] == "applied_voltage_V":
                row["status"] = "FAIL"
        return rows

    _rewrite(target, mutate)
    codes = _codes(accept(target))
    assert "DECISION_CONSISTENT" in codes
    assert "FAILURE_NOTES" in codes


def test_rollback_failures_do_not_force_rejected(tmp_path):
    """自适应回滚允许中途 FAIL：判决只看最后一次尝试（W5 的语义不能被 T14 判错）。"""
    directory = _write(NEWTON, tmp_path / "newton")
    result = accept(directory, expected_decision="ACCEPTED")
    assert result.ok, [finding.code for finding in result.failures]
    statuses = {
        row["status"]
        for row in csv.DictReader(
            open(directory / "comparison.csv", "r", encoding="utf-8-sig", newline="")
        )
        if row["record_type"] == "sweep" and row["quantity"] == "applied_voltage_V"
    }
    assert "FAIL" in statuses, "这个基准在中途确实回滚过"
