"""W4-B4 / §12.1–§12.4：三个固定文件的 schema 与内容契约。

只检查 writer 的可验收性质：22 列与顺序、失败偏压也留记录、缺失理论留空而不写 0、
报告十节齐全、PNG 是真实图片、SI 输出只改 writer 层。
"""

from __future__ import annotations

import csv
import dataclasses
import pathlib

import pytest

from minitcad.api import run, write_required_outputs
from minitcad.core.device import OutputConfig
from minitcad.io.config import config_fingerprint, load_config
from minitcad.io.results_csv import CSV_COLUMNS, columns_for, read_header

BASELINE = pathlib.Path(__file__).resolve().parents[1] / "examples" / "pn_baseline_0_to_030.yaml"


@pytest.fixture(scope="module")
def written(tmp_path_factory):
    directory = tmp_path_factory.mktemp("w4_outputs")
    result = run(load_config(BASELINE))
    paths = write_required_outputs(result, str(directory))
    return result, [pathlib.Path(path) for path in paths]


def _rows(path: pathlib.Path) -> list[dict[str, str]]:
    with open(path, "r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def test_three_files_exist_and_are_non_empty(written):
    _, paths = written
    for path in paths:
        assert path.is_file() and path.stat().st_size > 0, path


def test_csv_header_is_the_frozen_schema(written):
    _, paths = written
    assert read_header(paths[0]) == CSV_COLUMNS


def test_csv_records_every_attempt_and_expected_row_counts(written):
    result, paths = written
    rows = _rows(paths[0])
    biases = {float(row["bias_V"]) for row in rows if row["record_type"] == "sweep"}
    assert biases == {attempt.bias_V for attempt in result.attempts}

    nodes = [row for row in rows if row["record_type"] == "node"]
    edges = [row for row in rows if row["record_type"] == "edge"]
    m = result.device.mesh.node_count
    accepted = len(result.successful)
    assert len(nodes) == 11 * m * accepted
    assert len(edges) == 4 * (m - 1) * accepted


def test_csv_never_writes_fake_numbers(written):
    _, paths = written
    for row in _rows(paths[0]):
        for column in ("numeric_value", "theory_value", "reference_value"):
            cell = row[column]
            assert cell.lower() not in {"nan", "none", "null"}, (row["quantity"], column, cell)


def test_csv_leaves_theory_empty_when_there_is_no_closed_form(written):
    _, paths = written
    node_rows = [row for row in _rows(paths[0]) if row["record_type"] == "node"]
    assert node_rows
    for row in node_rows:
        assert row["theory_kind"] == "none"
        assert row["theory_value"] == "" and row["theory_rel_error"] == ""


def test_csv_zero_theory_uses_absolute_error_only(written):
    _, paths = written
    zero_current = [
        row
        for row in _rows(paths[0])
        if row["record_type"] == "sweep" and row["quantity"] == "Jterminal_A_cm2" and row["bias_V"] == "0.0"
    ]
    assert len(zero_current) == 1
    row = zero_current[0]
    assert row["theory_value"] == "0.0"
    assert row["theory_rel_error"] == ""
    assert "zero theory" in row["notes"]


def test_csv_benchmark_rows_are_exact_with_zero_error(written):
    _, paths = written
    benchmarks = [row for row in _rows(paths[0]) if row["record_type"] == "benchmark"]
    quantities = {row["quantity"] for row in benchmarks}
    assert {"contact_n_cm3", "contact_p_cm3", "contact_potential_V", "built_in_potential_V"} <= quantities
    for row in benchmarks:
        assert row["theory_kind"] == "exact"
        assert float(row["theory_rel_error"] or 0.0) == pytest.approx(0.0)


def test_positions_are_inside_the_device(written):
    result, paths = written
    positions = [
        float(row["position_cm"])
        for row in _rows(paths[0])
        if row["record_type"] == "node" and row["quantity"] == "phi_V" and row["bias_V"] == "0.0"
    ]
    assert len(positions) == result.device.mesh.node_count
    assert positions[0] == pytest.approx(0.0)
    assert positions[-1] == pytest.approx(result.config.length_cm, rel=1e-12)
    assert all(later > earlier for earlier, later in zip(positions, positions[1:]))


def test_convergence_report_has_the_ten_sections(written):
    result, paths = written
    text = paths[2].read_text(encoding="utf-8")
    for index in range(1, 11):
        assert f"## {index}." in text, index
    assert result.run_id in text
    assert config_fingerprint(result.config) in text
    assert f"**{result.decision}**" in text


def test_png_is_a_real_image(written):
    _, paths = written
    assert paths[1].read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"


def test_si_output_renames_only_the_writer_surface(written, tmp_path):
    result, cm_paths = written
    si_config = dataclasses.replace(result.config, output=OutputConfig(directory=str(tmp_path), unit_system="si"))
    si_result = run(si_config)
    si_paths = write_required_outputs(si_result, str(tmp_path))
    si_csv = pathlib.Path(si_paths[0])

    assert read_header(si_csv) == columns_for("si")
    cm_rows = _rows(cm_paths[0])
    si_rows = _rows(si_csv)

    def pick(rows, quantity, bias, kind="node", index="0"):
        return next(
            row for row in rows
            if row["record_type"] == kind and row["quantity"] == quantity
            and row["bias_V"] == bias and row["index"] == index
        )

    cm_node = pick(cm_rows, "n_cm3", "0.0")
    si_node = pick(si_rows, "n_m3", "0.0")
    assert float(si_node["numeric_value"]) == pytest.approx(float(cm_node["numeric_value"]) * 1e6, rel=1e-12)
    assert si_node["unit"] == "m^-3"
    assert float(si_node["position_m"]) == pytest.approx(float(cm_node["position_cm"]) * 1e-2, rel=1e-12)

    # 同一个无量纲解：SI 只是输出层换算，物理量必须一致。
    assert si_result.device.l0_cm == result.device.l0_cm
