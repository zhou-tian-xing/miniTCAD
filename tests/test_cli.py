"""W3-B3：CLI 骨架。validate 真跑，solve 明确未接线（第 13.1 节）。"""

from __future__ import annotations

import pathlib
import textwrap

import pytest

from minitcad import cli
from minitcad.io.results_csv import CSV_COLUMNS, write_rows

EXAMPLES = pathlib.Path(__file__).resolve().parents[1] / "examples"
CM_YAML = EXAMPLES / "pn_left_p_right_n.yaml"
SI_YAML = EXAMPLES / "pn_left_p_right_n_si.yaml"


def test_validate_yaml_prints_banner_and_summary(capsys):
    assert cli.main(["validate", str(CM_YAML)]) == 0
    out = capsys.readouterr().out
    assert "PN orientation: P(left, x=0) -> N(right, x=L)" in out
    assert "core units   : dimensionless" in out
    assert "input units  : cm" in out
    assert "nodes        : 201" in out
    assert "C_L / C_R" in out


def test_validate_si_yaml_reports_si_input(capsys):
    assert cli.main(["validate", str(SI_YAML)]) == 0
    assert "input units  : si" in capsys.readouterr().out


def test_validate_rejects_broken_config_without_traceback(tmp_path, capsys):
    broken = tmp_path / "broken.yaml"
    broken.write_text(
        textwrap.dedent(
            """
            device:
              type: pn
              length_cm: 1.0e-3
              area_cm2: 1.0e-4
              temperature_K: 300.0
            material:
              name: silicon
            regions:
              - {start_cm: 0.0, end_cm: 1.0e-3, donor_cm3: 1.0e16}
            mesh: {type: uniform, nodes: 11}
            normalization: {reference_diffusivity_cm2_s: 10.0}
            """
        ),
        encoding="utf-8",
    )
    assert cli.main(["validate", str(broken)]) == 2
    captured = capsys.readouterr()
    assert captured.out == ""
    assert "REJECTED:" in captured.err
    assert "units.system" in captured.err


BASELINE_YAML = EXAMPLES / "pn_baseline_0_to_030.yaml"
NEWTON_YAML = EXAMPLES / "pn_newton_adaptive_0_to_030.yaml"


def test_solve_writes_the_three_required_files(tmp_path, capsys):
    """W4-B4：solve 走 run() + writer，产出第 12 节的三文件（§12.4）。"""
    assert cli.main(["solve", str(BASELINE_YAML), "--out", str(tmp_path)]) == 0
    out = capsys.readouterr().out
    assert "decision: ACCEPTED" in out
    for name in ("comparison.csv", "comparison.png", "convergence.md"):
        assert (tmp_path / name).is_file(), name
    assert cli.main(["validate", str(tmp_path / "comparison.csv")]) == 0


def test_solve_refuses_to_overwrite_existing_outputs(tmp_path, capsys):
    assert cli.main(["solve", str(BASELINE_YAML), "--out", str(tmp_path)]) == 0
    capsys.readouterr()
    assert cli.main(["solve", str(BASELINE_YAML), "--out", str(tmp_path)]) == 2
    assert "already contains" in capsys.readouterr().err
    assert cli.main(["solve", str(BASELINE_YAML), "--out", str(tmp_path), "--overwrite"]) == 0


def test_solve_accepts_newton_with_adaptive_sweep(tmp_path, capsys):
    """W5：`method: newton` 已接线，CLI 走 `run()` 的同一条路径并写出三文件。"""
    config = tmp_path / "newton.yaml"
    config.write_text(NEWTON_YAML.read_text(encoding="utf-8"), encoding="utf-8")
    out = tmp_path / "out"
    assert cli.main(["solve", str(config), "--out", str(out)]) == 0
    captured = capsys.readouterr().out
    assert "decision: ACCEPTED" in captured
    assert "solver  : newton" in captured
    for name in ("comparison.csv", "comparison.png", "convergence.md"):
        assert (out / name).is_file()


def test_solve_failure_is_reported_with_bias(tmp_path, capsys):
    """强偏压 + 单次迭代：必须 REJECTED，并在 stderr/摘要里给出偏压与原因。"""
    text = BASELINE_YAML.read_text(encoding="utf-8")
    text = text.replace("max_iter: 60", "max_iter: 1").replace("stop_V: 0.3", "stop_V: 0.3")
    config = tmp_path / "fail.yaml"
    config.write_text(text, encoding="utf-8")
    assert cli.main(["solve", str(config), "--out", str(tmp_path / "out")]) == 2
    out = capsys.readouterr().out
    assert "decision: REJECTED" in out
    assert "reason" in out


def test_validate_missing_file_is_rejected(tmp_path, capsys):
    assert cli.main(["validate", str(tmp_path / "nope.yaml")]) == 2
    assert "no such file" in capsys.readouterr().err


def test_validate_csv_accepts_frozen_header(tmp_path, capsys):
    path = tmp_path / "comparison.csv"
    write_rows(path, [])
    assert cli.main(["validate", str(path)]) == 0
    out = capsys.readouterr().out
    assert f"header OK: {len(CSV_COLUMNS)} columns" in out


def test_validate_csv_rejects_broken_header(tmp_path, capsys):
    path = tmp_path / "comparison.csv"
    path.write_text("run_id,config_sha256,bias_V\n", encoding="utf-8-sig")
    assert cli.main(["validate", str(path)]) == 2
    assert "missing column" in capsys.readouterr().err


def test_argparse_requires_a_subcommand():
    with pytest.raises(SystemExit):
        cli.main([])
