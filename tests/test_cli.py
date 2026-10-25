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


def test_solve_is_honestly_not_wired(capsys):
    assert cli.main(["solve", str(CM_YAML)]) == 3
    captured = capsys.readouterr()
    assert captured.out == ""
    assert "not wired yet" in captured.err


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
