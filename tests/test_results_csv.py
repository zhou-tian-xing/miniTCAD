"""W3-B4：comparison.csv 的 22 列 long-form schema（第 12.1 节）。"""

from __future__ import annotations

import csv

import pytest

from minitcad.io.results_csv import (
    CSV_COLUMNS,
    check_header,
    format_row,
    header_line,
    read_header,
    write_rows,
)

FROZEN = (
    "run_id,config_sha256,solver,bias_V,record_type,quantity,location_kind,index,"
    "position_cm,iteration,unit,numeric_value,theory_value,theory_kind,"
    "theory_abs_error,theory_rel_error,reference_value,reference_kind,"
    "reference_abs_error,reference_rel_error,status,notes"
)


def test_header_matches_the_frozen_22_columns():
    assert header_line() == FROZEN
    assert len(CSV_COLUMNS) == 22


def test_missing_values_become_empty_fields_not_zero():
    row = format_row(
        run_id="r1",
        config_sha256="abc",
        solver="newton",
        bias_V=0.2,
        record_type="node",
        quantity="phi_V",
        location_kind="node",
        index=42,
        position_cm=2.5e-4,
        iteration=None,
        unit="V",
        numeric_value=0.5123,
        theory_value=None,
        theory_kind="none",
        theory_abs_error=None,
        theory_rel_error=None,
        reference_value=None,
        reference_kind="none",
        reference_abs_error=None,
        reference_rel_error=None,
        status="PASS",
        notes="bias point 0.2 V",
    )
    assert len(row) == 22
    assert row[CSV_COLUMNS.index("theory_value")] == ""
    assert row[CSV_COLUMNS.index("iteration")] == ""
    assert row[CSV_COLUMNS.index("numeric_value")] == "0.5123"


def test_booleans_are_written_as_one_zero():
    values = {column: None for column in CSV_COLUMNS}
    values["record_type"] = "sweep"
    values["quantity"] = "converged"
    values["numeric_value"] = True
    assert format_row(**values)[CSV_COLUMNS.index("numeric_value")] == "1"


def test_missing_or_extra_columns_are_rejected():
    with pytest.raises(ValueError):
        format_row(run_id="r1")
    values = {column: None for column in CSV_COLUMNS}
    values["density"] = 1e16
    with pytest.raises(ValueError):
        format_row(**values)


def test_write_then_read_header_round_trip(tmp_path):
    path = tmp_path / "comparison.csv"
    write_rows(path, [])
    assert read_header(path) == CSV_COLUMNS
    # utf-8-sig：文件头带 BOM，但解析后列名不含 BOM
    assert path.read_bytes().startswith(b"\xef\xbb\xbf")


def test_write_rows_quotes_fields_containing_commas(tmp_path):
    path = tmp_path / "comparison.csv"
    values = {column: None for column in CSV_COLUMNS}
    values.update(run_id="r1", status="PASS", notes="a,b")
    write_rows(path, [values])
    with open(path, encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.reader(handle))
    assert rows[1][CSV_COLUMNS.index("notes")] == "a,b"


def test_check_header_reports_order_and_missing():
    assert check_header(CSV_COLUMNS) == []
    problems = check_header(CSV_COLUMNS[:-1])
    assert problems and "missing column" in problems[0]
    swapped = list(CSV_COLUMNS)
    swapped[0], swapped[1] = swapped[1], swapped[0]
    assert check_header(swapped) == ["column order differs from the frozen long-form schema"]
