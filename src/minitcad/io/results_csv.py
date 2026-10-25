"""comparison.csv 的 long-form schema（第 12.1 节，22 列，不可改名）。

W3 只冻结表头与序列化器；真正的数据行由 W4 之后的 writer 填充。
缺失值一律写空字段，不写 0、不写字符串 NaN。
"""

from __future__ import annotations

import csv
import pathlib
from typing import Iterable, Mapping, Sequence

# 默认 cm 模式强制列，顺序即表头顺序（第 12.1 节）。
CSV_COLUMNS: tuple[str, ...] = (
    "run_id",
    "config_sha256",
    "solver",
    "bias_V",
    "record_type",
    "quantity",
    "location_kind",
    "index",
    "position_cm",
    "iteration",
    "unit",
    "numeric_value",
    "theory_value",
    "theory_kind",
    "theory_abs_error",
    "theory_rel_error",
    "reference_value",
    "reference_kind",
    "reference_abs_error",
    "reference_rel_error",
    "status",
    "notes",
)

RECORD_TYPES = ("node", "edge", "sweep", "benchmark", "iteration", "diagnostic")
THEORY_KINDS = ("exact", "approximation", "none")
REFERENCE_KINDS = ("refined_mesh", "external_solver", "none")
STATUSES = ("PASS", "FAIL", "N_A")


def header_line() -> str:
    """返回表头行（字符串形式，便于报告与测试直接比对）。"""
    return ",".join(CSV_COLUMNS)


def format_row(**values: object) -> list[str]:
    """按冻结顺序把关键字参数排成一行；None 写成空字段。

    缺列或多列都报错，避免悄悄写出偏离 schema 的 CSV。
    """
    missing = [column for column in CSV_COLUMNS if column not in values]
    extra = [column for column in values if column not in CSV_COLUMNS]
    if missing or extra:
        raise ValueError(f"missing={missing} extra={sorted(extra)}")
    row = []
    for column in CSV_COLUMNS:
        value = values[column]
        if value is None:
            row.append("")
        elif isinstance(value, bool):
            row.append("1" if value else "0")
        else:
            row.append(str(value))
    return row


def write_rows(path: str | pathlib.Path, rows: Iterable[Sequence[str] | Mapping[str, object]]) -> None:
    """写 UTF-8 with BOM（utf-8-sig）的 CSV，便于中文 Excel 直接打开。"""
    with open(path, "w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(CSV_COLUMNS)
        for row in rows:
            writer.writerow(format_row(**row) if isinstance(row, Mapping) else list(row))


def read_header(path: str | pathlib.Path) -> tuple[str, ...]:
    """读回表头；BOM 由 utf-8-sig 自动去掉。"""
    with open(path, "r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle)
        return tuple(next(reader, ()))


def check_header(header: Sequence[str]) -> list[str]:
    """对照冻结 schema 检查表头，返回问题列表；空列表表示通过。"""
    columns = list(header)
    if columns == list(CSV_COLUMNS):
        return []
    problems = []
    missing = [column for column in CSV_COLUMNS if column not in columns]
    extra = [column for column in columns if column not in CSV_COLUMNS]
    if missing:
        problems.append(f"missing column(s): {missing}")
    if extra:
        problems.append(f"unknown column(s): {extra}")
    if not missing and not extra:
        problems.append("column order differs from the frozen long-form schema")
    return problems
