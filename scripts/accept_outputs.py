"""W7-B3：三个固定输出文件的独立验收器（§12.4）。

只读文件，不 import `io.results_csv`/`io.report`/`io.plot`——验收器必须独立于被验收的
writer，否则 writer 改了字段名、验收器跟着一起改，等于没验。**本文件里硬编码了冻结的
schema**（22 列、各 `record_type` 的 quantity/unit、报告十节标题），这就是契约本身。

用法：

```bash
python scripts/accept_outputs.py runs/pn_baseline_0_to_030
python scripts/accept_outputs.py runs/pn_baseline_0_to_030 --expect REJECTED
python scripts/accept_outputs.py runs/foo --quiet   # 只输出失败项
```

退出码：0 = 全部通过；1 = 有失败项。
"""

from __future__ import annotations

import argparse
import csv
import pathlib
import struct
import sys
from dataclasses import dataclass, field

CSV_NAME = "comparison.csv"
PNG_NAME = "comparison.png"
MD_NAME = "convergence.md"
REQUIRED_FILES = (CSV_NAME, PNG_NAME, MD_NAME)

# §12.1 冻结的 22 列（cm 模式；SI 只把第 9 列改名为 position_m）。
CSV_COLUMNS = (
    "run_id", "config_sha256", "solver", "bias_V", "record_type", "quantity",
    "location_kind", "index", "position_cm", "iteration", "unit", "numeric_value",
    "theory_value", "theory_kind", "theory_abs_error", "theory_rel_error",
    "reference_value", "reference_kind", "reference_abs_error", "reference_rel_error",
    "status", "notes",
)
SI_POSITION_COLUMN = "position_m"

# 每个 record_type 允许出现的 (quantity, unit) 对（§12.1）。
ALLOWED = {
    "cm": {
        "node": {
            ("phi_V", "V"), ("n_cm3", "cm^-3"), ("p_cm3", "cm^-3"),
            ("ND_cm3", "cm^-3"), ("NA_cm3", "cm^-3"), ("net_doping_cm3", "cm^-3"),
            ("charge_C_cm3", "C/cm^3"), ("SRH_cm3_s", "cm^-3/s"),
            ("F_poisson_nd", "1"), ("F_electron_nd", "1"), ("F_hole_nd", "1"),
        },
        "edge": {
            ("E_V_cm", "V/cm"), ("Jn_A_cm2", "A/cm^2"), ("Jp_A_cm2", "A/cm^2"),
            ("Jtotal_A_cm2", "A/cm^2"),
        },
        "sweep": {
            ("applied_voltage_V", "V"), ("iteration_count", "1"), ("runtime_s", "s"),
            ("converged", "1"), ("Jterminal_A_cm2", "A/cm^2"), ("I_A", "A"),
        },
        "iteration": {
            ("residual_norm", "1"), ("update_norm", "1"), ("alpha", "1"),
        },
        "benchmark": {
            ("contact_n_cm3", "cm^-3"), ("contact_p_cm3", "cm^-3"),
            ("contact_potential_V", "V"), ("built_in_potential_V", "V"),
            ("equilibrium_np_ratio_nd", "1"), ("equilibrium_zero_current_A_cm2", "A/cm^2"),
        },
    },
    "si": {
        "node": {
            ("phi_V", "V"), ("n_m3", "m^-3"), ("p_m3", "m^-3"),
            ("ND_m3", "m^-3"), ("NA_m3", "m^-3"), ("net_doping_m3", "m^-3"),
            ("charge_C_m3", "C/m^3"), ("SRH_m3_s", "m^-3/s"),
            ("F_poisson_nd", "1"), ("F_electron_nd", "1"), ("F_hole_nd", "1"),
        },
        "edge": {
            ("E_V_m", "V/m"), ("Jn_A_m2", "A/m^2"), ("Jp_A_m2", "A/m^2"),
            ("Jtotal_A_m2", "A/m^2"),
        },
        "sweep": {
            ("applied_voltage_V", "V"), ("iteration_count", "1"), ("runtime_s", "s"),
            ("converged", "1"), ("Jterminal_A_m2", "A/m^2"), ("I_A", "A"),
        },
        "iteration": {
            ("residual_norm", "1"), ("update_norm", "1"), ("alpha", "1"),
        },
        "benchmark": {
            ("contact_n_m3", "m^-3"), ("contact_p_m3", "m^-3"),
            ("contact_potential_V", "V"), ("built_in_potential_V", "V"),
            ("equilibrium_np_ratio_nd", "1"), ("equilibrium_zero_current_A_m2", "A/m^2"),
        },
    },
}

# 报告必须有的十节（§12.3）。
REPORT_SECTIONS = (
    "## 1. Run summary", "## 2. Orientation", "## 3. Physics / units",
    "## 4. Solver settings", "## 5. Tolerances", "## 6. Bias table",
    "## 7. Validation", "## 8. Error policy", "## 9. Reproducibility",
    "## 10. Final decision",
)

PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"

# 理论值/误差列的浮点比较容差：CSV 是十进制文本，重算必然有末位差。
NUMERIC_TOL = 1e-12


@dataclass
class Finding:
    """一条检查结论；`ok=False` 时才需要看 `detail`。"""

    code: str
    ok: bool
    detail: str = ""


@dataclass
class Acceptance:
    """一次验收的结果。"""

    directory: pathlib.Path
    findings: list[Finding] = field(default_factory=list)

    @property
    def failures(self) -> list[Finding]:
        return [item for item in self.findings if not item.ok]

    @property
    def ok(self) -> bool:
        return not self.failures

    def add(self, code: str, ok: bool, detail: str = "") -> None:
        self.findings.append(Finding(code=code, ok=ok, detail=detail))


def _float(text: str) -> float | None:
    try:
        return float(text)
    except (TypeError, ValueError):
        return None


def _check_png(path: pathlib.Path, result: Acceptance) -> None:
    """PNG 必须可打开且非空：只读文件头，不依赖 Pillow/matplotlib。"""
    data = path.read_bytes()
    if not data.startswith(PNG_SIGNATURE) or len(data) < 33:
        result.add("PNG", False, f"{path.name} is not a valid PNG ({len(data)} bytes)")
        return
    width, height = struct.unpack(">II", data[16:24])
    result.add("PNG", width > 0 and height > 0, f"{width}x{height}, {len(data)} bytes")


def _check_report(path: pathlib.Path, result: Acceptance, expected: str | None) -> str | None:
    text = path.read_text(encoding="utf-8-sig")
    missing = [section for section in REPORT_SECTIONS if section not in text]
    result.add("REPORT_SECTIONS", not missing, f"missing sections: {missing}" if missing else "")

    decision = None
    for line in text.splitlines():
        if line.startswith("- decision:"):
            decision = line.split("**")[1] if "**" in line else line.split(":", 1)[1].strip()
            break
    result.add("REPORT_DECISION", decision in {"ACCEPTED", "REJECTED"}, f"decision={decision!r}")
    if expected is not None:
        result.add(
            "REPORT_EXPECTED",
            decision == expected,
            f"report says {decision!r}, expected {expected!r}",
        )
    return decision


def _check_header(rows: list[dict], columns: list[str], result: Acceptance) -> str:
    cm_header = list(CSV_COLUMNS)
    si_header = [SI_POSITION_COLUMN if name == "position_cm" else name for name in CSV_COLUMNS]
    if columns == cm_header:
        system = "cm"
    elif columns == si_header:
        system = "si"
    else:
        result.add("HEADER", False, f"unexpected columns: {columns}")
        return "cm"
    result.add("HEADER", True, f"{len(columns)} columns, {system}")
    if rows:
        result.add("HEADER_ROWS", all(set(row) == set(columns) for row in rows))
    return system


def _check_metadata(rows: list[dict], result: Acceptance) -> None:
    run_ids = {row["run_id"] for row in rows}
    shas = {row["config_sha256"] for row in rows}
    solvers = {row["solver"] for row in rows}
    result.add("RUN_ID", len(run_ids) == 1, f"{sorted(run_ids)}")
    result.add("CONFIG_SHA", len(shas) == 1 and len(next(iter(shas), "")) == 64, f"{sorted(shas)}")
    result.add("SOLVER", len(solvers) == 1, f"{sorted(solvers)}")


def _check_units(rows: list[dict], system: str, result: Acceptance) -> None:
    bad = []
    for row in rows:
        allowed = ALLOWED[system].get(row["record_type"])
        if allowed is None:
            bad.append(f"{row['record_type']}: unknown record_type")
        elif (row["quantity"], row["unit"]) not in allowed:
            bad.append(f"{row['record_type']}: {row['quantity']} [{row['unit']}]")
    unique = sorted(set(bad))
    result.add("QUANTITY_UNIT", not unique, f"{unique[:6]} (and {len(unique)} kinds)" if unique else "")


def _check_theory(rows: list[dict], result: Acceptance) -> None:
    """`theory_kind` 与理论/误差列必须自洽（§12.1、§12.4）。"""
    problems = []
    for row in rows:
        kind = row["theory_kind"]
        numeric = _float(row["numeric_value"])
        theory = _float(row["theory_value"])
        abs_error = _float(row["theory_abs_error"])
        rel_error = _float(row["theory_rel_error"])
        if kind == "none":
            if any(value is not None for value in (theory, abs_error, rel_error)):
                problems.append(f"{row['quantity']}: theory_kind=none but values present")
            continue
        if kind not in {"exact", "reference"}:
            problems.append(f"{row['quantity']}: unknown theory_kind {kind!r}")
            continue
        if numeric is None or theory is None:
            problems.append(f"{row['quantity']}: theory_kind={kind} needs numeric and theory")
            continue
        if abs_error is None or abs(abs_error - abs(numeric - theory)) > NUMERIC_TOL * max(
            1.0, abs(numeric), abs(theory)
        ):
            problems.append(f"{row['quantity']}: abs_error inconsistent ({abs_error} vs {numeric - theory})")
        if theory != 0.0:
            expected = abs(numeric - theory) / abs(theory)
            if rel_error is None or abs(rel_error - expected) > NUMERIC_TOL * max(1.0, expected):
                problems.append(f"{row['quantity']}: rel_error inconsistent ({rel_error} vs {expected})")
        elif rel_error is not None:
            problems.append(f"{row['quantity']}: zero theory must leave rel_error empty")
    unique = sorted(set(problems))
    result.add("THEORY_CONSISTENCY", not unique, f"{unique[:6]} (and {len(unique)} kinds)" if unique else "")


def _check_counts(rows: list[dict], system: str, result: Acceptance) -> None:
    """每个收敛偏压：node 行 `M` 个、edge 行 `M-1` 个，每种 quantity 各一套。"""
    position = "position_cm" if system == "cm" else SI_POSITION_COLUMN
    node_quantities = {quantity for quantity, _ in ALLOWED[system]["node"]}
    edge_quantities = {quantity for quantity, _ in ALLOWED[system]["edge"]}

    passed_biases = sorted(
        {
            float(row["bias_V"])
            for row in rows
            if row["record_type"] == "sweep" and row["status"] == "PASS"
        }
    )
    if not passed_biases:
        result.add("DATA_COMPLETE", False, "no PASS bias point in the sweep rows")
        return

    sizes = {}
    problems = []
    for bias in passed_biases:
        selected = [row for row in rows if row["record_type"] in {"node", "edge"} and float(row["bias_V"]) == bias]
        node_rows = [row for row in selected if row["record_type"] == "node"]
        edge_rows = [row for row in selected if row["record_type"] == "edge"]
        if not node_rows or not edge_rows:
            problems.append(f"{bias:g} V: missing node/edge rows")
            continue
        node_count = len({row[position] for row in node_rows if row["quantity"] == "phi_V"})
        edge_count = len({row[position] for row in edge_rows if row["quantity"] == "E_V_cm" or row["quantity"] == "E_V_m"})
        sizes[bias] = (node_count, edge_count)
        if edge_count != node_count - 1:
            problems.append(f"{bias:g} V: {node_count} nodes but {edge_count} edges")
        for quantity in node_quantities:
            count = len({row[position] for row in node_rows if row["quantity"] == quantity})
            if count != node_count:
                problems.append(f"{bias:g} V: {quantity} has {count} rows, expected {node_count}")
        for quantity in edge_quantities:
            count = len({row[position] for row in edge_rows if row["quantity"] == quantity})
            if count != edge_count:
                problems.append(f"{bias:g} V: {quantity} has {count} rows, expected {edge_count}")

    uniform_sizes = {size for size in sizes.values()}
    if len(uniform_sizes) != 1:
        problems.append(f"node counts differ between bias points: {sorted(uniform_sizes)}")
    result.add("DATA_COMPLETE", not problems, f"{sorted(problems)[:4]}" if problems else "")


def _check_iterations(rows: list[dict], result: Acceptance) -> None:
    """迭代行按 `(bias, iteration)` 分组，每个迭代号必须同时有 residual/update/alpha。"""
    groups: dict[tuple[str, str], set[str]] = {}
    for row in rows:
        if row["record_type"] != "iteration":
            continue
        key = (row["bias_V"], row["iteration"])
        groups.setdefault(key, set()).add(row["quantity"])
    incomplete = [key for key, quantities in groups.items() if len(quantities) != 3]
    result.add(
        "ITERATION_ROWS",
        bool(groups) and not incomplete,
        f"{len(groups)} iteration indices" + (f", incomplete: {incomplete[:3]}" if incomplete else ""),
    )


def _check_decision(rows: list[dict], report_decision: str | None, result: Acceptance) -> None:
    """最终判决必须等于**最后一次尝试**的结果（§12.4）。

    自适应扫压允许中途 `FAIL`（回滚后继续），所以“有 FAIL 就是 REJECTED”是错的；
    判据是最后那次尝试成功与否：成功到终点 → ACCEPTED，停在失败点 → REJECTED。
    `applied_voltage_V` 行按尝试顺序写出，最后一条就是最后一次尝试。
    """
    if report_decision is None:
        return
    attempts = [
        row
        for row in rows
        if row["record_type"] == "sweep" and row["quantity"] == "applied_voltage_V"
    ]
    if not attempts:
        result.add("DECISION_CONSISTENT", False, "no applied_voltage_V sweep rows")
        return
    last = attempts[-1]
    expected = "ACCEPTED" if last["status"] == "PASS" else "REJECTED"
    result.add(
        "DECISION_CONSISTENT",
        report_decision == expected,
        f"last attempt {last['bias_V']} V is {last['status']}, report says {report_decision}",
    )


def _check_failure_reasons(rows: list[dict], result: Acceptance) -> None:
    """失败尝试必须留下原因：§12.4 要求 FAIL 记录可追溯，不能只是缺数据。"""
    silent = [
        row["bias_V"]
        for row in rows
        if row["record_type"] == "sweep"
        and row["quantity"] == "applied_voltage_V"
        and row["status"] == "FAIL"
        and not row["notes"]
    ]
    result.add(
        "FAILURE_NOTES",
        not silent,
        f"FAIL at {sorted(set(silent))} V has no notes" if silent else "",
    )


def accept(directory: str | pathlib.Path, *, expected_decision: str | None = None) -> Acceptance:
    """对目录里的三个文件跑一遍 §12.4 的检查。"""
    root = pathlib.Path(directory)
    result = Acceptance(directory=root)

    missing = [name for name in REQUIRED_FILES if not (root / name).is_file()]
    result.add("FILES", not missing, f"missing: {missing}" if missing else "")
    if missing:
        return result
    empty = [name for name in REQUIRED_FILES if (root / name).stat().st_size == 0]
    result.add("FILE_NOT_EMPTY", not empty, f"empty: {empty}" if empty else "")

    _check_png(root / PNG_NAME, result)
    report_decision = _check_report(root / MD_NAME, result, expected_decision)

    with open(root / CSV_NAME, "r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        columns = list(reader.fieldnames or ())
        rows = list(reader)
    result.add("CSV_ROWS", bool(rows), f"{len(rows)} rows")
    if not rows:
        return result

    system = _check_header(rows, columns, result)
    _check_metadata(rows, result)
    _check_units(rows, system, result)
    _check_theory(rows, result)
    _check_counts(rows, system, result)
    _check_iterations(rows, result)
    _check_decision(rows, report_decision, result)
    _check_failure_reasons(rows, result)
    return result


def _print(result: Acceptance, quiet: bool) -> None:
    if not quiet:
        print(f"accepting {result.directory}")
    for finding in result.findings:
        if quiet and finding.ok:
            continue
        mark = "PASS" if finding.ok else "FAIL"
        detail = f"  {finding.detail}" if finding.detail else ""
        print(f"  [{mark}] {finding.code}{detail}")
    total = len(result.findings)
    failed = len(result.failures)
    print(f"{'ACCEPTED' if result.ok else 'REJECTED'}: {total - failed}/{total} checks passed")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Accept the three fixed output files (section 12.4).")
    parser.add_argument("directory", help="directory holding comparison.csv/.png/convergence.md")
    parser.add_argument("--expect", choices=("ACCEPTED", "REJECTED"), default=None)
    parser.add_argument("--quiet", action="store_true", help="only report failures")
    args = parser.parse_args(argv)

    result = accept(args.directory, expected_decision=args.expect)
    _print(result, args.quiet)
    return 0 if result.ok else 1


if __name__ == "__main__":
    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
    raise SystemExit(main())
