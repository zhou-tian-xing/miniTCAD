"""`comparison.csv` 的 long-form schema 与数据行构造（第 12.1 节，22 列）。

表头在 W3 冻结；W4 起本模块把 `RunResult` 摊平成行：
每个尝试偏压都有 `sweep` 行（失败也写），已收敛点另有 `node`/`edge`/`iteration` 行，
器件级解析量写 `benchmark` 行。缺失值一律空字段，不写 `0`、不写字符串 `NaN`。

writer 只转换**已恢复的物理量**：cm 模式直接用 cm 工程单位；`output.unit_system: si`
时只改 `position_m` 列名、量名后缀、`unit` 与数值因子，列序不动（§12.1）。
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

# 每个量在 cm 模式下的 unit；SI 模式只换后缀与单位表，不换列序。
QUANTITY_UNITS_CM: dict[str, str] = {
    "phi_V": "V",
    "n_cm3": "cm^-3",
    "p_cm3": "cm^-3",
    "ND_cm3": "cm^-3",
    "NA_cm3": "cm^-3",
    "net_doping_cm3": "cm^-3",
    "charge_C_cm3": "C/cm^3",
    "SRH_cm3_s": "cm^-3/s",
    "F_poisson_nd": "1",
    "F_electron_nd": "1",
    "F_hole_nd": "1",
    "E_V_cm": "V/cm",
    "Jn_A_cm2": "A/cm^2",
    "Jp_A_cm2": "A/cm^2",
    "Jtotal_A_cm2": "A/cm^2",
    "I_A": "A",
    "Jterminal_A_cm2": "A/cm^2",
    "applied_voltage_V": "V",
    "iteration_count": "1",
    "runtime_s": "s",
    "converged": "1",
    "residual_norm": "1",
    "update_norm": "1",
    "alpha": "1",
    "contact_n_cm3": "cm^-3",
    "contact_p_cm3": "cm^-3",
    "contact_potential_V": "V",
    "equilibrium_np_ratio_nd": "1",
    "equilibrium_zero_current_A_cm2": "A/cm^2",
    "built_in_potential_V": "V",
}

# cm 量名 -> (SI 量名, cm→SI 数值因子)；未列出的量在 SI 下不变（§12.1 单位映射）。
_SI_QUANTITY: dict[str, tuple[str, float]] = {
    "n_cm3": ("n_m3", 1.0e6),
    "p_cm3": ("p_m3", 1.0e6),
    "ND_cm3": ("ND_m3", 1.0e6),
    "NA_cm3": ("NA_m3", 1.0e6),
    "net_doping_cm3": ("net_doping_m3", 1.0e6),
    "contact_n_cm3": ("contact_n_m3", 1.0e6),
    "contact_p_cm3": ("contact_p_m3", 1.0e6),
    "charge_C_cm3": ("charge_C_m3", 1.0e6),
    "SRH_cm3_s": ("SRH_m3_s", 1.0e6),
    "E_V_cm": ("E_V_m", 1.0e2),
    "Jn_A_cm2": ("Jn_A_m2", 1.0e4),
    "Jp_A_cm2": ("Jp_A_m2", 1.0e4),
    "Jtotal_A_cm2": ("Jtotal_A_m2", 1.0e4),
    "Jterminal_A_cm2": ("Jterminal_A_m2", 1.0e4),
    "equilibrium_zero_current_A_cm2": ("equilibrium_zero_current_A_m2", 1.0e4),
}

_SI_UNIT: dict[str, str] = {
    "cm^-3": "m^-3",
    "C/cm^3": "C/m^3",
    "cm^-3/s": "m^-3/s",
    "V/cm": "V/m",
    "A/cm^2": "A/m^2",
    "cm": "m",
}


def columns_for(unit_system: str = "cm") -> tuple[str, ...]:
    """该单位制下的表头；SI 只把 `position_cm` 改名为 `position_m`（第 12.1 节）。"""
    if unit_system == "si":
        return tuple("position_m" if column == "position_cm" else column for column in CSV_COLUMNS)
    return CSV_COLUMNS


def header_line(unit_system: str = "cm") -> str:
    """返回表头行（字符串形式，便于报告与测试直接比对）。"""
    return ",".join(columns_for(unit_system))


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


def write_rows(
    path: str | pathlib.Path,
    rows: Iterable[Sequence[str] | Mapping[str, object]],
    *,
    unit_system: str = "cm",
) -> None:
    """写 UTF-8 with BOM（utf-8-sig）的 CSV，便于中文 Excel 直接打开。

    `unit_system="si"` 时只把表头的 `position_cm` 改成 `position_m`；
    数值换算在行构造时已完成（`_scale`），这里不再改数。
    """
    with open(path, "w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(columns_for(unit_system))
        for row in rows:
            writer.writerow(format_row(**row) if isinstance(row, Mapping) else list(row))


def read_header(path: str | pathlib.Path) -> tuple[str, ...]:
    """读回表头；BOM 由 utf-8-sig 自动去掉。"""
    with open(path, "r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle)
        return tuple(next(reader, ()))


def check_header(header: Sequence[str]) -> list[str]:
    """对照冻结 schema 检查表头，返回问题列表；空列表表示通过。

    合法的表头有两种：cm 模式与 SI 模式（仅第 9 列名不同）。校验器不知道文件
    是哪种模式，所以两种都接受，其余 21 列与顺序必须完全一致。
    """
    columns = list(header)
    if columns in (list(columns_for("cm")), list(columns_for("si"))):
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


# ---------------------------------------------------------------------------
# 数据行构造：把 RunResult 摊平成 §12.1 的 long-form 行
# ---------------------------------------------------------------------------


def _blank_row(**values: object) -> dict:
    """22 列齐全的行；缺省列写空字段（`format_row` 会拒绝缺键）。"""
    row = {column: None for column in CSV_COLUMNS}
    row.update(values)
    return row


def _scale(system: str, quantity: str, value: float | None) -> float | None:
    if value is None or system == "cm":
        return value
    return value * _SI_QUANTITY.get(quantity, (quantity, 1.0))[1]


def _name(system: str, quantity: str) -> str:
    return quantity if system == "cm" else _SI_QUANTITY.get(quantity, (quantity, 1.0))[0]


def _unit(system: str, quantity: str) -> str:
    unit = QUANTITY_UNITS_CM[quantity]
    return unit if system == "cm" else _SI_UNIT.get(unit, unit)


def _node_rows(run, bias_fields, system: str) -> list[dict]:
    """11 个结点量 × M 个结点；逐点无闭式解，`theory_kind=none`。"""
    node = bias_fields.node
    quantities = {
        "phi_V": node.phi_V,
        "n_cm3": node.n_cm3,
        "p_cm3": node.p_cm3,
        "ND_cm3": node.nd_cm3,
        "NA_cm3": node.na_cm3,
        "net_doping_cm3": node.net_doping_cm3,
        "charge_C_cm3": node.charge_C_cm3,
        "SRH_cm3_s": node.srh_cm3_s,
        "F_poisson_nd": node.f_poisson_nd,
        "F_electron_nd": node.f_electron_nd,
        "F_hole_nd": node.f_hole_nd,
    }
    positions = run.device.mesh.x_nodes_cm
    factor = 1e-2 if system == "si" else 1.0
    rows = []
    for quantity, values in quantities.items():
        for index, value in enumerate(values):
            rows.append(
                _blank_row(
                    run_id=run.run_id,
                    config_sha256=run.config_sha256,
                    solver=run.config.solver.method,
                    bias_V=bias_fields.bias_V,
                    record_type="node",
                    quantity=_name(system, quantity),
                    location_kind="node",
                    index=index,
                    position_cm=float(positions[index]) * factor,
                    unit=_unit(system, quantity),
                    numeric_value=_scale(system, quantity, float(value)),
                    theory_kind="none",
                    status="PASS",
                    notes="no closed form under bias",
                )
            )
    return rows


def _edge_rows(run, bias_fields, system: str) -> list[dict]:
    edge = bias_fields.edge
    quantities = {
        "E_V_cm": edge.E_V_cm,
        "Jn_A_cm2": edge.jn_A_cm2,
        "Jp_A_cm2": edge.jp_A_cm2,
        "Jtotal_A_cm2": edge.jtotal_A_cm2,
    }
    factor = 1e-2 if system == "si" else 1.0
    rows = []
    for quantity, values in quantities.items():
        for index, value in enumerate(values):
            rows.append(
                _blank_row(
                    run_id=run.run_id,
                    config_sha256=run.config_sha256,
                    solver=run.config.solver.method,
                    bias_V=bias_fields.bias_V,
                    record_type="edge",
                    quantity=_name(system, quantity),
                    location_kind="edge",
                    index=index,
                    position_cm=float(edge.position_cm[index]) * factor,
                    unit=_unit(system, quantity),
                    numeric_value=_scale(system, quantity, float(value)),
                    theory_kind="none",
                    status="PASS",
                    notes="edge midpoint",
                )
            )
    return rows


def _sweep_rows(run, attempt, bias_fields, system: str) -> list[dict]:
    """每个尝试偏压都要有记录；失败点只写状态与原因，物理量留空（§12.1）。"""
    converged = bias_fields is not None
    status = "PASS" if converged else "FAIL"
    reason = None if converged else (attempt.result.reason if attempt.result else "NO_RESULT")
    base = dict(
        run_id=run.run_id,
        config_sha256=run.config_sha256,
        solver=run.config.solver.method,
        bias_V=float(attempt.bias_V),
        record_type="sweep",
        status=status,
    )
    rows = [
        _blank_row(**base, quantity="applied_voltage_V", unit="V",
                   numeric_value=float(attempt.bias_V), theory_value=float(attempt.bias_V),
                   theory_kind="exact", theory_abs_error=0.0, theory_rel_error=0.0, notes=reason),
        _blank_row(**base, quantity="iteration_count", unit="1",
                   numeric_value=float(attempt.result.iterations) if attempt.result else None),
        _blank_row(**base, quantity="runtime_s", unit="s", numeric_value=float(attempt.runtime_s)),
        _blank_row(**base, quantity="converged", unit="1",
                   numeric_value=1.0 if converged else 0.0),
    ]
    if not converged:
        return rows

    j_terminal = bias_fields.terminal_current_density_A_cm2()
    current = bias_fields.terminal_current_A(run.config.area_cm2)
    # 0 V 的精确理论是零电流：写绝对误差，相对误差留空（§12.1 理论/误差规则）。
    zero_bias = abs(float(attempt.bias_V)) < 1e-12
    rows.append(
        _blank_row(
            **base,
            quantity=_name(system, "Jterminal_A_cm2"),
            unit=_unit(system, "Jterminal_A_cm2"),
            numeric_value=_scale(system, "Jterminal_A_cm2", j_terminal),
            theory_value=0.0 if zero_bias else None,
            theory_kind="exact" if zero_bias else "none",
            theory_abs_error=abs(_scale(system, "Jterminal_A_cm2", j_terminal)) if zero_bias else None,
            theory_rel_error=None,
            notes="zero theory: relative error undefined" if zero_bias else None,
        )
    )
    rows.append(
        _blank_row(
            **base,
            quantity="I_A",
            unit="A",
            numeric_value=current,
            theory_value=0.0 if zero_bias else None,
            theory_kind="exact" if zero_bias else "none",
            theory_abs_error=abs(current) if zero_bias else None,
            theory_rel_error=None,
            notes="I_A = area_cm2 * Jterminal_A_cm2",
        )
    )
    return rows


def _iteration_rows(run, attempt) -> list[dict]:
    history = attempt.result.iteration_history if attempt.result else ()
    rows = []
    for record in history:
        for quantity in ("residual_norm", "update_norm", "alpha"):
            value = record.get(quantity)
            rows.append(
                _blank_row(
                    run_id=run.run_id,
                    config_sha256=run.config_sha256,
                    solver=run.config.solver.method,
                    bias_V=float(attempt.bias_V),
                    record_type="iteration",
                    quantity=quantity,
                    unit="1",
                    iteration=int(record.get("iteration", 0)),
                    numeric_value=None if value is None else float(value),
                    status="PASS" if record.get("status") == "ACCEPTED" else "FAIL",
                    notes=record.get("status"),
                )
            )
    return rows


def _benchmark_rows(run, system: str) -> list[dict]:
    """只依赖器件的精确量：接触浓度/接触势、内建电势、质量作用与零电流基线。

    `contact_n_cm3`/`contact_p_cm3`/`contact_potential_V` 是**按端**的量：
    同一 quantity 出现两次，用 `index`（0 或 M-1）与位置区分左右接触（§12.1 的 location 列）。
    """
    from minitcad.simulation.diagnostics import benchmark_fields

    values = benchmark_fields(run.config, run.device)
    m = run.device.mesh.node_count
    x = run.device.mesh.x_nodes_cm
    scale_position = 1e-2 if system == "si" else 1.0
    base = dict(
        run_id=run.run_id,
        config_sha256=run.config_sha256,
        solver=run.config.solver.method,
        bias_V=0.0,
        record_type="benchmark",
        status="PASS",
    )

    rows = []
    contacts = (
        (0, "left", x[0], values["contact_n_left_cm3"], values["contact_p_left_cm3"],
         values["contact_potential_left_V"]),
        (m - 1, "right", x[-1], values["contact_n_right_cm3"], values["contact_p_right_cm3"],
         values["contact_potential_right_V"]),
    )
    for index, side, position_cm, n_cm3, p_cm3, phi_V in contacts:
        location = dict(location_kind="node", index=index, position_cm=float(position_cm) * scale_position)
        for quantity, value in (("contact_n_cm3", n_cm3), ("contact_p_cm3", p_cm3)):
            # 接触浓度是精确理论值：数值列写同一值，误差为 0。
            rows.append(
                _blank_row(
                    **base,
                    **location,
                    quantity=_name(system, quantity),
                    unit=_unit(system, quantity),
                    numeric_value=_scale(system, quantity, value),
                    theory_value=_scale(system, quantity, value),
                    theory_kind="exact",
                    theory_abs_error=0.0,
                    theory_rel_error=0.0,
                    notes=f"{side} contact: n=e^a, p=e^-a, np=1",
                )
            )
        rows.append(
            _blank_row(
                **base,
                **location,
                quantity="contact_potential_V",
                unit="V",
                numeric_value=phi_V,
                theory_value=phi_V,
                theory_kind="exact",
                theory_abs_error=0.0,
                theory_rel_error=0.0,
                notes=f"{side} contact: a(C)=asinh(C/2)",
            )
        )

    for quantity, value, note in (
        ("built_in_potential_V", values["built_in_potential_V"], "V_bi = V_T*(a_R - a_L)"),
    ):
        rows.append(
            _blank_row(
                **base,
                quantity=quantity,
                unit="V",
                numeric_value=value,
                theory_value=value,
                theory_kind="exact",
                theory_abs_error=0.0,
                theory_rel_error=0.0,
                notes=note,
            )
        )
    rows.append(
        _blank_row(
            **base,
            quantity="equilibrium_np_ratio_nd",
            unit="1",
            numeric_value=values["equilibrium_np_ratio_nd"],
            theory_value=1.0,
            theory_kind="exact",
            theory_abs_error=0.0,
            theory_rel_error=0.0,
        )
    )
    rows.append(
        _blank_row(
            **base,
            quantity=_name(system, "equilibrium_zero_current_A_cm2"),
            unit=_unit(system, "equilibrium_zero_current_A_cm2"),
            numeric_value=values["equilibrium_zero_current_A_cm2"],
            theory_value=0.0,
            theory_kind="exact",
            theory_abs_error=0.0,
            theory_rel_error=None,
            notes="equilibrium: zero current is exact",
        )
    )
    return rows


def build_rows(run) -> list[dict]:
    """把 `RunResult` 摊平成 long-form 行；`run` 需含 `config/device/attempts`。"""
    from minitcad.simulation.diagnostics import fields_for_attempt

    system = run.config.output.unit_system
    rows: list[dict] = []
    for attempt in run.attempts:
        bias_fields = fields_for_attempt(run.config, run.device, attempt)
        rows.extend(_sweep_rows(run, attempt, bias_fields, system))
        rows.extend(_iteration_rows(run, attempt))
        if bias_fields is not None:
            rows.extend(_node_rows(run, bias_fields, system))
            rows.extend(_edge_rows(run, bias_fields, system))
    rows.extend(_benchmark_rows(run, system))
    return rows


def write_comparison_csv(run, path: str | pathlib.Path) -> pathlib.Path:
    """写 `comparison.csv` 并返回路径；表头/编码/换行由 `write_rows` 统一。"""
    target = pathlib.Path(path)
    write_rows(target, build_rows(run), unit_system=run.config.output.unit_system)
    return target
