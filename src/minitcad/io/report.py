"""`convergence.md`（第 12.3 节）：十个固定章节的收敛与验收报告。

报告不复制 CSV 的逐点数据：逐次迭代在 `comparison.csv` 的 `iteration` 行里，
这里只放每偏压汇总、验收指标与复现信息。所有 `PASS` 都由本文件第 7 节列出的
实测数值支撑；未通过的项直说未通过，不写“基本通过”。

指标口径（都用无量纲残差/物理量换算后的量，见 `simulation.diagnostics`）：

* `np-1`：0 V 平衡点的质量作用偏差；`np = n_nd * p_nd`，理论值 1。
* `zero_current`：0 V 的终端总电流密度 `|J_n+J_p|`，理论值 0。
* `current_spread`：同一收敛偏压下 `max|J_total - median(J_total)| / max(1, |median|)`；
  跨结点的总电流应沿空间恒定，这是离散格式守恒性的直接检查（T09）。
* `contact_error`：两端 `u` 与 `contact_targets` 的差、`n,p` 与 `e^{±a}` 的差的最大值。
"""

from __future__ import annotations

import platform
import pathlib
import sys

import numpy as np
import yaml

from minitcad.io.config import config_payload
from minitcad.simulation.diagnostics import fields_for_attempt

NP_TOL = 1e-6
ZERO_CURRENT_TOL = 1e-6
CONTACT_TOL = 1e-10
CURRENT_SPREAD_TOL = 1e-3


def validation_metrics(run) -> dict[str, float]:
    """从已收敛状态直接算验收指标；没有收敛点就返回空字典。"""
    fields = [fields_for_attempt(run.config, run.device, attempt) for attempt in run.successful]
    fields = [item for item in fields if item is not None]
    if not fields:
        return {}

    zero = next((item for item in fields if abs(item.bias_V) < 1e-12), None)
    np_error = float("nan")
    zero_current = float("nan")
    if zero is not None:
        np_error = float(np.max(np.abs(zero.node.n_cm3 * zero.node.p_cm3 / run.device.ni_cm3**2 - 1.0)))
        zero_current = abs(zero.terminal_current_density_A_cm2())

    # 总电流沿空间恒定是离散守恒性的检查：用 |J| 的最大值归一化。
    # 0 V 时 J 恒为零，此时 ptp 只是舍入噪声，用 1e-9 A/cm^2 的绝对下限吸收，
    # 否则会把它报成“守恒失败”。
    spread = 0.0
    for item in fields:
        j = item.edge.jtotal_A_cm2
        denominator = max(float(np.max(np.abs(j))), 1e-9)
        spread = max(spread, float(np.ptp(j)) / denominator)

    targets = run.device.contact_targets(0.0)
    state = run.successful[0].result.state
    contact_error = max(
        abs(state.u[0] - targets[0]),
        abs(state.u[-1] - targets[3]),
        abs(np.log(state.n[0]) - np.log(targets[1])),
        abs(-np.log(state.p[0]) - (-np.log(targets[2]))),
        abs(np.log(state.n[-1]) - np.log(targets[4])),
        abs(-np.log(state.p[-1]) - (-np.log(targets[5]))),
    )
    residual = max(attempt.result.residual_norm for attempt in run.successful)
    update = max(attempt.result.update_norm for attempt in run.successful)
    return {
        "np_error": np_error,
        "zero_current": zero_current,
        "current_spread": spread,
        "contact_error": float(contact_error),
        "max_residual": float(residual),
        "max_update": float(update),
    }


def write_convergence_md(run, path: str | pathlib.Path) -> pathlib.Path:
    """写 `convergence.md` 并返回路径。"""
    target = pathlib.Path(path)
    device = run.device
    scales = device.scales
    metrics = validation_metrics(run)
    unit_system = run.config.output.unit_system

    lines: list[str] = []
    add = lines.append

    add(f"# convergence.md — {run.run_id}")
    add("")
    add("## 1. Run summary")
    add("")
    add(f"- run_id: `{run.run_id}`")
    add(f"- config_sha256: `{run.config_sha256}`")
    add(f"- platform: {platform.platform()}")
    add(f"- python: {sys.version.split()[0]}")
    add(f"- decision: **{run.decision}**")
    add(f"- attempts: {len(run.attempts)} (accepted: {len(run.successful)})")
    if run.failure_reason is not None:
        add(f"- failure_reason: {run.failure_reason}")
    add("")

    add("## 2. Orientation")
    add("")
    add("P(left, x=0) -> N(right, x=L); V_app = V_left - V_right; +V = forward bias; right = 0 V")
    add("")

    add("## 3. Physics / units")
    add("")
    add(f"- input_unit_system: `{run.config.input_unit_system}`")
    add("- core_unit_system: `dimensionless`")
    add(f"- output_unit_system: `{unit_system}`")
    add(f"- material: {device.material.name}, T = {device.temperature_k} K, V_T = {scales.v0_v:.9g} V")
    add(f"- parameter_policy: {device.material_parameter_policy}")
    add(
        f"- scales (cm): eps_s = {device.material.epsilon_f_per_cm:.9g} F/cm, "
        f"n_i = {scales.n0_cm3:.9g} cm^-3, mu_0 = {scales.mu0_cm2_per_vs:.9g} cm^2/(V s)"
    )
    add(
        f"- scales (cm): L_0 = {scales.l0_cm:.9g} cm, D_0 = {scales.d0_cm2_s:.9g} cm^2/s, "
        f"t_0 = {scales.t0_s:.9g} s, J_0 = {scales.j0_a_cm2:.9g} A/cm^2, R_0 = {scales.r0_cm3_s:.9g} cm^-3/s"
    )
    if unit_system == "si":
        add(
            f"- scales (si, writer only): L_0 = {scales.l0_cm * 1e-2:.9g} m, "
            f"n_i = {scales.n0_cm3 * 1e6:.9g} m^-3, J_0 = {scales.j0_a_cm2 * 1e4:.9g} A/m^2"
        )
    add("")

    add("## 4. Solver settings")
    add("")
    solver = run.config.solver
    add(
        f"- method: {solver.method}; jacobian: {solver.jacobian}; max_iter: {solver.max_iter}; "
        f"alpha_min: {solver.alpha_min}; n_floor: {solver.n_floor}"
    )
    add("")

    add("## 5. Tolerances")
    add("")
    add(f"- tol_residual: {solver.tol_residual}; tol_update: {solver.tol_update} (both scaled)")
    add(
        f"- acceptance: |np-1| <= {NP_TOL}; |J(0 V)| <= {ZERO_CURRENT_TOL}; "
        f"contact error <= {CONTACT_TOL}; current spread <= {CURRENT_SPREAD_TOL}"
    )
    add("")

    add("## 6. Bias table")
    add("")
    add("| V_app [V] | status | iterations | final Rnorm | final Dnorm | runtime [s] | reason |")
    add("|---|---|---|---|---|---|---|")
    for attempt in run.attempts:
        result = attempt.result
        add(
            f"| {attempt.bias_V:.6g} | {attempt.status} | "
            f"{result.iterations if result else '-'} | "
            f"{result.residual_norm:.3e} | {result.update_norm:.3e} | "
            f"{attempt.runtime_s:.3f} | {result.reason if result else 'NO_RESULT'} |"
            if result
            else f"| {attempt.bias_V:.6g} | {attempt.status} | - | - | - | {attempt.runtime_s:.3f} | NO_RESULT |"
        )
    add("")

    add("## 7. Validation")
    add("")
    if metrics:
        add("| item | measured | threshold | result |")
        add("|---|---|---|---|")
        add(
            f"| mass action at 0 V (max abs np-1) | {metrics['np_error']:.3e} | {NP_TOL:g} | "
            f"{'PASS' if metrics['np_error'] <= NP_TOL else 'FAIL'} |"
        )
        add(
            f"| zero current at 0 V (A/cm^2) | {metrics['zero_current']:.3e} | {ZERO_CURRENT_TOL:g} | "
            f"{'PASS' if metrics['zero_current'] <= ZERO_CURRENT_TOL else 'FAIL'} |"
        )
        add(
            f"| contact neutrality (max abs) | {metrics['contact_error']:.3e} | {CONTACT_TOL:g} | "
            f"{'PASS' if metrics['contact_error'] <= CONTACT_TOL else 'FAIL'} |"
        )
        add(
            f"| current conservation (relative spread) | {metrics['current_spread']:.3e} | "
            f"{CURRENT_SPREAD_TOL:g} | {'PASS' if metrics['current_spread'] <= CURRENT_SPREAD_TOL else 'FAIL'} |"
        )
        add(
            f"| max scaled residual | {metrics['max_residual']:.3e} | {solver.tol_residual:g} | "
            f"{'PASS' if metrics['max_residual'] <= solver.tol_residual else 'FAIL'} |"
        )
        add(
            f"| max scaled update | {metrics['max_update']:.3e} | {solver.tol_update:g} | "
            f"{'PASS' if metrics['max_update'] <= solver.tol_update else 'FAIL'} |"
        )
        add("")
        add(
            "Bernoulli / SRH / SG / Poisson 的独立 oracle（T01–T07）由 "
            "`python -m pytest tests/test_oracle_*.py -q` 覆盖，不在本报告重复。"
        )
    else:
        add("No accepted bias point; validation metrics are undefined (N_A).")
    add("")

    add("## 8. Error policy")
    add("")
    add("- `theory_kind=exact` 仅用于接触浓度/接触势、内建电势、0 V 的质量作用与零电流。")
    add("- 有偏压的逐点 `phi, n, p` 无通用闭式解：`theory_kind=none`，理论列留空。")
    add("- 零理论值（0 V 电流）写绝对误差，相对误差留空（`zero theory: relative error undefined`）。")
    add("")

    add("## 9. Reproducibility")
    add("")
    add("```bash")
    add(f"python -m minitcad.cli solve <config.yaml> --out {run.config.output.directory}")
    add("```")
    add("")
    add("```yaml")
    for line in yaml.safe_dump(config_payload(run.config), sort_keys=False, allow_unicode=True).splitlines():
        add(line)
    add("```")
    add("")

    add("## 10. Final decision")
    add("")
    add(f"**{run.decision}**")
    if run.failure_reason is not None:
        add(f"- reason: {run.failure_reason}")
    add("")

    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return target
