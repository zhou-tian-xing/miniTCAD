"""`comparison.png`（第 12.2 节）：一张多面板图，数值 / 理论 / 误差分开呈现。

面板与数据都来自 `comparison.csv` 的同一份展开结果（`simulation.diagnostics`），
所以图和表不会各写一套物理量。没有理论曲线的偏压只画数值，并在图注里写明
"Theory: N/A"，不臆造理论（§12.2 的边界）。

图只在内存里用 `Agg` 后端绘制，不依赖桌面环境。
"""

from __future__ import annotations

import pathlib

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402  (backend 必须先设置)
import numpy as np  # noqa: E402

from minitcad.simulation.diagnostics import fields_for_attempt  # noqa: E402

ORIENTATION_NOTE = "P(left, x=0) -> N(right, x=L) | V_app = V_left - V_right | +V = forward | right = 0 V"


def _accepted_fields(run) -> list:
    """所有已收敛偏压的展开结果，按偏压升序。"""
    fields = []
    for attempt in run.successful:
        expanded = fields_for_attempt(run.config, run.device, attempt)
        if expanded is not None:
            fields.append(expanded)
    return fields


def write_comparison_png(run, path: str | pathlib.Path) -> pathlib.Path:
    """画 6 个面板并写 PNG；返回路径。"""
    target = pathlib.Path(path)
    fields = _accepted_fields(run)
    unit_cm = run.config.output.unit_system == "cm"
    x_unit = "x [cm]" if unit_cm else "x [m]"
    j_unit = "J [A/cm$^2$]" if unit_cm else "J [A/m$^2$]"
    n_unit = "n, p [cm$^{-3}$]" if unit_cm else "n, p [m$^{-3}$]"

    fig, axes = plt.subplots(3, 2, figsize=(13, 11), constrained_layout=True)
    fig.suptitle(
        f"miniTCAD — {run.config.solver.method} — {run.decision}"
        f"{'' if run.accepted else ' (see convergence.md)'}\n{ORIENTATION_NOTE}",
        fontsize=10,
    )

    # (1) I-V：只画已接受偏压，失败的偏压点不插值、不连线。
    ax = axes[0, 0]
    if fields:
        voltages = [f.bias_V for f in fields]
        currents = [f.terminal_current_A(run.config.area_cm2) for f in fields]
        ax.plot(voltages, currents, "o-", lw=1.5)
    ax.set_xlabel("V_app [V]")
    ax.set_ylabel("I [A]")
    ax.set_title("I-V (accepted biases only)")
    ax.grid(alpha=0.3)

    if not fields:
        for ax in axes.flat[1:]:
            ax.axis("off")
        target.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(target, dpi=150)
        plt.close(fig)
        return target

    representative = [fields[0], fields[len(fields) // 2], fields[-1]]
    scale_position = 1e-2 if not unit_cm else 1.0

    # (2) phi(x) 若干偏压
    ax = axes[0, 1]
    for f in representative:
        ax.plot(run.device.mesh.x_nodes_cm * scale_position, f.node.phi_V, label=f"{f.bias_V:.3g} V")
    ax.set_xlabel(x_unit)
    ax.set_ylabel("phi [V]")
    ax.set_title("Potential profiles")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # (3) n, p 对数轴
    ax = axes[1, 0]
    x = run.device.mesh.x_nodes_cm * scale_position
    for f in [fields[-1]]:
        ax.semilogy(x, f.node.n_cm3, label="n")
        ax.semilogy(x, f.node.p_cm3, label="p")
    ax.set_xlabel(x_unit)
    ax.set_ylabel(n_unit)
    ax.set_title(f"Carriers at {fields[-1].bias_V:.3g} V")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")

    # (4) E(x) 与 J(x)，总电流沿空间应接近常数
    ax = axes[1, 1]
    edge = fields[-1].edge
    axe = ax.twinx()
    ax.plot(edge.position_cm * scale_position, edge.jn_A_cm2, label="$J_n$")
    ax.plot(edge.position_cm * scale_position, edge.jp_A_cm2, label="$J_p$")
    ax.plot(edge.position_cm * scale_position, edge.jtotal_A_cm2, "--", label="$J_{total}$")
    axe.plot(edge.position_cm * scale_position, edge.E_V_cm, ":", color="0.4", label="$E$")
    ax.set_xlabel(x_unit)
    ax.set_ylabel(j_unit)
    axe.set_ylabel("E [V/cm]" if unit_cm else "E [V/m]")
    lines = ax.get_lines() + axe.get_lines()
    ax.legend(lines, [line.get_label() for line in lines], fontsize=8)
    ax.set_title("Current density and field")
    ax.grid(alpha=0.3)

    # (5) 迭代：残差 vs 迭代步，按偏压区分
    ax = axes[2, 0]
    for f in fields:
        records = [r for r in f.iteration_history if r.get("residual_norm") is not None
                   and np.isfinite(r["residual_norm"])]
        if records:
            ax.semilogy([r["iteration"] for r in records], [r["residual_norm"] for r in records],
                        "o-", ms=3, lw=1, label=f"{f.bias_V:.3g} V")
    ax.set_xlabel("outer iteration")
    ax.set_ylabel("scaled residual")
    ax.set_title("Convergence history")
    if ax.get_legend_handles_labels()[0]:
        ax.legend(fontsize=7)
    ax.grid(alpha=0.3, which="both")

    # (6) 验收摘要：每偏压的迭代与终端电流，加上报告第 7 节的指标（§12.2 第 6 条）
    from minitcad.io.report import CURRENT_SPREAD_TOL, NP_TOL, ZERO_CURRENT_TOL, validation_metrics

    ax = axes[2, 1]
    ax.axis("off")
    metrics = validation_metrics(run)
    rows = [("check", "measured", "threshold", "result")]
    if metrics:
        rows.extend(
            [
                ("|np-1| @ 0 V", f"{metrics['np_error']:.2e}", f"{NP_TOL:g}",
                 "PASS" if metrics["np_error"] <= NP_TOL else "FAIL"),
                ("|J| @ 0 V [A/cm2]", f"{metrics['zero_current']:.2e}", f"{ZERO_CURRENT_TOL:g}",
                 "PASS" if metrics["zero_current"] <= ZERO_CURRENT_TOL else "FAIL"),
                ("current spread", f"{metrics['current_spread']:.2e}", f"{CURRENT_SPREAD_TOL:g}",
                 "PASS" if metrics["current_spread"] <= CURRENT_SPREAD_TOL else "FAIL"),
                ("accepted biases", f"{len(fields)}/{len(run.attempts)}", "-",
                 "PASS" if run.accepted else "FAIL"),
            ]
        )
    else:
        rows.append(("no accepted point", "N/A", "-", "FAIL"))
    table = ax.table(cellText=rows, loc="center", cellLoc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1, 1.5)
    ax.set_title("Acceptance summary (convergence.md section 7)", fontsize=10)

    target.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(target, dpi=150)
    plt.close(fig)
    return target
