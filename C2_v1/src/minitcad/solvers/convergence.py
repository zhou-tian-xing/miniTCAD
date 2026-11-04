"""Gummel 的收敛判据：scaled residual 与 scaled update（第 8.2 节）。

两个量都除以各自的量级，因此可以和 `tol_residual`、`tol_update` 直接比较，
不受掺杂浓度（可达 `1e6` 量级）与电流量级的影响。

**scaled residual**：`F(X)` 的三类行分别除以该行在物理上“不可能更小”的尺度——

* Poisson 行除 `max(1, |C̄_i|, n_i, p_i)`：电荷密度的量级；
* 连续性行除 `max(1, |j|/w, |r|)`：通量散度与复合率的量级；
* 接触行除 `max(1, |target|)`：Dirichlet 目标值的量级。

**scaled update**：`Δu / max(1,|u_old|,|u_new|)`、`Δn / max(floor, n_old, n_new)`、
`Δp / max(floor, p_old, p_new)` 的最大值。载流子用 `floor` 防止 `n` 极小时把
相对步长放大成假的“未收敛”。

`floor` 取 `n_floor`（默认 `1e-15`），与 `NormalizedDevice` 的量级一致：
`1e-15` 无量纲对应 `1e-5 cm^-3`，低于任何有物理意义的浓度。
"""

from __future__ import annotations

import numpy as np

from minitcad.core.contact import contact_targets_qf
from minitcad.core.device import NormalizedDevice
from minitcad.core.state import QFState, State
from minitcad.discretization.residual import ResidualEvaluation


def scaled_residual_qf(
    qf: QFState,
    device: NormalizedDevice,
    bias_V: float,
    evaluation: ResidualEvaluation,
    physical: State | None = None,
) -> float:
    """`max_i |F_i| / scale_i`；`physical` 缺省时由 `qf.to_physical()` 重建。

    `scale` 的三类构造见模块文档。`F` 与 `scale` 一定是同形状 `(M,3)`：
    `evaluation.f` 是结点交错的一维向量，先 `reshape`。
    """
    if physical is None:
        physical = qf.to_physical()
    m = qf.node_count
    f = evaluation.f.reshape(m, 3)
    flux = evaluation.fluxes
    rate = evaluation.recombination_nd
    w = device.mesh.w_interior_nd

    interior = slice(1, -1)
    scale = np.empty((m, 3), dtype=np.float64)
    scale[interior, 0] = np.maximum.reduce(
        (np.ones(m - 2), np.abs(device.cbar_nd), physical.n[interior], physical.p[interior])
    )
    # 边通量有 M-1 条，内部结点 i 用到 j_{i-1/2} 与 j_{i+1/2}，即数组下标 i-1 与 i。
    scale[interior, 1] = np.maximum.reduce(
        (np.ones(m - 2), np.abs(flux.jn_edges[1:]) / w, np.abs(flux.jn_edges[:-1]) / w, np.abs(rate[interior]))
    )
    scale[interior, 2] = np.maximum.reduce(
        (np.ones(m - 2), np.abs(flux.jp_edges[1:]) / w, np.abs(flux.jp_edges[:-1]) / w, np.abs(rate[interior]))
    )

    targets = contact_targets_qf(
        device=device, bias_V=bias_V, ref_n=qf.ref_n, ref_p=qf.ref_p
    )
    scale[0] = np.maximum(1.0, np.abs(targets[:3]))
    scale[-1] = np.maximum(1.0, np.abs(targets[3:]))

    value = float(np.max(np.abs(f) / scale))
    if not np.isfinite(value):
        raise FloatingPointError("NONFINITE_SCALED_RESIDUAL")
    return value


def scaled_update(old: State, new: State, floor: float) -> float:
    """三个物理量的相对变化的最大值；载流子分母用 `floor` 兜底。"""
    su = np.maximum.reduce((np.ones(old.node_count), np.abs(old.u), np.abs(new.u)))
    sn = np.maximum.reduce((np.full(old.node_count, floor), old.n, new.n))
    sp = np.maximum.reduce((np.full(old.node_count, floor), old.p, new.p))
    return float(
        max(
            np.max(np.abs(new.u - old.u) / su),
            np.max(np.abs(new.n - old.n) / sn),
            np.max(np.abs(new.p - old.p) / sp),
        )
    )


def is_converged(residual_norm: float, update_norm: float, tol_residual: float, tol_update: float) -> bool:
    """两个判据必须同时满足（第 8.2 节）：只看残差会把停滞的迭代误判为收敛。"""
    return residual_norm <= tol_residual and update_norm <= tol_update
