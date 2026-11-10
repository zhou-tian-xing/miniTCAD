"""两种求解器共用的收敛判据：scaled residual 与 scaled update（第 8.2 节）。

两个量都除以各自的量级，因此可以和 `tol_residual`、`tol_update` 直接比较，
不受掺杂浓度（可达 `1e6` 量级）与电流量级的影响。

Gummel 与 Newton 必须共用本模块，否则“同一物理状态在两种方法下收敛性不同”
会把求解器的差异和判据的差异混在一起，Gate 2 的双求解器对照就失去意义。

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
from numpy.typing import NDArray

from minitcad.core.contact import contact_targets_qf
from minitcad.core.device import NormalizedDevice
from minitcad.core.state import QFState, State
from minitcad.discretization.residual import ResidualEvaluation


def residual_scales_qf(
    qf: QFState,
    device: NormalizedDevice,
    bias_V: float,
    evaluation: ResidualEvaluation,
    physical: State | None = None,
) -> NDArray[np.float64]:
    """残差的行尺度，展平成长度 `3M` 的向量，行序与 `evaluation.f` 一致。

    Newton 需要逐行尺度而不只是它们的最大值：它把 Jacobian 的每一行除以该行尺度
    再解线性方程（行缩放），否则 Poisson 行（量级 `1e6`）会把连续性行（量级 `1`）
    的修正压在求解精度以下。`scaled_residual_qf` 只是对同一尺度的 `max` 归约，
    两者必须共用本函数，避免判定用的尺度与求解用的尺度悄悄分叉。

    `physical` 缺省时由 `qf.to_physical()` 重建；求解器在热路径上已经持有物理状态
    时传入，可以省掉一次 `exp`。
    """
    if physical is None:
        physical = qf.to_physical()
    m = qf.node_count
    if m != device.mesh.node_count or m != physical.node_count:
        raise ValueError("state/device node_count mismatch")
    if (
        evaluation.f.shape != (3 * m,)
        or evaluation.recombination_nd.shape != (m,)
        or evaluation.fluxes.jn_edges.shape != (m - 1,)
        or evaluation.fluxes.jp_edges.shape != (m - 1,)
    ):
        raise ValueError("invalid QF evaluation shapes")

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

    if not np.isfinite(scale).all() or np.any(scale <= 0.0):
        raise FloatingPointError("NONFINITE_RESIDUAL_SCALE")
    return scale.reshape(3 * m)


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
    scale = residual_scales_qf(qf, device, bias_V, evaluation, physical)
    value = float(np.max(np.abs(evaluation.f) / scale))
    if not np.isfinite(value):
        raise FloatingPointError("NONFINITE_SCALED_RESIDUAL")
    return value


def scaled_update_physical(
    old: State,
    new: State,
    n_floor: float = 1e-15,
) -> float:
    """三个物理量的相对变化的最大值；载流子分母用 `n_floor` 兜底。

    名字里的 `physical` 是必要的：判据作用在 `State(u,n,p)` 上，而不是 QF 迭代变量
    `(u, eta_n, eta_p)` 上。对 `eta` 直接取相对变化没有物理含义（`eta` 的零点由锚点
    `ref_n/ref_p` 决定，换一个锚点就换一个数值），所以统一先 `to_physical()`。
    """
    if old.node_count != new.node_count:
        raise ValueError("old/new node_count mismatch")
    if not np.isfinite(n_floor) or n_floor <= 0.0:
        raise ValueError("n_floor must be finite and positive")
    m = old.node_count
    su = np.maximum.reduce((np.ones(m), np.abs(old.u), np.abs(new.u)))
    sn = np.maximum.reduce((np.full(m, n_floor), old.n, new.n))
    sp = np.maximum.reduce((np.full(m, n_floor), old.p, new.p))
    value = float(
        max(
            np.max(np.abs(new.u - old.u) / su),
            np.max(np.abs(new.n - old.n) / sn),
            np.max(np.abs(new.p - old.p) / sp),
        )
    )
    if not np.isfinite(value):
        raise FloatingPointError("NONFINITE_PHYSICAL_UPDATE")
    return value


def scaled_update(old: State, new: State, floor: float) -> float:
    """`scaled_update_physical` 的旧名（W4 起 Gummel 在用），保留以免下游导入失效。"""
    return scaled_update_physical(old, new, n_floor=floor)


def is_converged(residual_norm: float, update_norm: float, tol_residual: float, tol_update: float) -> bool:
    """两个判据必须同时满足（第 8.2 节）：只看残差会把停滞的迭代误判为收敛。"""
    return residual_norm <= tol_residual and update_norm <= tol_update
