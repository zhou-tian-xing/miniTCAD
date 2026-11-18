"""两种求解器共用的收敛判据（第 8.2 节）：scaled residual 与 scaled update。

两个量都除以自身的量级，才能和 `tol_residual`、`tol_update` 直接比较——掺杂浓度到
`1e6` 量级，绝对残差没有意义。Gummel 与 Newton 共用本模块，否则双求解器对照会把
“方法差异”和“判据差异”混在一起。

* scaled residual：三类行分别除 `max(1,|C̄|,n,p)`、`max(1,|j|/w,|r|)`、`max(1,|target|)`；
* scaled update：`Δu`、`Δn`、`Δp` 的相对变化取最大值，载流子分母用 `n_floor` 兜底，
  防止 `n` 极小时把相对步长放大成假的“未收敛”。
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
    """残差的行尺度，展平成长度 `3M`，行序与 `evaluation.f` 一致。

    Newton 的逐行缩放与判据的 `max` 归约都取自这里，避免两边悄悄分叉。
    `physical` 缺省时由 `qf.to_physical()` 重建，热路径上传入可省一次 `exp`。
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
    """`max_i |F_i| / scale_i`；`physical` 缺省时由 `qf.to_physical()` 重建。"""
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
    """`Δu`、`Δn`、`Δp` 的相对变化取最大值；载流子分母用 `n_floor` 兜底。

    名字里的 `physical` 是必要的：判据作用在 `State(u,n,p)` 上。`η` 的零点由锚点决定，
    对它取相对变化没有物理含义。
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
    """`scaled_update_physical` 的旧名，保留以免下游导入失效。"""
    return scaled_update_physical(old, new, n_floor=floor)


def is_converged(residual_norm: float, update_norm: float, tol_residual: float, tol_update: float) -> bool:
    """两个判据必须同时满足：只看残差会把停滞的迭代误判为收敛。"""
    return residual_norm <= tol_residual and update_norm <= tol_update
