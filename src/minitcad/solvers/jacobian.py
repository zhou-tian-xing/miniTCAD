"""Newton 用的 Jacobian：九色稀疏中心差分。

用差分而不是解析求导：残差里的 `expm1`、Bernoulli 和分段 SG 通量都要与
`discretization/` 逐字同步，手写导数一旦落后于残差，收敛性会悄悄退化。差分永远同源。

九色：控制体离散后第 `r` 行只依赖结点 `k-1, k, k+1`，按 `k % 3` 分色后同色结点相距
≥3 个网格间距，扰动互不干扰。于是 `3` 个分量 × `3` 种颜色 = 9 次整向量求值就拿到
全部 `3M` 列，而不是 `6M` 次。

每列只取 `k-1..k+1` 这 9 行：整向量差分的其余行理论上是 0，截掉可以去掉抵消噪声。
"""

from __future__ import annotations

import numpy as np
from scipy.sparse import csc_matrix, coo_matrix

from minitcad.core.device import NormalizedDevice
from minitcad.core.state import QFState
from minitcad.discretization.residual import residual_qf

# 扰动点落进指数溢出区是正常的，按“该方向不可用”处理。
_PERTURBATION_FAILURES = (ValueError, FloatingPointError, OverflowError, ArithmeticError)


def jacobian_qf(
    state: QFState,
    device: NormalizedDevice,
    bias_V: float,
    *,
    rel_step: float = 1e-5,
) -> csc_matrix:
    """返回 `dF/dz`，`z = QFState.pack_qf()`，形状 `(3M, 3M)` 稀疏 CSC。

    步长 `delta_j = rel_step * max(1, |z_j|)`。中心差分的截断误差 `O(delta^2)`、
    舍入误差 `O(eps/delta)` 在 `delta ~ 1e-7` 处平衡；取 `1e-5` 偏向抑制舍入，
    Newton 只需要方向够准。
    """
    m = state.node_count
    if m != device.mesh.node_count:
        raise ValueError("QFState/device node_count mismatch")
    if isinstance(rel_step, (bool, np.bool_)) or not np.isfinite(rel_step) or rel_step <= 0.0:
        raise ValueError("rel_step must be positive and finite")

    z = state.pack_qf()
    n_var = z.size
    f0 = residual_qf(state=state, device=device, bias_V=bias_V)
    if f0.shape != (n_var,) or not np.isfinite(f0).all():
        raise FloatingPointError("invalid baseline QF residual")

    steps = rel_step * np.maximum(1.0, np.abs(z))
    if not np.isfinite(steps).all() or np.any(steps <= 0.0):
        raise ValueError("invalid Jacobian steps")

    rows: list[int] = []
    cols: list[int] = []
    data: list[float] = []

    def evaluate(vector: np.ndarray) -> np.ndarray:
        qf = QFState.unpack_qf(
            z=vector, node_count=m, ref_n=state.ref_n, ref_p=state.ref_p
        )
        out = residual_qf(state=qf, device=device, bias_V=bias_V)
        if out.shape != (n_var,) or not np.isfinite(out).all():
            raise FloatingPointError("nonfinite perturbed residual")
        return out

    def local_derivative(column: int, step: float) -> np.ndarray:
        """单列差分；中心差分不可用时退化为单边，步长减半重试。"""
        for _ in range(12):
            if z[column] + step == z[column] or z[column] - step == z[column]:
                # 步长小到加不上去，再减半也是同一个数。
                break
            plus = None
            minus = None
            zp = z.copy()
            zm = z.copy()
            zp[column] += step
            zm[column] -= step
            try:
                plus = evaluate(vector=zp)
            except _PERTURBATION_FAILURES:
                pass
            try:
                minus = evaluate(vector=zm)
            except _PERTURBATION_FAILURES:
                pass
            if plus is not None and minus is not None:
                return (plus - minus) / (2.0 * step)
            if plus is not None:
                return (plus - f0) / step
            if minus is not None:
                return (f0 - minus) / step
            step *= 0.5
        raise FloatingPointError(f"cannot differentiate QF coordinate {column}")

    for component in range(3):
        for node_color in range(3):
            columns = 3 * np.arange(node_color, m, 3) + component
            if columns.size == 0:
                continue
            zp = z.copy()
            zm = z.copy()
            zp[columns] += steps[columns]
            zm[columns] -= steps[columns]
            try:
                fp = evaluate(zp)
                fm = evaluate(zm)
            except _PERTURBATION_FAILURES:
                fp = None
                fm = None
            for col in columns:
                node = col // 3
                start = 3 * max(0, node - 1)
                end = 3 * min(m, node + 2)
                if fp is not None and fm is not None:
                    local = (fp[start:end] - fm[start:end]) / (2.0 * steps[col])
                else:
                    local = local_derivative(col, steps[col])[start:end]
                if not np.isfinite(local).all():
                    raise FloatingPointError(f"nonfinite Jacobian column {col}")
                for row, value in zip(range(start, end), local):
                    if value != 0.0:
                        rows.append(row)
                        cols.append(int(col))
                        data.append(float(value))

    jac = coo_matrix((data, (rows, cols)), shape=(n_var, n_var)).tocsc()
    if not np.isfinite(jac.data).all():
        raise FloatingPointError("nonfinite sparse QF Jacobian")
    return jac
