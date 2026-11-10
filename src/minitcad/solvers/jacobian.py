"""Newton 用的 `(M,3)` 未知量 Jacobian：九色稀疏中心差分（W5-A6 → B 的 T10 校验对象）。

**为什么是差分而不是解析求导**：残差 `F(z)` 里含 `expm1`、Bernoulli 函数与分段选择的
SG 通量，手写解析 Jacobian 的每个分支都要和 `discretization/` 保持逐字同步；一旦残差改了
而 Jacobian 没跟上，Newton 的二次收敛会悄悄退化成线性，错误极难定位。差分版本永远与残差
同源，代价只是若干次残差评估。

**为什么可以“九色”**：列 `j` 是结点 `k = j//3` 的第 `c = j%3` 个分量。控制体离散后，
残差第 `r` 行只依赖结点 `k-1, k, k+1` 的未知量（见 `discretization/residual.py` 的装配），
所以按 `k % 3` 分色后，**同色结点两两相距至少 3 个网格间距**，扰动一列不会污染另一列的
差分区间。于是 9 次整向量求值（3 个分量 × 3 种颜色）就能拿到全部 `3M` 列，而不是 `6M` 次。

**为什么只取局部 3 个结点块**：整向量差分 `fp - fm` 在数值上只对 `k-1..k+1` 这 9 行非零，
其余行是 0/0 的抵消噪声（理论上恰好为 0，浮点上是 `1e-16` 量级的垃圾）。截取局部块既
去掉噪声，也避免把稀疏矩阵填满。

**边界列的退化**：若某个扰动点让 `exp` 溢出，整向量求值失败，该列退回单体中心差分
（`local_derivative`）；仍失败则把步长减半重试，最多 12 次后显式失败——不静默填 0。
"""

from __future__ import annotations

import numpy as np
from scipy.sparse import csc_matrix, coo_matrix

from minitcad.core.device import NormalizedDevice
from minitcad.core.state import QFState
from minitcad.discretization.residual import residual_qf

# 扰动点落在指数溢出区是正常的：这里把它当成“该方向不可用”，由调用方退化处理。
_PERTURBATION_FAILURES = (ValueError, FloatingPointError, OverflowError, ArithmeticError)


def jacobian_qf(
    state: QFState,
    device: NormalizedDevice,
    bias_V: float,
    *,
    rel_step: float = 1e-5,
) -> csc_matrix:
    """返回 `dF/dz`，`z` 是 `QFState.pack_qf()` 的 `3M` 向量，形状 `(3M, 3M)` 稀疏 CSC。

    `rel_step` 是相对步长：`delta_j = rel_step * max(1, |z_j|)`。取 `1e-5` 而不是机器
    epsilon 量级，是因为 `exp` 的舍入让残差的相对噪声约在 `1e-14`；中心差分的截断误差
    为 `O(delta^2)`、舍入误差为 `O(eps/delta)`，两者在 `delta ~ 1e-7` 处平衡。这里取
    `1e-5` 偏向抑制舍入，代价是截断误差略大——Newton 只需要方向足够准，不需要 Jacobian
    的每个元素都准（线搜索会用真实残差复核）。
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
                # 步长已小到加不上去：再减半只会得到同一个数，直接失败。
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
