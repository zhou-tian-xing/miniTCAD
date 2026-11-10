"""T10 的独立 Jacobian oracle：对**同一个残差**做中心差分，与 A 的九色差分互相验证。

被微分的对象是 `discretization.residual.residual_qf`，也就是 Newton 实际使用的那一个
`F(z)`；差分方向是 `z = QFState.pack_qf()` 的 `3M` 个分量。这里刻意不 import
`solvers.jacobian`：T10 要回答的是“`jacobian_qf` 是不是 `F` 的真导数”，而不是
“两段代码写得像不像”。

**为什么用中心差分**：单边差分的截断误差是 `O(h)`，中心差分是 `O(h^2)`。要压到 `1e-4`
的相对精度，单边必须取 `h ~ 1e-4`，那时舍入误差 `eps*|F|/h ~ 1e-12` 还不算大，但截断误差
已经吃掉 `1e-4`；中心差分取 `h ~ 1e-6` 就能同时满足两者。因此默认 `rel_step = 1e-6`，
比 `jacobian_qf` 内部的 `1e-5` 小一个量级——两者的误差来源不同，故意错开。

**越界怎么办**：QF 变量经 `exp` 得到 `n,p`，`n` 或 `p` 溢出成 `inf` 时残差不是有限值。
这种情况**不静默混用两种差分**：某列的 `+h` 侧失败就只算 `-h` 侧并标记 `ONE_SIDED_PLUS`，
两侧都失败就把该列标记 `FAILED` 并抛异常（除非调用方显式允许）。`fs` 里边的
`ONE_SIDED` 计数会出现在测试的诊断输出里：一旦 T10 的误差放大，首先要怀疑的就是
“差分退化了”，而不是“A 的 Jacobian 错了”。
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from minitcad.core.device import NormalizedDevice
from minitcad.core.state import QFState
from minitcad.discretization.residual import residual_qf

CENTRAL = "CENTRAL"
ONE_SIDED_PLUS = "ONE_SIDED_PLUS"
ONE_SIDED_MINUS = "ONE_SIDED_MINUS"
FAILED = "FAILED"

# 试探点溢出属于“该方向不可用”，不是程序错误。
_TRIAL_FAILURES = (ValueError, FloatingPointError, OverflowError, ArithmeticError)


class FdOracleError(ArithmeticError):
    """两侧试探都失败、无法给出该列导数。"""


@dataclass(frozen=True, slots=True)
class FdJacobian:
    """差分 Jacobian 与逐列差分方式。

    `matrix` 是稠密的 `(3M, 3M)`：oracle 只用于验证，节点数在 `10^2` 量级时稠密矩阵
    只有几 MB，换来的是“可以整列做统计”的方便。生产路径不会用到它。
    """

    matrix: NDArray[np.float64]
    statuses: tuple[str, ...]
    steps: NDArray[np.float64]

    def one_sided_columns(self) -> tuple[int, ...]:
        return tuple(
            index
            for index, status in enumerate(self.statuses)
            if status in (ONE_SIDED_PLUS, ONE_SIDED_MINUS)
        )

    def failed_columns(self) -> tuple[int, ...]:
        return tuple(index for index, status in enumerate(self.statuses) if status == FAILED)


def _trial(state: QFState, device: NormalizedDevice, bias_V: float, z: np.ndarray) -> np.ndarray | None:
    """在打包向量 `z` 上求值；失败（溢出/非有限）返回 `None`。"""
    try:
        trial = QFState.unpack_qf(
            z=z, node_count=state.node_count, ref_n=state.ref_n, ref_p=state.ref_p
        )
        # 显式走一次 to_physical：n 或 p 非正、非有限都在这里被拦下，
        # 残差本身的有限性检查是第二道关。
        physical = trial.to_physical()
        if not (np.isfinite(physical.n).all() and np.isfinite(physical.p).all()):
            return None
        if np.any(physical.n <= 0.0) or np.any(physical.p <= 0.0):
            return None
        return residual_qf(state=trial, device=device, bias_V=bias_V)
    except _TRIAL_FAILURES:
        return None


def fd_jacobian(
    state: QFState,
    device: NormalizedDevice,
    bias_V: float,
    *,
    rel_step: float = 1e-6,
    allow_failed: bool = False,
) -> FdJacobian:
    """逐列中心差分，返回 `dF/dz`。

    步长按变量尺度取：`h_j = rel_step * max(1, |z_j|)`，与 `jacobian_qf` 的尺度规则一致
    （尺度规则是“接口”的一部分，不跟着它走就没法区分“步长不同”与“导数不同”）。
    """
    if isinstance(rel_step, (bool, np.bool_)) or not np.isfinite(rel_step) or rel_step <= 0.0:
        raise ValueError("rel_step must be positive and finite")

    z = state.pack_qf()
    n_var = z.size
    f0 = _trial(state, device, bias_V, z)
    if f0 is None:
        raise FdOracleError("baseline residual is not finite at the given state")

    steps = rel_step * np.maximum(1.0, np.abs(z))
    matrix = np.full((n_var, n_var), np.nan, dtype=np.float64)
    statuses: list[str] = []

    for column in range(n_var):
        step = float(steps[column])
        plus_z = z.copy()
        minus_z = z.copy()
        plus_z[column] += step
        minus_z[column] -= step
        plus = _trial(state, device, bias_V, plus_z)
        minus = _trial(state, device, bias_V, minus_z)

        if plus is not None and minus is not None:
            matrix[:, column] = (plus - minus) / (2.0 * step)
            statuses.append(CENTRAL)
        elif plus is not None:
            matrix[:, column] = (plus - f0) / step
            statuses.append(ONE_SIDED_PLUS)
        elif minus is not None:
            matrix[:, column] = (f0 - minus) / step
            statuses.append(ONE_SIDED_MINUS)
        else:
            statuses.append(FAILED)

    result = FdJacobian(matrix=matrix, statuses=tuple(statuses), steps=steps)
    if result.failed_columns() and not allow_failed:
        raise FdOracleError(
            f"both one-sided trials failed for columns {result.failed_columns()}; "
            "cannot produce a finite Jacobian"
        )
    return result


def jacobian_vector_product(
    state: QFState,
    device: NormalizedDevice,
    bias_V: float,
    direction: NDArray[np.float64],
    *,
    rel_step: float = 1e-6,
) -> tuple[NDArray[np.float64], FdJacobian]:
    """`J_FD @ direction`，同时返回完整的 `FdJacobian` 以便检查逐列差分方式。"""
    fd = fd_jacobian(state, device, bias_V, rel_step=rel_step)
    return fd.matrix @ direction, fd


def relative_difference(
    analytic: NDArray[np.float64], reference: NDArray[np.float64]
) -> float:
    """`‖a − r‖∞ / max(1, ‖r‖∞)`——T10 的判据（§14.3）。"""
    numerator = float(np.linalg.norm(np.asarray(analytic) - np.asarray(reference), ord=np.inf))
    denominator = max(1.0, float(np.linalg.norm(np.asarray(reference), ord=np.inf)))
    return numerator / denominator
