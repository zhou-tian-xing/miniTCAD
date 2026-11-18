"""T10 的独立 Jacobian oracle：对同一个残差做中心差分。

被微分的是 `discretization.residual.residual_qf`（Newton 实际使用的 `F(z)`），
方向是 `QFState.pack_qf()` 的 `3M` 个分量。本模块不 import `solvers.jacobian`，
否则 T10 就成了“用 A 的实现检验 A 的实现”。

越界时不混用差分：某列 `+h` 失败只算 `-h` 一侧并标 `ONE_SIDED_PLUS`，
两侧都失败标 `FAILED`（默认抛异常）。
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

_TRIAL_FAILURES = (ValueError, FloatingPointError, OverflowError, ArithmeticError)


class FdOracleError(ArithmeticError):
    """两侧试探都失败，无法给出该列导数。"""


@dataclass(frozen=True, slots=True)
class FdJacobian:
    """差分 Jacobian 与逐列差分方式。"""

    matrix: NDArray[np.float64]
    statuses: tuple[str, ...]

    def one_sided_columns(self) -> tuple[int, ...]:
        return tuple(
            index
            for index, status in enumerate(self.statuses)
            if status in (ONE_SIDED_PLUS, ONE_SIDED_MINUS)
        )

    def failed_columns(self) -> tuple[int, ...]:
        return tuple(index for index, status in enumerate(self.statuses) if status == FAILED)


def _trial(
    state: QFState, device: NormalizedDevice, bias_V: float, z: np.ndarray
) -> np.ndarray | None:
    """在打包向量 `z` 上求值；溢出或 `n/p` 非正时返回 `None`。"""
    try:
        trial = QFState.unpack_qf(
            z=z, node_count=state.node_count, ref_n=state.ref_n, ref_p=state.ref_p
        )
        physical = trial.to_physical()
        # `exp` 下溢会给出 0，负值不可能；两种情况都按“该方向不可用”处理。
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

    步长 `h_j = rel_step * max(1, |z_j|)`，与 `jacobian_qf` 用同一尺度规则；
    取 `1e-6` 比实现里的 `1e-5` 小一个量级，让两者的误差来源错开。
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

    result = FdJacobian(matrix=matrix, statuses=tuple(statuses))
    if result.failed_columns() and not allow_failed:
        raise FdOracleError(
            f"both one-sided trials failed for columns {result.failed_columns()}; "
            "cannot produce a finite Jacobian"
        )
    return result


def relative_difference(
    analytic: NDArray[np.float64], reference: NDArray[np.float64]
) -> float:
    """`‖a − r‖∞ / max(1, ‖r‖∞)`，即 T10 的判据。"""
    numerator = float(np.linalg.norm(np.asarray(analytic) - np.asarray(reference), ord=np.inf))
    denominator = max(1.0, float(np.linalg.norm(np.asarray(reference), ord=np.inf)))
    return numerator / denominator
