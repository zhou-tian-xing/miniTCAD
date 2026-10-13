"""求解状态和结果。

State 存电势 u 和载流子 n、p（都是无量纲，n、p 必须为正），数组只读。
状态可以打包成一维向量 X = [u0, n0, p0, u1, n1, p1, ...]，长度 3M，供求解器使用。
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray


def _readonly(values: object, name: str) -> NDArray[np.float64]:
    """转成 float64 一维数组、只读，并检查都是有限值。"""
    array = np.array(values, dtype=np.float64, copy=True)
    if array.ndim != 1:
        raise ValueError(f"{name} must be a 1D array")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain finite values only")
    array.setflags(write=False)
    return array


@dataclass(frozen=True, slots=True)
class State:
    """某次求解的状态：电势 u 和载流子 n、p，长度都是 M。"""

    u: NDArray[np.float64]
    n: NDArray[np.float64]
    p: NDArray[np.float64]

    def __post_init__(self) -> None:
        u = _readonly(self.u, "u")
        n = _readonly(self.n, "n")
        p = _readonly(self.p, "p")
        if not (u.size == n.size == p.size):
            raise ValueError("u, n, p must share the same length M")
        if np.any(n <= 0.0) or np.any(p <= 0.0):
            raise ValueError("normalized n and p must be strictly positive")
        object.__setattr__(self, "u", u)
        object.__setattr__(self, "n", n)
        object.__setattr__(self, "p", p)

    # 结点总数 M
    @property
    def node_count(self) -> int:
        return int(self.u.size)

    def pack(self) -> NDArray[np.float64]:
        """打包成一维向量 [u0, n0, p0, u1, n1, p1, ...]，长度 3M。"""
        m = self.u.size
        x = np.empty(3 * m, dtype=np.float64)
        x[0::3] = self.u
        x[1::3] = self.n
        x[2::3] = self.p
        return x

    @staticmethod
    def unpack(x: object, node_count: int) -> "State":
        """pack 的逆运算：把长度 3M 的向量还原成 State。"""
        flat = np.asarray(x, dtype=np.float64)
        if flat.ndim != 1 or flat.size != 3 * node_count:
            raise ValueError(f"x must be 1D with length 3*{node_count}={3 * node_count}")
        return State(u=flat[0::3].copy(), n=flat[1::3].copy(), p=flat[2::3].copy())


@dataclass(frozen=True, slots=True)
class Fluxes:
    """每条边上的电流密度，长度都是 M-1。"""

    jn_edges: NDArray[np.float64]
    jp_edges: NDArray[np.float64]

    def __post_init__(self) -> None:
        jn = _readonly(self.jn_edges, "jn_edges")
        jp = _readonly(self.jp_edges, "jp_edges")
        if jn.size != jp.size:
            raise ValueError("jn_edges and jp_edges must share the same length M-1")
        object.__setattr__(self, "jn_edges", jn)
        object.__setattr__(self, "jp_edges", jp)


@dataclass(frozen=True, slots=True)
class SolveResult:
    """一次求解的结果。converged=False 时 state 可能是 None，sweep 不能拿它当初值。"""

    bias_V: float
    state: State | None
    converged: bool
    iterations: int
    residual_norm: float
    update_norm: float
    reason: str
    iteration_history: tuple[dict, ...] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        if not np.isfinite(self.bias_V):
            raise ValueError("bias_V must be finite")
        if self.iterations < 0:
            raise ValueError("iterations must be >= 0")
        if not self.reason:
            raise ValueError("reason must be a non-empty status code")
