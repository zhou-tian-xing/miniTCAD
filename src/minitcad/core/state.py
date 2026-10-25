"""求解状态和结果。

公开状态 ``State`` 存电势 u 和载流子 n、p（都是无量纲，n、p 必须为正），数组只读，
可打包成 ``X = [u0, n0, p0, ...]``（长度 3M）。

A 的 QF 重构（2026-10-21）新增求解器内部的 ``QFState(u, eta_n, eta_p, ref_n, ref_p)``：
``n = exp(u + ref_n + eta_n)``、``p = exp(-u - ref_p - eta_p)``，用于 QF 残差与
更高精度的 SRH/SG；``PhysicsDiagnostics`` 则携带收敛状态上的通量、复合率与最终残差。
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
class QFState:
    """求解器内部的准费米能级状态（A 的 QF 重构后的数据结构）。

    与 ``State`` 的关系（无量纲）::

        n = exp( u + ref_n + eta_n)
        p = exp(-u - ref_p - eta_p)

    ``eta_n/eta_p`` 分别是平移后的电子/空穴准费米变量；``ref_n/ref_p`` 是
    每次求解固定的参考平移量，用来把接触边界归一到 0、减少指数运算的抵消。
    ``eta`` 在接触边界上是 Dirichlet 未知量，这也是 QF 版残差能避开
    载流子浓度相减损失、提高精度的原因。

    数组与 ``State`` 一样只读；节点数 ``M >= 3``。``pack_qf`` 的布局为
    ``[u0, eta_n0, eta_p0, ..., u_{M-1}, eta_n_{M-1}, eta_p_{M-1}]``，长度 3M。
    """

    u: NDArray[np.float64]
    eta_n: NDArray[np.float64]
    eta_p: NDArray[np.float64]
    ref_n: float = 0.0
    ref_p: float = 0.0

    def __post_init__(self) -> None:
        u = _readonly(self.u, "u")
        eta_n = _readonly(self.eta_n, "eta_n")
        eta_p = _readonly(self.eta_p, "eta_p")
        if u.size < 3 or eta_n.size != u.size or eta_p.size != u.size:
            raise ValueError("QFState requires equal lengths M >= 3")
        if not np.isfinite([self.ref_n, self.ref_p]).all():
            raise ValueError("QF reference values must be finite")
        object.__setattr__(self, "u", u)
        object.__setattr__(self, "eta_n", eta_n)
        object.__setattr__(self, "eta_p", eta_p)

    @property
    def node_count(self) -> int:
        return int(self.u.size)

    def to_physical(self) -> State:
        """还原成物理状态 ``State(u, n, p)``；溢出或非正值直接失败。"""
        with np.errstate(over="raise", invalid="raise"):
            n = np.exp(self.u + self.ref_n + self.eta_n)
            p = np.exp(-self.u - self.ref_p - self.eta_p)
        return State(u=self.u, n=n, p=p)

    @classmethod
    def from_physical(
        cls,
        state: State,
        *,
        ref_n: float = 0.0,
        ref_p: float = 0.0,
    ) -> "QFState":
        """从 ``State`` 构造 QF 状态；``ref`` 必须与后续求解使用的参考一致。"""
        return cls(
            u=state.u,
            eta_n=np.log(state.n) - state.u - ref_n,
            eta_p=-np.log(state.p) - state.u - ref_p,
            ref_n=ref_n,
            ref_p=ref_p,
        )

    def pack_qf(self) -> NDArray[np.float64]:
        """打包成 ``[u, eta_n, eta_p]`` 结点交错向量，长度 3M。"""
        z = np.empty(3 * self.node_count, dtype=np.float64)
        z[0::3] = self.u
        z[1::3] = self.eta_n
        z[2::3] = self.eta_p
        return z

    @classmethod
    def unpack_qf(
        cls,
        z: object,
        node_count: int,
        *,
        ref_n: float = 0.0,
        ref_p: float = 0.0,
    ) -> "QFState":
        """``pack_qf`` 的逆运算。"""
        packed = np.asarray(z, dtype=np.float64)
        if packed.ndim != 1 or packed.size != 3 * node_count:
            raise ValueError(f"z must be 1D with length 3*{node_count}={3 * node_count}")
        return cls(
            u=packed[0::3].copy(),
            eta_n=packed[1::3].copy(),
            eta_p=packed[2::3].copy(),
            ref_n=ref_n,
            ref_p=ref_p,
        )


@dataclass(frozen=True, slots=True)
class PhysicsDiagnostics:
    """收敛状态上的派生物理量（全部无量纲），供 writer/验证器消费。

    ``fluxes`` 为 M-1 条边的电流，``recombination_nd`` 为 M 个结点的 SRH 速率，
    ``residual_qf_nd`` 是最终 QF 残差（长度 3M）。QF 迭代变量本身不对外暴露。
    """

    fluxes: Fluxes
    recombination_nd: NDArray[np.float64]
    residual_qf_nd: NDArray[np.float64]

    def __post_init__(self) -> None:
        if not isinstance(self.fluxes, Fluxes):
            raise ValueError("fluxes must be a Fluxes instance")
        rate = _readonly(self.recombination_nd, "recombination_nd")
        residual = _readonly(self.residual_qf_nd, "residual_qf_nd")
        m = rate.size
        if m < 3 or self.fluxes.jn_edges.size != m - 1 or residual.size != 3 * m:
            raise ValueError("invalid physics diagnostics shapes")
        object.__setattr__(self, "recombination_nd", rate)
        object.__setattr__(self, "residual_qf_nd", residual)


@dataclass(frozen=True, slots=True)
class SolveResult:
    """一次求解的结果。converged=False 时 state 可能是 None，sweep 不能拿它当初值。

    ``state`` 始终是物理状态 ``State(u,n,p)``；``physics`` 是 QF 求解器在最终状态
    上复算的派生量（W3 起由 A 的 ``evaluate_qf`` 提供），失败时为 ``None``。
    """

    bias_V: float
    state: State | None
    converged: bool
    iterations: int
    residual_norm: float
    update_norm: float
    reason: str
    iteration_history: tuple[dict, ...] = field(default_factory=tuple)
    physics: PhysicsDiagnostics | None = None

    def __post_init__(self) -> None:
        if not np.isfinite(self.bias_V):
            raise ValueError("bias_V must be finite")
        if self.iterations < 0:
            raise ValueError("iterations must be >= 0")
        if not self.reason:
            raise ValueError("reason must be a non-empty status code")
        if self.physics is not None and not isinstance(self.physics, PhysicsDiagnostics):
            raise ValueError("physics must be PhysicsDiagnostics or None")
        if self.physics is not None:
            if self.state is None or not self.converged:
                raise ValueError("physics is only valid for a converged state")
            if self.physics.residual_qf_nd.size != 3 * self.state.node_count:
                raise ValueError("physics/state node_count mismatch")
