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
    return np.frombuffer(array.tobytes(), dtype=np.float64)


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
class PhysicsDiagnostics:
    """
    Physical numerical diagnostics at the final accepted state.

    Values are dimensionless.
    No quasi-Fermi iteration variables are exposed.
    """

    fluxes: Fluxes
    recombination_nd: NDArray[np.float64]
    residual_qf_nd: NDArray[np.float64]

    def __post_init__(self) -> None:
        rate = _readonly(
            self.recombination_nd,
            "recombination_nd",
        )
        residual = _readonly(
            self.residual_qf_nd,
            "residual_qf_nd",
        )

        m = rate.size
        if (
            m < 3
            or self.fluxes.jn_edges.size != m - 1
            or residual.size != 3 * m
        ):
            raise ValueError("Invalid physics diagnostics shapes")

        object.__setattr__(self, "recombination_nd", rate)
        object.__setattr__(self, "residual_qf_nd", residual)

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
    physics : PhysicsDiagnostics | None = None

    def __post_init__(self) -> None:
        if not np.isfinite(self.bias_V):
            raise ValueError("bias_V must be finite")
        if self.iterations < 0:
            raise ValueError("iterations must be >= 0")
        if not self.reason:
            raise ValueError("reason must be a non-empty status code")

@dataclass(frozen=True, slots=True)
class QFState:
    """
    Solver-internal quasi-Fermi state.

    u: normalized electrostatic potential
    eta_n: shifted electron quasi-Fermi variable
    eta_p: shifted hole quasi-Fermi variable

    n = exp(u + ref_n + eta_n)
    p = exp(-u - ref_p - eta_p)
    """

    u: NDArray[np.float64]
    eta_n: NDArray[np.float64]
    eta_p: NDArray[np.float64]
    ref_n: float = 0.0
    ref_p: float = 0.0

    def __post_init__(self) -> None:
        u = _readonly(self.u, "u")
        en = _readonly(self.eta_n, "eta_n")
        ep = _readonly(self.eta_p, "eta_p")

        if u.size < 3 or en.size != u.size or ep.size != u.size:
            raise ValueError("QFState requires equal lengths M >= 3")

        if not np.isfinite([self.ref_n, self.ref_p]).all():
            raise ValueError("QF reference values must be finite")

        object.__setattr__(self, "u", u)
        object.__setattr__(self, "eta_n", en)
        object.__setattr__(self, "eta_p", ep)

    @property
    def node_count(self) -> int:
        return int(self.u.size)

    def to_physical(self) -> State:
        """Recover normalized carrier densities."""
        with np.errstate(over="raise", invalid="raise"):
            n = np.exp(self.u + self.ref_n + self.eta_n)
            p = np.exp(-self.u - self.ref_p - self.eta_p)

        # State verifies finite, strictly positive n and p.
        return State(u=self.u, n=n, p=p)

    @classmethod
    def from_physical(
        cls,
        state: State,
        *,
        ref_n: float = 0.0,
        ref_p: float = 0.0,
    ) -> "QFState":
        return cls(
            u=state.u,
            eta_n=np.log(state.n) - state.u - ref_n,
            eta_p=-np.log(state.p) - state.u - ref_p,
            ref_n=ref_n,
            ref_p=ref_p,
        )

    def pack_qf(self) -> NDArray[np.float64]:
        """Interleaved [u0, eta_n0, eta_p0, ...]."""
        z = np.empty(3 * self.node_count, dtype=np.float64)
        z[0::3] = self.u
        z[1::3] = self.eta_n
        z[2::3] = self.eta_p
        return z

    @classmethod
    def unpack_qf(
        cls,
        z: NDArray[np.float64],
        node_count: int,
        *,
        ref_n: float = 0.0,
        ref_p: float = 0.0,
    ) -> "QFState":
        z = np.asarray(z, dtype=np.float64)

        if z.ndim != 1 or z.size != 3 * node_count:
            raise ValueError("Invalid QF vector shape")

        return cls(
            u=z[0::3],
            eta_n=z[1::3],
            eta_p=z[2::3],
            ref_n=ref_n,
            ref_p=ref_p,
        )