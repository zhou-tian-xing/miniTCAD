"""一维网格。

Mesh1D 保存结点坐标，以及间距和控制体宽度。
坐标有两套：_cm 是物理单位，_nd 是无量纲（除以 L0）。数组都设成只读，防止被改坏。
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

MIN_NODES = 3


# 冻结array，关闭写权限
def _frozen(values: object, *, name: str) -> NDArray[np.float64]:
    array = np.array(values, dtype=np.float64, copy=True)
    if array.ndim != 1:
        raise ValueError(f"{name} must be a 1D array")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain finite values only")
    array.setflags(write=False)
    return array


@dataclass(frozen=True, slots=True)
class Mesh1D:
    """一维结点网格，结点 i = 0..M-1，坐标从 0 递增到 L。

    x_* 是结点坐标 (M,)，h_* 是相邻间距 (M-1,)，w_* 是内部结点的控制体宽度 (M-2,)。
    """

    x_nodes_cm: NDArray[np.float64]
    h_edges_cm: NDArray[np.float64]
    w_interior_cm: NDArray[np.float64]
    x_nodes_nd: NDArray[np.float64]
    h_edges_nd: NDArray[np.float64]
    w_interior_nd: NDArray[np.float64]

    def __post_init__(self) -> None:
        x_cm = _frozen(self.x_nodes_cm, name="x_nodes_cm")
        h_cm = _frozen(self.h_edges_cm, name="h_edges_cm")
        w_cm = _frozen(self.w_interior_cm, name="w_interior_cm")
        x_nd = _frozen(self.x_nodes_nd, name="x_nodes_nd")
        h_nd = _frozen(self.h_edges_nd, name="h_edges_nd")
        w_nd = _frozen(self.w_interior_nd, name="w_interior_nd")
        object.__setattr__(self, "x_nodes_cm", x_cm)
        object.__setattr__(self, "h_edges_cm", h_cm)
        object.__setattr__(self, "w_interior_cm", w_cm)
        object.__setattr__(self, "x_nodes_nd", x_nd)
        object.__setattr__(self, "h_edges_nd", h_nd)
        object.__setattr__(self, "w_interior_nd", w_nd)

        # 一些检查
        m = x_cm.size
        if m < MIN_NODES:
            raise ValueError(f"a mesh needs at least {MIN_NODES} nodes, got {m}")
        if h_cm.size != m - 1:
            raise ValueError(f"h_edges_cm must have length M-1={m - 1}")
        if w_cm.size != m - 2:
            raise ValueError(f"w_interior_cm must have length M-2={m - 2}")
        if h_nd.size != m - 1 or w_nd.size != m - 2 or x_nd.size != m:
            raise ValueError("cm and nd meshes must have the same shape")
        if x_cm[0] != 0.0:
            raise ValueError(f"x_nodes_cm must start at 0, got {x_cm[0]}")
        if not np.all(np.diff(x_cm) > 0.0):
            raise ValueError("x_nodes_cm must be strictly increasing")
        if not np.allclose(h_cm, np.diff(x_cm), rtol=0.0, atol=0.0):
            raise ValueError("h_edges_cm must equal diff(x_nodes_cm)")
        # 控制体宽度必须是相邻两段间距的一半之和（第 5.1 节），错取相邻边会静默改物理。
        if not np.allclose(w_cm, 0.5 * (h_cm[:-1] + h_cm[1:]), rtol=1e-12, atol=0.0):
            raise ValueError("w_interior_cm must equal adjacent edge-halfwidth sum")

        # cm 与 nd 两套坐标必须只差同一个比例 L0
        scale = h_nd[0] / h_cm[0]
        if not np.isfinite(scale) or scale <= 0.0:
            raise ValueError("normalization scale derived from the mesh must be > 0")
        same_scale = (
            np.allclose(h_nd, h_cm * scale, rtol=1e-12, atol=0.0)
            and np.allclose(w_nd, w_cm * scale, rtol=1e-12, atol=0.0)
            and np.allclose(x_nd, x_cm * scale, rtol=1e-12, atol=0.0)
        )
        if not same_scale:
            raise ValueError("cm and nd meshes are not related by a single scale L0")

    # 结点总数 M
    @property
    def node_count(self) -> int:
        return int(self.x_nodes_cm.size)

    # 物理长度 [cm]
    @property
    def length_cm(self) -> float:
        return float(self.x_nodes_cm[-1])

    # 无量纲长度
    @property
    def length_nd(self) -> float:
        return float(self.x_nodes_nd[-1])

    # 由网格反推 L0 [cm]，用来交叉核对
    def normalization_scale(self) -> float:
        return float(self.h_edges_cm[0] / self.h_edges_nd[0])


# 输入节点位置，自动创建mesh
def _build(x_nodes_cm: NDArray[np.float64], l0_cm: float) -> Mesh1D:
    x = np.asarray(x_nodes_cm, dtype=np.float64)
    h = np.diff(x)
    w = 0.5 * (h[:-1] + h[1:])
    return Mesh1D(
        x_nodes_cm=x,
        h_edges_cm=h,
        w_interior_cm=w,
        x_nodes_nd=x / l0_cm,
        h_edges_nd=h / l0_cm,
        w_interior_nd=w / l0_cm,
    )


def uniform(length_cm: float, nodes: int, l0_cm: float) -> Mesh1D:
    """生成均匀网格：把 [0, length_cm] 等分成 nodes-1 段。

    nodes 是结点总数 M。l0_cm 用来同时算出无量纲网格。
    """
    if not np.isfinite(length_cm) or length_cm <= 0.0:
        raise ValueError("length_cm must be finite and > 0")
    if not np.isfinite(l0_cm) or l0_cm <= 0.0:
        raise ValueError("l0_cm must be finite and > 0")
    if isinstance(nodes, bool) or not isinstance(nodes, (int, np.integer)):
        raise ValueError(f"nodes must be an integer, got {nodes!r}")
    if nodes < MIN_NODES:
        raise ValueError(f"nodes must be >= {MIN_NODES}, got {nodes}")
    return _build(np.linspace(0.0, length_cm, int(nodes)), l0_cm)


def segmented(
    length_cm: float,
    breakpoints_cm: object,
    per_segment_nodes: object,
    l0_cm: float,
) -> Mesh1D:
    """生成分段网格：每段各自等分。

    breakpoints_cm 是分段点，从 0 到 length_cm 递增；per_segment_nodes[k] 是第 k 段的
    结点数（含两端，至少 2）。
    """
    if not np.isfinite(length_cm) or length_cm <= 0.0:
        raise ValueError("length_cm must be finite and > 0")
    if not np.isfinite(l0_cm) or l0_cm <= 0.0:
        raise ValueError("l0_cm must be finite and > 0")

    breaks = np.array(breakpoints_cm, dtype=np.float64)
    raw_counts = np.asarray(per_segment_nodes)
    if raw_counts.ndim != 1 or not np.issubdtype(raw_counts.dtype, np.integer):
        raise ValueError("per_segment_nodes must contain integer node counts")
    counts = raw_counts.astype(np.int64)
    if breaks.ndim != 1 or breaks.size < 2:
        raise ValueError("breakpoints_cm must list at least two coordinates")
    if counts.shape != (breaks.size - 1,):
        raise ValueError("per_segment_nodes must have one entry per segment")
    if not np.all(np.isfinite(breaks)):
        raise ValueError("breakpoints_cm must be finite")
    if breaks[0] != 0.0 or breaks[-1] != length_cm:
        raise ValueError("breakpoints_cm must start at 0 and end at length_cm")
    if not np.all(np.diff(breaks) > 0.0):
        raise ValueError("breakpoints_cm must be strictly increasing")
    if np.any(counts < 2):
        raise ValueError("each segment needs at least 2 nodes")

    # 生成每两个间断点之间的网格，然后拼在一起，注意去除重复点
    pieces = []
    for k in range(breaks.size - 1):
        piece = np.linspace(breaks[k], breaks[k + 1], int(counts[k]))
        pieces.append(piece if k == 0 else piece[1:])
    return _build(np.concatenate(pieces), l0_cm)


def refined(
    length_cm: float,
    junction_cm: float,
    nodes: int,
    l0_cm: float,
    *,
    refine: float = 20.0,
) -> Mesh1D:
    """结区加密网格：结点**总数固定**，结区一侧最密，向两端等比变疏。

    `refine` 是**远端单元与结区单元的长度比**（无量纲，`refine=1` 即均匀）。用比例而不是
    公比来参数化，是因为公比的绝对值依赖段数：`M ≈ 200` 时公比 `1.15` 会给出 1000 倍的
    跨度，而 `1.03` 才对应 20 倍——写成比例就不会选错。

    结区落在结点上（加密的意义所在）。两侧各是一条等比数列：某侧长 `l`、`n` 段、公比 `r`，
    首段（贴结区）为 `h0 = l (r-1) / (r^n - 1)`，`h_k = h0 r^k`，故 `sum h_k = l`；
    `r` 由 `r^(n-1) = refine` 解出。

    段数 `nodes - 1` 按结区位置分配，四舍五入到整数，两侧都至少一段。所以给定 `nodes`
    一定得到 `nodes` 个结点，可与同结点数的 `uniform` 直接比较（§15 实验 C 要求总结点数不变）。

    参数非法一律显式失败：`refine < 1`、结区不在 `(0, length_cm)` 内、结点数太少都报错。
    """
    if not np.isfinite(length_cm) or length_cm <= 0.0:
        raise ValueError("length_cm must be finite and > 0")
    if not np.isfinite(l0_cm) or l0_cm <= 0.0:
        raise ValueError("l0_cm must be finite and > 0")
    if isinstance(nodes, bool) or not isinstance(nodes, (int, np.integer)):
        raise ValueError(f"nodes must be an integer, got {nodes!r}")
    if nodes < MIN_NODES:
        raise ValueError(f"nodes must be >= {MIN_NODES}, got {nodes}")
    if not np.isfinite(refine) or refine < 1.0:
        raise ValueError(f"refine must be finite and >= 1, got {refine}")
    if not np.isfinite(junction_cm) or not (0.0 < junction_cm < length_cm):
        raise ValueError(
            f"junction_cm must lie strictly inside (0, {length_cm}), got {junction_cm}"
        )

    if np.isclose(refine, 1.0, rtol=0.0, atol=1e-12):
        # 不加密就是均匀网格：直接复用 uniform，让两者逐位相同。
        return uniform(length_cm, int(nodes), l0_cm)

    cells = int(nodes) - 1
    left_cells = int(round(cells * junction_cm / length_cm))
    left_cells = max(1, min(left_cells, cells - 1))
    right_cells = cells - left_cells

    def offsets(length: float, count: int) -> NDArray[np.float64]:
        """距结区一端的 `count+1` 个偏移量，间距从这里向远端按公比增长。"""
        if count == 1:
            return np.array([0.0, length])
        growth = refine ** (1.0 / (count - 1))
        first = length * (growth - 1.0) / (growth**count - 1.0)
        if not np.isfinite(first) or first <= 0.0:
            raise ValueError(
                f"refine={refine} over {count} cells cannot resolve the first step; "
                "reduce refine or the node count"
            )
        return np.concatenate(
            ([0.0], np.cumsum(first * growth ** np.arange(count, dtype=np.float64)))
        )

    # 两侧都以“贴结区的一段”为最小段，所以从结区向外生成再映射回去。
    left = junction_cm - offsets(junction_cm, left_cells)
    right = junction_cm + offsets(length_cm - junction_cm, right_cells)
    x = np.concatenate((left[::-1], right[1:]))
    x[0] = 0.0
    x[-1] = length_cm
    return _build(x, l0_cm)
