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
    if int(nodes) < MIN_NODES:
        raise ValueError(f"nodes must be >= {MIN_NODES}, got {nodes}")
    return _build(np.linspace(0.0, length_cm, int(nodes)), l0_cm)


def segmented(
    length_cm: float,
    breakpoints_cm: object,
    per_segment_nodes: object,
    l0_cm: float,
) -> Mesh1D:
    """生成分段网格：每段各自等分，方便在结区加密。

    breakpoints_cm 是分段点，从 0 到 length_cm 递增。
    per_segment_nodes[k] 是第 k 段的结点数（含两端，至少 2）。
    """
    if not np.isfinite(length_cm) or length_cm <= 0.0:
        raise ValueError("length_cm must be finite and > 0")
    if not np.isfinite(l0_cm) or l0_cm <= 0.0:
        raise ValueError("l0_cm must be finite and > 0")

    breaks = np.array(breakpoints_cm, dtype=np.float64)
    counts = np.array(per_segment_nodes, dtype=np.int64)
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
