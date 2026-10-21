"""分段掺杂：把器件切成几段，每段掺杂恒定，并算出控制体平均值。

掺杂用净施主浓度 C = N_D - N_A 表示 [cm^-3]。C 为负是 P 型，为正时 N 型。
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite

import numpy as np
from numpy.typing import NDArray


@dataclass(frozen=True, slots=True)
class Region:
    """一段掺杂恒定的区域，用 [start_cm, end_cm] 表示。"""

    start_cm: float
    end_cm: float
    material: str
    acceptor_cm3: float
    donor_cm3: float

    def __post_init__(self) -> None:
        if not self.material:
            raise ValueError("region material must be non-empty")
        if not (isfinite(self.start_cm) and isfinite(self.end_cm)):
            raise ValueError("region bounds must be finite")
        if self.end_cm <= self.start_cm:
            raise ValueError(
                f"region must have end_cm > start_cm, got [{self.start_cm}, {self.end_cm}]"
            )
        for name, value in (("acceptor_cm3", self.acceptor_cm3), ("donor_cm3", self.donor_cm3)):
            if not isfinite(value) or value < 0.0:
                raise ValueError(f"{name} must be finite and >= 0")

    @property
    def net_doping_cm3(self) -> float:
        """净掺杂 C = N_D - N_A [cm^-3]。"""
        return self.donor_cm3 - self.acceptor_cm3


class DopingProfile:
    """由若干首尾相接的 Region 拼成的整条掺杂剖面。

    构造时就检查：必须从 0 开始、到 length_cm 结束、段与段之间既不重叠也不留缝。
    """

    def __init__(self, regions: object, length_cm: float) -> None:
        if not isfinite(length_cm) or length_cm <= 0.0:
            raise ValueError("length_cm must be finite and > 0")
        items = tuple(regions)
        if not items:
            raise ValueError("at least one region is required")

        ordered = sorted(items, key=lambda region: region.start_cm)
        if ordered[0].start_cm != 0.0:
            raise ValueError(
                f"regions must start at 0, first region starts at {ordered[0].start_cm}"
            )
        if ordered[-1].end_cm != length_cm:
            raise ValueError(
                f"regions must end at {length_cm}, last region ends at {ordered[-1].end_cm}"
            )
        for left, right in zip(ordered, ordered[1:]):
            if right.start_cm < left.end_cm:
                raise ValueError(
                    f"regions overlap at {right.start_cm}: previous region ends at {left.end_cm}"
                )
            if right.start_cm > left.end_cm:
                raise ValueError(
                    f"gap between regions from {left.end_cm} to {right.start_cm}"
                )

        self._regions = ordered
        self._length_cm = float(length_cm)
        # 每个结点 x 落在哪一段，用段右端点做二分查找即可。
        self._ends_cm = np.array([region.end_cm for region in ordered], dtype=np.float64)
        self._values_cm3 = np.array([region.net_doping_cm3 for region in ordered], dtype=np.float64)

    @property
    def regions(self) -> tuple[Region, ...]:
        return self._regions

    @property
    def length_cm(self) -> float:
        return self._length_cm

    def net_doping_cm3(self, x_cm: object) -> NDArray[np.float64]:
        """给出 x 处的净掺杂 [cm^-3]，x 可以是标量或数组。

        x 正好落在两段交界处时取右边那一段的值。
        """
        x = np.asarray(x_cm, dtype=np.float64)
        if not np.all(np.isfinite(x)):
            raise ValueError("x_cm must be finite")
        index = np.searchsorted(self._ends_cm, x, side="right")
        index = np.clip(index, 0, self._values_cm3.size - 1)
        return self._values_cm3[index]

    def _integral_cm3(self, a: float, b: float) -> float:
        """净掺杂在 [a, b] 上的积分 [cm^-2]，按段累加。"""
        total = 0.0
        for region in self._regions:
            lo = max(a, region.start_cm)
            hi = min(b, region.end_cm)
            if hi > lo:
                total += region.net_doping_cm3 * (hi - lo)
        return total

    def cv_averaged_net_doping_cm3(
        self, x_nodes_cm: object, h_edges_cm: object
    ) -> NDArray[np.float64]:
        """每个内部结点的控制体平均净掺杂，长度 M-2。

        结点 i 的控制体是 [x_i - h_{i-1}/2, x_i + h_i/2]，平均值等于该区间上的
        积分除以宽度。这样跨在结区上的控制体会自动得到两段的加权平均，
        而不是错拿交界处某一侧的值。
        """
        x = np.asarray(x_nodes_cm, dtype=np.float64)
        h = np.asarray(h_edges_cm, dtype=np.float64)
        if x.ndim != 1 or h.ndim != 1 or h.size != x.size - 1:
            raise ValueError("expect 1D x_nodes_cm and h_edges_cm with len(h) == len(x) - 1")

        out = np.empty(x.size - 2, dtype=np.float64)
        for i in range(1, x.size - 1):
            left = x[i] - 0.5 * h[i - 1]
            right = x[i] + 0.5 * h[i]
            out[i - 1] = self._integral_cm3(left, right) / (right - left)
        return out

    @property
    def net_doping_left_cm3(self) -> float:
        """左端（x=0）所在段的净掺杂。"""
        return self._regions[0].net_doping_cm3

    @property
    def net_doping_right_cm3(self) -> float:
        """右端（x=L）所在段的净掺杂。"""
        return self._regions[-1].net_doping_cm3

    @property
    def is_left_p_right_n(self) -> bool:
        """是否左端 P 型、右端 N 型。"""
        return self.net_doping_left_cm3 < 0.0 < self.net_doping_right_cm3
