"""物理量：热电压、归一化尺度与 SRH 复合。

``thermal_voltage``/``normalization`` 由 B 从 ``C2_v1/PhysicsBase`` 搬运并冻结接口；
``recombination`` 是 A 的 W3 交付（中能级 SRH，含 QF 版）。
"""

from minitcad.physics.normalization import (
    NormalizationScales,
    build_normalization_scales,
)
from minitcad.physics.recombination import (
    srh_recombination_nd,
    srh_recombination_qf,
)

__all__ = [
    "NormalizationScales",
    "build_normalization_scales",
    "srh_recombination_nd",
    "srh_recombination_qf",
]
