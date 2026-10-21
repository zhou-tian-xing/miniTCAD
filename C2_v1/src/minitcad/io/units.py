"""B 的 cm↔SI 边界换算小接口：转换物理量，不触碰 nd State/F(X)。"""
from __future__ import annotations

import numpy as np

# 类型标签不猜数值量级；后续 YAML 适配器只能按已声明的 units.system 调用。
_CM_TO_SI = {
    "length": 1.0e-2,
    "area": 1.0e-4,
    "density": 1.0e6,
    "mobility": 1.0e-4,
    "diffusivity": 1.0e-4,
    "electric_field": 1.0e2,
    "current_density": 1.0e4,
    "recombination": 1.0e6,
    "charge_density": 1.0e6,
    "permittivity": 1.0e2,
    "voltage": 1.0,
    "current": 1.0,
    "lifetime": 1.0,
    "temperature": 1.0,
}


def convert_units(value: object, quantity: str, source: str, target: str):
    """显式标记 cm/si，返回目标单位的标量或同形数组。"""
    if source not in ("cm", "si") or target not in ("cm", "si"):
        raise ValueError("units.system must be cm or si")
    if quantity not in _CM_TO_SI:
        raise ValueError(f"unsupported quantity: {quantity!r}")
    arr = np.asarray(value, dtype=np.float64)
    if not np.isfinite(arr).all():
        raise ValueError("unit conversion input must be finite")
    multiplier = _CM_TO_SI[quantity] if source == "cm" and target == "si" else 1.0
    if source == "si" and target == "cm":
        multiplier = 1.0 / _CM_TO_SI[quantity]
    result = arr * multiplier
    if not np.isfinite(result).all():
        raise ValueError("converted quantity is not finite")
    return float(result) if result.ndim == 0 else result
