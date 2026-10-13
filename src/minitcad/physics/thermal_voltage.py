"""热电势 V_T = k_B * T / q。"""

from __future__ import annotations

from math import isfinite

from minitcad.core.constants import k_b_j_per_k, q


def thermal_voltage(temperature_k: float) -> float:
    """返回热电势 [V]；温度非正或非有限就报错。"""
    if not isfinite(temperature_k) or temperature_k <= 0.0:
        raise ValueError("temperature_k must be finite and > 0")
    return k_b_j_per_k * temperature_k / q
