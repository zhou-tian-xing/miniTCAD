#归一化因子，手写全部后为适配新接口vibecoding
from __future__ import annotations
from dataclasses import dataclass
from math import isfinite, sqrt
import numpy as np
from minitcad.core.constants import k_b_j_per_k, q
from minitcad.core.material import Material, material_parameter_policy
@dataclass(frozen=True, slots=True)
class NormalizationScales:
    v0_v: float
    n0_cm3: float
    l0_cm: float
    mu0_cm2_per_vs: float
    d0_cm2_s: float
    t0_s: float
    e0_v_cm: float
    j0_a_cm2: float
    r0_cm3_s: float
def build_normalization_scales(
    material: Material,
    temperature_k: float,
    *,
    reference_diffusivity_cm2_s: float | None = None,
    reference_mobility_cm2_per_Vs: float | None = None,
) -> NormalizationScales:
    """在实际 T 和 ni(T) 下计算 L0、D0、t0、E0、J0、R0。

    D0 必须显式选定：配置 D0 或参考 μ0 二选一，不能暗用 300 K 的常量。
    """
    material_parameter_policy(material, temperature_k)
    if (reference_diffusivity_cm2_s is None) == (reference_mobility_cm2_per_Vs is None):
        raise ValueError("give exactly one D0 or reference mobility")
    v0 = k_b_j_per_k * temperature_k / q
    n0 = material.intrinsic_density_cm3
    epsilon = material.epsilon_f_per_cm
    if reference_diffusivity_cm2_s is not None:
        d0 = float(reference_diffusivity_cm2_s)
    else:
        d0 = float(reference_mobility_cm2_per_Vs) * v0
    if not isfinite(d0) or d0 <= 0.0:
        raise ValueError("reference diffusion scale must be finite and > 0")
    l0_squared = epsilon * v0 / (q * n0)
    if not isfinite(l0_squared) or l0_squared <= 0.0:
        raise ValueError("invalid normalization length scale")
    l0 = sqrt(l0_squared)
    scales = NormalizationScales(
        v0_v=v0,
        n0_cm3=n0,
        l0_cm=l0,
        mu0_cm2_per_vs=d0 / v0,
        d0_cm2_s=d0,
        t0_s=l0_squared / d0,
        e0_v_cm=v0 / l0,
        j0_a_cm2=q * d0 * n0 / l0,
        r0_cm3_s=d0 * n0 / l0_squared,
    )
    if not all(isfinite(value) and value > 0.0 for value in vars_as_tuple(scales)):
        raise ValueError("normalization scales exceed finite positive numeric range")
    return scales
def vars_as_tuple(scales: NormalizationScales) -> tuple[float, ...]:
    return (
        scales.v0_v, scales.n0_cm3, scales.l0_cm,
        scales.mu0_cm2_per_vs, scales.d0_cm2_s, scales.t0_s,
        scales.e0_v_cm, scales.j0_a_cm2, scales.r0_cm3_s,
    )
def _factor_cm(quantity: str, scales: NormalizationScales) -> float:
    factors = {
        "x_cm": scales.l0_cm,
        "phi_V": scales.v0_v,
        "density_cm3": scales.n0_cm3,
        "net_doping_cm3": scales.n0_cm3,
        "mobility_cm2_per_Vs": scales.mu0_cm2_per_vs,
        "diffusivity_cm2_s": scales.d0_cm2_s,
        "E_V_cm": scales.e0_v_cm,
        "J_A_cm2": scales.j0_a_cm2,
        "R_cm3_s": scales.r0_cm3_s,
        "tau_s": scales.t0_s,
    }
    try:
        return factors[quantity]
    except KeyError:
        raise ValueError(f"unsupported cm quantity: {quantity!r}") from None
def _checked_values(value: object) -> np.ndarray:
    data = np.asarray(value, dtype=np.float64)
    if not np.all(np.isfinite(data)):
        raise ValueError("physical/normalized input must be finite")
    return data
def normalize_cm(value: object, quantity: str, scales: NormalizationScales):
    """cm → nd；标量返回 float，数组保留形状。"""
    data = _checked_values(value) / _factor_cm(quantity, scales)
    if not np.all(np.isfinite(data)):
        raise ValueError("nonfinite normalized quantity")
    return float(data) if data.ndim == 0 else data
def denormalize_cm(value: object, quantity: str, scales: NormalizationScales):
    """nd → cm；恢复真正有单位的输出时显式调用。"""
    data = _checked_values(value) * _factor_cm(quantity, scales)
    if not np.all(np.isfinite(data)):
        raise ValueError("nonfinite restored physical quantity")
    return float(data) if data.ndim == 0 else data