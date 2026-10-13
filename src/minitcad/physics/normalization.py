"""归一化尺度：把 cm 单位换算成无量纲时用到的一组基准量。

算出来的 L0、V_T、n_i、D0 等就是各个物理量除以的基准。
D0 由调用方二选一给出：直接给 D0，或者给参考迁移率 mu0 再算 D0 = mu0 * V_T。
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite, sqrt

from minitcad.core.constants import q
from minitcad.core.material import Material
from minitcad.physics.thermal_voltage import thermal_voltage


@dataclass(frozen=True, slots=True)
class NormalizationScales:
    """一次运行固定不变的一组 cm 基准尺度。"""

    v0_v: float
    n0_cm3: float
    l0_cm: float
    mu0_cm2_per_vs: float
    d0_cm2_s: float
    t0_s: float
    e0_v_per_cm: float
    j0_a_cm2: float
    r0_cm3_s: float
    qr0_a_cm3: float


def build_normalization_scales(
    material: Material,
    temperature_k: float,
    *,
    reference_diffusivity_cm2_s: float | None = None,
    reference_mobility_cm2_per_Vs: float | None = None,
) -> NormalizationScales:
    """算出归一化尺度。D0 必须由两种给法之一明确给出。"""
    if (reference_diffusivity_cm2_s is None) == (reference_mobility_cm2_per_Vs is None):
        raise ValueError(
            "give exactly one of reference_diffusivity_cm2_s or reference_mobility_cm2_per_Vs"
        )

    v0 = thermal_voltage(temperature_k)
    n0 = material.intrinsic_density_cm3
    epsilon = material.epsilon_f_per_cm

    if reference_diffusivity_cm2_s is not None:
        d0 = float(reference_diffusivity_cm2_s)
    else:
        d0 = float(reference_mobility_cm2_per_Vs) * v0

    for name, value in (("n0", n0), ("epsilon", epsilon), ("d0", d0)):
        if not isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be finite and > 0")

    mu0 = d0 / v0
    l0_squared = epsilon * v0 / (q * n0)
    if not isfinite(l0_squared) or l0_squared <= 0.0:
        raise ValueError("invalid length scale")

    l0 = sqrt(l0_squared)
    t0 = l0_squared / d0
    e0 = v0 / l0
    j0 = q * d0 * n0 / l0
    r0 = d0 * n0 / l0_squared
    qr0 = q * r0

    if not all(isfinite(s) and s > 0.0 for s in (l0, t0, e0, j0, r0, qr0, mu0)):
        raise ValueError("normalization scales exceed the supported numeric range")

    return NormalizationScales(
        v0_v=v0,
        n0_cm3=n0,
        l0_cm=l0,
        mu0_cm2_per_vs=mu0,
        d0_cm2_s=d0,
        t0_s=t0,
        e0_v_per_cm=e0,
        j0_a_cm2=j0,
        r0_cm3_s=r0,
        qr0_a_cm3=qr0,
    )
