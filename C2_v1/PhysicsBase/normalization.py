#归一化因子
from dataclasses import dataclass
from math import isfinite, sqrt
from constants import q
from materials import Material
from PhysicsBase.IntrinsicCarrier import IntrinsicCarrierModel
from PhysicsBase.thermal_voltage import thermal_voltage
@dataclass(frozen=True, slots=True)
class NormalizationScales:
    v0_v: float
    n0_per_cm3: float
    l0_cm: float
    mu0_cm2_per_vs: float
    d0_cm2_per_s: float
    t0_s: float
    e0_v_per_cm: float
    j0_a_per_cm2: float
    r0_per_cm3_s: float
    qr0_a_per_cm3: float
def build_normalization_scales(
    material: Material,
    temperature_k: float,
    intrinsic_model: IntrinsicCarrierModel,
    mu0_cm2_per_vs: float = 1000.0,
) -> NormalizationScales:
    if not isfinite(mu0_cm2_per_vs) or mu0_cm2_per_vs <= 0.0:
        raise ValueError(
            "mu0_cm2_per_vs must be positive and finite"
        )
    v0 = thermal_voltage(temperature_k)
    n0 = intrinsic_model.intrinsic_concentration(
        material,
        temperature_k,
    )
    epsilon = material.epsilon_s
    if not isfinite(n0) or n0 <= 0.0:
        raise ValueError(
            "intrinsic concentration must be positive and finite"
        )

    if not isfinite(epsilon) or epsilon <= 0.0:
        raise ValueError(
            "material permittivity must be positive and finite"
        )
    charge_density = q * n0
    if charge_density == 0.0:
        raise ValueError(
            "intrinsic concentration is too small for this numeric range"
        )
    l0_squared = epsilon * v0 / charge_density
    d0 = mu0_cm2_per_vs * v0
    if not isfinite(l0_squared) or l0_squared <= 0.0:
        raise ValueError("invalid length scale")
    if not isfinite(d0) or d0 <= 0.0:
        raise ValueError("invalid diffusion scale")
    l0 = sqrt(l0_squared)
    t0 = l0_squared / d0
    e0 = v0 / l0
    j0 = q * d0 * n0 / l0
    r0 = d0 * n0 / l0_squared
    qr0 = q * r0
    if not all(
        isfinite(s) and s > 0.0
        for s in (l0, t0, e0, j0, r0, qr0)
    ):
        raise ValueError(
            "normalization scales exceed the supported numeric range"
        )
    return NormalizationScales(
        v0_v=v0,
        n0_per_cm3=n0,
        l0_cm=l0,
        mu0_cm2_per_vs=mu0_cm2_per_vs,
        d0_cm2_per_s=d0,
        t0_s=t0,
        e0_v_per_cm=e0,
        j0_a_per_cm2=j0,
        r0_per_cm3_s=r0,
        qr0_a_per_cm3=qr0,
    )