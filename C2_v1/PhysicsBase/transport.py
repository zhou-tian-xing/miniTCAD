#扩散迁移率由Einstein关系给出
from dataclasses import dataclass
from PhysicsBase.thermal_voltage import thermal_voltage
from PhysicsBase.mobility import Mobilities
@dataclass(frozen = True,slots = True)
class DiffusionCoefficients:
    electron_cm2_per_s: float
    hole_cm2_per_s: float
def einstein_diffusion(
        mobilities : Mobilities,
        Temperature : float
) -> DiffusionCoefficients:
    mu = mobilities.electron_cm2_per_vs
    mp = mobilities.hole_cm2_per_vs
    vt = thermal_voltage(Temperature)
    De = mu * vt
    Dp = mp * vt
    return DiffusionCoefficients(
        electron_cm2_per_s=De,
        hole_cm2_per_s=Dp
    )