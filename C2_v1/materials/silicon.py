#硅模型
from materials.base import Material
SILICON = Material(
    name="Silicon",
    symbol="Si",
    res=11.7,
    reference_temperature_k=300.0,
    intrinsic_concentration_ref_per_cm3=1.0e10,
    electron_mobility_ref_cm2_per_vs=1400.0,
    hole_mobility_ref_cm2_per_vs=450.0,
)