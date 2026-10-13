"""材料参数。

每种材料给一组参数，都是 cm 工程单位（能量用 eV）。
绝对介电常数 epsilon_f_per_cm = epsilon_relative * epsilon_0，epsilon_0 取 F/cm。
参数标定的温度和运行时温度不一致时直接报错，不做外推。
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isclose, isfinite

from minitcad.core.constants import epsilon_0_f_per_cm

POLICY_EXPLICIT_AT_T = "explicit_at_T"
POLICY_NAMED_MODEL = "named_model"


class MaterialTemperatureMismatch(ValueError):
    """材料参数标定的温度和运行温度对不上。"""


@dataclass(frozen=True, slots=True)
class Material:
    """一种材料在某个温度下的一组参数。"""

    name: str
    symbol: str
    epsilon_relative: float
    reference_temperature_k: float
    intrinsic_density_cm3: float
    mobility_n_cm2_per_Vs: float
    mobility_p_cm2_per_Vs: float
    lifetime_n_s: float
    lifetime_p_s: float
    bandgap_eV: float

    def __post_init__(self) -> None:
        if not self.name or not self.symbol:
            raise ValueError("name and symbol must be non-empty")
        positive = (
            ("epsilon_relative", self.epsilon_relative),
            ("reference_temperature_k", self.reference_temperature_k),
            ("intrinsic_density_cm3", self.intrinsic_density_cm3),
            ("mobility_n_cm2_per_Vs", self.mobility_n_cm2_per_Vs),
            ("mobility_p_cm2_per_Vs", self.mobility_p_cm2_per_Vs),
            ("lifetime_n_s", self.lifetime_n_s),
            ("lifetime_p_s", self.lifetime_p_s),
            ("bandgap_eV", self.bandgap_eV),
        )
        for field_name, value in positive:
            if not isfinite(value) or value <= 0.0:
                raise ValueError(f"{field_name} must be finite and > 0")

    # 绝对介电常数 [F/cm]
    @property
    def epsilon_f_per_cm(self) -> float:
        return self.epsilon_relative * epsilon_0_f_per_cm


def material_parameter_policy(
    material: Material,
    temperature_k: float,
    *,
    tolerance_k: float = 1.0e-9,
) -> str:
    """检查材料温度和运行温度是否一致。

    一致就返回参数来源标记；不一致就报 MaterialTemperatureMismatch。
    """
    if not isfinite(temperature_k) or temperature_k <= 0.0:
        raise ValueError("temperature_k must be finite and > 0")
    if not isclose(temperature_k, material.reference_temperature_k, rel_tol=0.0, abs_tol=tolerance_k):
        raise MaterialTemperatureMismatch(
            f"{material.name}: parameters are calibrated at "
            f"{material.reference_temperature_k} K but this run uses {temperature_k} K; "
            "supply parameters at the run temperature instead of extrapolating"
        )
    return POLICY_EXPLICIT_AT_T


# 参考材料：数值取自 C2_v1/materials/silicon.py（已实测通过），寿命与带隙按示例配置补齐。
SILICON = Material(
    name="silicon",
    symbol="Si",
    epsilon_relative=11.7,
    reference_temperature_k=300.0,
    intrinsic_density_cm3=1.0e10,
    mobility_n_cm2_per_Vs=1400.0,
    mobility_p_cm2_per_Vs=450.0,
    lifetime_n_s=1.0e-6,
    lifetime_p_s=2.0e-6,
    bandgap_eV=1.12,
)

_REGISTRY: dict[str, Material] = {
    "si": SILICON,
    "silicon": SILICON,
}


def get_material(name: str) -> Material:
    """按名字查参考材料；名字不认识就报错。"""
    key = name.strip().lower()
    try:
        return _REGISTRY[key]
    except KeyError:
        raise ValueError(f"unknown material: {name!r}; known: {sorted(_REGISTRY)}") from None
