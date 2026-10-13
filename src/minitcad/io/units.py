"""cm 和 SI 之间的单位换算。

配置里用 units.system 声明输入是 cm 还是 si。si 输入先在这里换成 cm，
之后走同一条流程。字段名后缀要和声明的单位制对上：
cm 用 _cm/_cm2/_cm3/_cm2_per_Vs/_cm2_s，si 用 _m/_m2/_m3/_m2_per_Vs/_m2_s。
两种后缀不能混用，也不能出现没有单位的物理字段（例如裸 density）。
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite

CM = "cm"
SI = "si"

# 每个物理量从 cm 换成 SI 要乘的系数
CM_TO_SI_FACTOR: dict[str, float] = {
    "length": 1.0e-2,
    "area": 1.0e-4,
    "density": 1.0e6,
    "mobility": 1.0e-4,
    "diffusivity": 1.0e-4,
    "electric_field": 1.0e2,
    "current_density": 1.0e4,
    "recombination_rate": 1.0e6,
    "charge_density": 1.0e6,
    "permittivity": 1.0e2,
    # 下列量在 cm 与 SI 下数值不变（V、A、s、K）。
    "potential": 1.0,
    "voltage": 1.0,
    "current": 1.0,
    "lifetime": 1.0,
    "temperature": 1.0,
}

CM_SUFFIXES = ("_cm", "_cm2", "_cm3", "_cm2_per_Vs", "_cm2_s")
SI_SUFFIXES = ("_m", "_m2", "_m3", "_m2_per_Vs", "_m2_s")


class UnitError(ValueError):
    """单位制声明、字段后缀或换算不符合约定。"""


@dataclass(frozen=True, slots=True)
class UnitSystem:
    """输入单位制，只能是 cm 或 si。"""

    name: str

    def __post_init__(self) -> None:
        if self.name not in (CM, SI):
            raise UnitError(f"units.system must be 'cm' or 'si', got {self.name!r}")

    @classmethod
    def parse(cls, value: object) -> "UnitSystem":
        if not isinstance(value, str):
            raise UnitError(f"units.system must be a string, got {type(value).__name__}")
        return cls(value.strip().lower())


def _factor(quantity: str) -> float:
    try:
        return CM_TO_SI_FACTOR[quantity]
    except KeyError:
        raise UnitError(f"unknown quantity {quantity!r}; known: {sorted(CM_TO_SI_FACTOR)}") from None


def cm_to_si(quantity: str, value: float) -> float:
    """把 cm 单位的值换成 SI。"""
    if not isfinite(value):
        raise UnitError("value must be finite")
    return value * _factor(quantity)


def si_to_cm(quantity: str, value: float) -> float:
    """把 SI 的值换成 cm 单位。"""
    if not isfinite(value):
        raise UnitError("value must be finite")
    return value / _factor(quantity)


def convert(quantity: str, value: float, from_system: str, to_system: str) -> float:
    """在 cm 和 si 之间换算；两边相同时原样返回。"""
    source = UnitSystem.parse(from_system)
    target = UnitSystem.parse(to_system)
    if source.name == target.name:
        return value
    if source.name == CM:
        return cm_to_si(quantity, value)
    return si_to_cm(quantity, value)


def require_suffixed_field(field_name: str, system: str | UnitSystem) -> str:
    """检查字段名后缀和声明的单位制一致，返回命中的后缀。

    三种情况报错：后缀属于另一种单位制、字段没有后缀、字段名不是物理量。
    """
    declared = system.name if isinstance(system, UnitSystem) else UnitSystem.parse(system).name
    wanted = CM_SUFFIXES if declared == CM else SI_SUFFIXES
    other_name = SI if declared == CM else CM
    other = SI_SUFFIXES if declared == CM else CM_SUFFIXES

    for suffix in wanted:
        if field_name.endswith(suffix):
            return suffix
    for suffix in other:
        if field_name.endswith(suffix):
            raise UnitError(
                f"field {field_name!r} uses a {other_name} suffix but units.system is "
                f"{declared!r}; do not mix unit systems in one config"
            )
    raise UnitError(
        f"field {field_name!r} carries no unit suffix; declare units explicitly"
    )
