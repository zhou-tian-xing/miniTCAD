#电迁移率
from materials import Material
from dataclasses import dataclass
from typing import final
from math import isfinite
from abc import ABC,abstractmethod
@dataclass(frozen = True,slots = True)
class Mobilities:
    electron_cm2_per_vs: float
    hole_cm2_per_vs: float
class MobilityModel(ABC):
    @final
    def evaluate(
            self,
            material:Material,
            temperature_k:float,
            doping_per_cm3: float | None = None,
            electric_field_v_per_cm: float | None = None
    ) -> Mobilities:
        if ((temperature_k <= 0) or (not isfinite(temperature_k))):
            raise ValueError(
                "temperature_k must be finite and > 0"
            )
        if (doping_per_cm3 is not None) and (
            (not isfinite(doping_per_cm3) or doping_per_cm3 < 0)):
             raise ValueError(
                "doping_per_cm3 must be finite and >= 0"
            )
        if (electric_field_v_per_cm is not None) and not isfinite(electric_field_v_per_cm):
            raise ValueError(
                "electric_field_v_per_cm must be finite"
            )
        result = self._calculate(
            material,
            temperature_k,
            doping_per_cm3,
            electric_field_v_per_cm          
        )
        if not isinstance(result, Mobilities):
            raise TypeError(
                "mobility model must return Mobilities"
            )

        if (
            not isfinite(result.electron_cm2_per_vs)
            or result.electron_cm2_per_vs <= 0
            or not isfinite(result.hole_cm2_per_vs)
            or result.hole_cm2_per_vs <= 0
        ):
            raise ValueError(
                "calculated mobilities must be finite and > 0"
            )

        return result
    @abstractmethod
    def _calculate(
            self,
            material:Material,
            temperature_k:float,
            doping_per_cm3: float | None = None,
            electric_field_v_per_cm: float | None = None
    ) -> Mobilities:
        ...
class ConstantMobilityModel(MobilityModel):
    def _calculate(
        self,
        material:Material,
        temperature_k:float,
        doping_per_cm3: float | None = None,
        electric_field_v_per_cm: float | None = None
    ) -> Mobilities:
        return Mobilities(
            electron_cm2_per_vs = material.electron_mobility_ref_cm2_per_vs,
            hole_cm2_per_vs = material.hole_mobility_ref_cm2_per_vs
        )