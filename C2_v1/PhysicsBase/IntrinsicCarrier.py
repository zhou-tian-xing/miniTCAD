#计算本征载流子浓度ni
from materials import Material
from typing import Protocol
class IntrinsicCarrierModel (Protocol):
    def intrinsic_concentration(
            self,
            material:Material,
            temperature_k:float
    ) -> float:
        ...
class ReferenceIntrinsicCarrierModel:
    def intrinsic_concentration(
            self,
            material:Material,
            temperature_k:float
    )-> float:
        if (temperature_k != material.reference_temperature_k):
            raise ValueError(
                "Reference intrinsic concentration is only valid "
                f"at {material.reference_temperature_k} K; "
                f"received {temperature_k} K."
            )
        return material.intrinsic_concentration_ref_per_cm3
class ConstantIntrinsicCarrierModel:
    def intrinsic_concentration(
            self,
            material:Material,
            temperature_k:float
    )-> float:
        return material.intrinsic_concentration_ref_per_cm3