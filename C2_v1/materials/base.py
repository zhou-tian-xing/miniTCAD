#给出物理模型的统一定义
from dataclasses import dataclass
from constants import epsilon_0
@dataclass(frozen = True,slots = True)
class Material:
    name: str
    symbol: str
    epsilon_r:float#相对介电常数
    reference_temperature_k: float#各个参考量的参考温度
    intrinsic_concentration_ref_per_cm3: float#本征载流子浓度
    electron_mobility_ref_cm2_per_vs: float#电子迁移率
    hole_mobility_ref_cm2_per_vs: float#空穴迁移率
    @property
    def epsilon_s(self):#绝对介电常数
        return self.epsilon_r * epsilon_0