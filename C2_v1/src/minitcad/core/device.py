"""由配置构造归一化器件。

DeviceConfig 存 YAML/API 里的 cm 字段；NormalizedDevice 是给求解器用的，
存归一化网格、掺杂、迁移率、寿命，以及这次运行的 cm 基准尺度（L0、n_i、D0、J0、R0）。
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite

import numpy as np
from numpy.typing import NDArray

from minitcad.core import mesh as mesh_module
from minitcad.core.contact import contact_targets, validate_pn_orientation
from minitcad.core.material import Material, material_parameter_policy
from minitcad.core.mesh import Mesh1D
from minitcad.core.region import DopingProfile, Region
from minitcad.physics.normalization import NormalizationScales, build_normalization_scales

DEVICE_PN = "pn"
DEVICE_GENERIC_1D = "generic_1d"
ORIENTATION_P_LEFT_N_RIGHT = "p_left_n_right"


@dataclass(frozen=True, slots=True)
class MeshSpec:
    """YAML 里的网格声明。uniform 用 nodes；segmented 用断点和每段结点数。"""

    kind: str
    nodes: int | None = None
    breakpoints_cm: tuple[float, ...] = ()
    per_segment_nodes: tuple[int, ...] = ()

    def __post_init__(self) -> None:
        if self.kind not in ("uniform", "segmented"):
            raise ValueError(f"unknown mesh kind: {self.kind!r}")
        if self.kind == "uniform":
            if (self.nodes is None or isinstance(self.nodes, bool) or
                    not isinstance(self.nodes, (int, np.integer)) or self.nodes < mesh_module.MIN_NODES):
                raise ValueError(f"uniform mesh needs integer nodes >= {mesh_module.MIN_NODES}")
        else:
            if len(self.breakpoints_cm) < 2:
                raise ValueError("segmented mesh needs at least two breakpoints")
            if len(self.per_segment_nodes) != len(self.breakpoints_cm) - 1:
                raise ValueError("per_segment_nodes must have one entry per segment")
            if any(isinstance(v, bool) or not isinstance(v, (int, np.integer)) for v in self.per_segment_nodes):
                raise ValueError("per_segment_nodes must contain integer node counts")

    def build(self, length_cm: float, l0_cm: float) -> Mesh1D:
        """按声明生成网格。"""
        if self.kind == "uniform":
            return mesh_module.uniform(length_cm, int(self.nodes), l0_cm)
        return mesh_module.segmented(length_cm, self.breakpoints_cm, self.per_segment_nodes, l0_cm)


@dataclass(frozen=True, slots=True)
class NormalizationConfig:
    """D0 的两种给法，只能选一个：直接给 D0，或者给参考迁移率 μ0 再算 D0 = μ0 * V_T。"""

    reference_diffusivity_cm2_s: float | None = None
    reference_mobility_cm2_per_Vs: float | None = None

    def __post_init__(self) -> None:
        both_missing = self.reference_diffusivity_cm2_s is None and self.reference_mobility_cm2_per_Vs is None
        both_given = self.reference_diffusivity_cm2_s is not None and self.reference_mobility_cm2_per_Vs is not None
        if both_missing or both_given:
            raise ValueError(
                "give exactly one of reference_diffusivity_cm2_s or reference_mobility_cm2_per_Vs"
            )

    def scales(self, material: Material, temperature_k: float) -> NormalizationScales:
        """算出这次的归一化尺度。"""
        return build_normalization_scales(
            material,
            temperature_k,
            reference_diffusivity_cm2_s=self.reference_diffusivity_cm2_s,
            reference_mobility_cm2_per_Vs=self.reference_mobility_cm2_per_Vs,
        )


@dataclass(frozen=True, slots=True)
class DeviceConfig:
    """一个器件的配置，字段都是 cm 工程单位。"""

    length_cm: float
    area_cm2: float | None
    temperature_k: float
    material: Material
    regions: tuple[Region, ...]
    mesh: MeshSpec
    normalization: NormalizationConfig
    device_type: str = DEVICE_PN
    orientation: str = ORIENTATION_P_LEFT_N_RIGHT
    right_voltage_V: float = 0.0

    def __post_init__(self) -> None:
        if not isfinite(self.length_cm) or self.length_cm <= 0.0:
            raise ValueError("length_cm must be finite and > 0")
        if self.area_cm2 is not None and (not isfinite(self.area_cm2) or self.area_cm2 <= 0.0):
            raise ValueError("area_cm2 must be finite and > 0 when provided")
        if not isfinite(self.temperature_k) or self.temperature_k <= 0.0:
            raise ValueError("temperature_k must be finite and > 0")
        if not isfinite(self.right_voltage_V) or self.right_voltage_V != 0.0:
            raise ValueError("right terminal must be grounded at 0 V")
        if self.device_type not in (DEVICE_PN, DEVICE_GENERIC_1D):
            raise ValueError(f"unknown device type: {self.device_type!r}")
        if self.device_type == DEVICE_GENERIC_1D and self.orientation == ORIENTATION_P_LEFT_N_RIGHT:
            raise ValueError("generic_1d cannot claim p_left_n_right orientation")
        if self.device_type == DEVICE_PN:
            if self.orientation != ORIENTATION_P_LEFT_N_RIGHT:
                raise ValueError(
                    f"device.type pn requires orientation {ORIENTATION_P_LEFT_N_RIGHT!r}, "
                    f"got {self.orientation!r}"
                )
            validate_pn_orientation(self.doping_profile)

    @property
    def doping_profile(self) -> DopingProfile:
        """由 regions 拼出整条掺杂剖面。"""
        return DopingProfile(self.regions, self.length_cm)


@dataclass(frozen=True, slots=True)
class NormalizedDevice:
    """归一化之后的器件，求解器直接用它。"""

    mesh: Mesh1D
    scales: NormalizationScales
    material: Material
    temperature_k: float
    area_cm2: float | None
    external_bias_reference_v: float
    net_doping_nd: NDArray[np.float64]        # (M,) 结点净掺杂除以 n_i
    cbar_nd: NDArray[np.float64]              # (M-2,) 控制体平均净掺杂除以 n_i
    mu_n_nd: float
    mu_p_nd: float
    tau_n_nd: float
    tau_p_nd: float
    material_parameter_policy: str

    def __post_init__(self) -> None:
        for name, count in (("net_doping_nd", self.mesh.node_count),
                            ("cbar_nd", self.mesh.node_count - 2)):
            array = np.asarray(getattr(self, name), dtype=np.float64)
            if array.shape != (count,) or not np.isfinite(array).all():
                raise ValueError(f"{name} must be finite, shape ({count},)")
            array = np.frombuffer(array.tobytes(), dtype=np.float64)
            object.__setattr__(self, name, array)
        for name in ("mu_n_nd", "mu_p_nd", "tau_n_nd", "tau_p_nd"):
            value = getattr(self, name)
            if not isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be finite and > 0")

    @property
    def thermal_voltage_v(self) -> float:
        return self.scales.v0_v

    @property
    def net_doping_left_nd(self) -> float:
        return float(self.net_doping_nd[0])

    @property
    def net_doping_right_nd(self) -> float:
        return float(self.net_doping_nd[-1])

    @property
    def ni_cm3(self) -> float:
        return self.scales.n0_cm3

    @property
    def l0_cm(self) -> float:
        return self.scales.l0_cm

    @property
    def d0_cm2_s(self) -> float:
        return self.scales.d0_cm2_s

    def contact_targets(self, bias_V: float) -> tuple[float, float, float, float, float, float]:
        """见 core.contact.contact_targets。"""
        return contact_targets(self, bias_V)


def build_device(config: DeviceConfig) -> NormalizedDevice:
    """由 DeviceConfig 构造 NormalizedDevice。"""
    policy = material_parameter_policy(config.material, config.temperature_k)
    scales = config.normalization.scales(config.material, config.temperature_k)
    mesh = config.mesh.build(config.length_cm, scales.l0_cm)
    profile = config.doping_profile
    if config.mesh.kind == "segmented":
        for region in profile.regions[:-1]:
            if region.end_cm not in config.mesh.breakpoints_cm:
                raise ValueError(f"segmented mesh must include junction {region.end_cm} cm")

    ni = scales.n0_cm3
    net_doping_nd = profile.net_doping_cm3(mesh.x_nodes_cm) / ni
    cbar_nd = profile.cv_averaged_net_doping_cm3(mesh.x_nodes_cm, mesh.h_edges_cm) / ni

    return NormalizedDevice(
        mesh=mesh,
        scales=scales,
        material=config.material,
        temperature_k=config.temperature_k,
        area_cm2=config.area_cm2,
        external_bias_reference_v=config.right_voltage_V,
        net_doping_nd=net_doping_nd,
        cbar_nd=cbar_nd,
        mu_n_nd=config.material.mobility_n_cm2_per_Vs / scales.mu0_cm2_per_vs,
        mu_p_nd=config.material.mobility_p_cm2_per_Vs / scales.mu0_cm2_per_vs,
        tau_n_nd=config.material.lifetime_n_s / scales.t0_s,
        tau_p_nd=config.material.lifetime_p_s / scales.t0_s,
        material_parameter_policy=policy,
    )