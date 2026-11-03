"""器件侧数据结构：网格、掺杂、材料、接触、状态、归一化器件。

公开 API（W2 冻结口径，W3 起不得无通知改名）：

* 网格：``MIN_NODES``、``Mesh1D``、``uniform``、``segmented``
* 掺杂：``Region``、``DopingProfile``
* 材料：``Material``、``SILICON``、``get_material``、``MaterialTemperatureMismatch``
* 接触：``contact_targets``、``contact_targets_qf``、``contact_potential_nd``、
  ``contact_concentrations_nd``、``validate_pn_orientation``、``INVALID_PN_ORIENTATION``
* 状态：``State``、``QFState``、``Fluxes``、``PhysicsDiagnostics``、``SolveResult``
* 器件：``DeviceConfig``、``NormalizedDevice``、``MeshSpec``、``NormalizationConfig``、``build_device``
"""

from minitcad.core.constants import epsilon_0_f_per_cm, k_b_j_per_k, q
from minitcad.core.mesh import MIN_NODES, Mesh1D, segmented, uniform
from minitcad.core.region import DopingProfile, Region
from minitcad.core.material import (
    POLICY_EXPLICIT_AT_T,
    POLICY_NAMED_MODEL,
    SILICON,
    Material,
    MaterialTemperatureMismatch,
    get_material,
    material_parameter_policy,
)
from minitcad.core.contact import (
    INVALID_PN_ORIENTATION,
    InvalidPnOrientation,
    contact_concentrations_nd,
    contact_potential_nd,
    contact_targets,
    contact_targets_qf,
    validate_pn_orientation,
)
from minitcad.core.state import (
    Fluxes,
    PhysicsDiagnostics,
    QFState,
    SolveResult,
    State,
)
from minitcad.core.device import (
    DEVICE_GENERIC_1D,
    DEVICE_PN,
    ORIENTATION_P_LEFT_N_RIGHT,
    DeviceConfig,
    MeshSpec,
    NormalizationConfig,
    NormalizedDevice,
    build_device,
)

__all__ = [
    "q",
    "k_b_j_per_k",
    "epsilon_0_f_per_cm",
    "MIN_NODES",
    "Mesh1D",
    "uniform",
    "segmented",
    "Region",
    "DopingProfile",
    "Material",
    "SILICON",
    "get_material",
    "material_parameter_policy",
    "MaterialTemperatureMismatch",
    "POLICY_EXPLICIT_AT_T",
    "POLICY_NAMED_MODEL",
    "contact_targets",
    "contact_targets_qf",
    "contact_potential_nd",
    "contact_concentrations_nd",
    "validate_pn_orientation",
    "InvalidPnOrientation",
    "INVALID_PN_ORIENTATION",
    "State",
    "QFState",
    "Fluxes",
    "PhysicsDiagnostics",
    "SolveResult",
    "DeviceConfig",
    "NormalizedDevice",
    "MeshSpec",
    "NormalizationConfig",
    "build_device",
    "DEVICE_PN",
    "DEVICE_GENERIC_1D",
    "ORIENTATION_P_LEFT_N_RIGHT",
]
