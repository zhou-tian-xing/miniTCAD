"""把求解结果展开成第 12.1 节需要的物理量（cm 工程单位）。

writer 只做“已恢复物理量”的序列化与可选 SI 换算，不重算物理，也不碰无量纲内核。
本模块是那一层“恢复”：无量纲 `State`/`Fluxes` 乘以归一化尺度，得到
`phi_V`、`n_cm3`、`E_V_cm`、`J_A_cm2` 等可直接写进 CSV/PNG 的量。

`F_*_nd` 直接取求解器最终残差（`PhysicsDiagnostics.residual_qf_nd`），这是求解器
自己报告的收敛证据；不要在这里用另一套公式重算，否则 writer 会掩盖求解器错误。
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from minitcad.core.constants import q
from minitcad.core.device import DeviceConfig, NormalizedDevice
from minitcad.core.state import SolveResult
from minitcad.simulation.run_result import Attempt

STATUS_PASS = "PASS"


@dataclass(frozen=True, slots=True)
class NodeFields:
    """结点量，长度都是 `M`。"""

    phi_V: np.ndarray
    n_cm3: np.ndarray
    p_cm3: np.ndarray
    nd_cm3: np.ndarray
    na_cm3: np.ndarray
    net_doping_cm3: np.ndarray
    charge_C_cm3: np.ndarray
    srh_cm3_s: np.ndarray
    f_poisson_nd: np.ndarray
    f_electron_nd: np.ndarray
    f_hole_nd: np.ndarray


@dataclass(frozen=True, slots=True)
class EdgeFields:
    """边量，长度都是 `M-1`；`position_cm` 取边中点。"""

    position_cm: np.ndarray
    E_V_cm: np.ndarray
    jn_A_cm2: np.ndarray
    jp_A_cm2: np.ndarray
    jtotal_A_cm2: np.ndarray


@dataclass(frozen=True, slots=True)
class BiasFields:
    """一个已收敛偏压点的全部可写量。"""

    bias_V: float
    node: NodeFields
    edge: EdgeFields
    iterations: int
    runtime_s: float
    reason: str
    iteration_history: tuple[dict, ...]

    def terminal_current_density_A_cm2(self) -> float:
        """终端总电流密度取最右端点，即 N 侧欧姆接触所在的一条边。"""
        return float(self.edge.jtotal_A_cm2[-1])

    def terminal_current_A(self, area_cm2: float) -> float:
        return float(area_cm2) * self.terminal_current_density_A_cm2()


def node_doping_cm3(config: DeviceConfig, x_cm: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """按结点位置取所在段的 `N_D`、`N_A`；结区阶跃点取右段（与 `DopingProfile` 一致）。"""
    regions = tuple(sorted(config.regions, key=lambda region: region.start_cm))
    ends = np.array([region.end_cm for region in regions], dtype=np.float64)
    index = np.clip(np.searchsorted(ends, x_cm, side="right"), 0, len(regions) - 1)
    donors = np.array([region.donor_cm3 for region in regions], dtype=np.float64)
    acceptors = np.array([region.acceptor_cm3 for region in regions], dtype=np.float64)
    return donors[index], acceptors[index]


def fields_for_attempt(
    config: DeviceConfig, device: NormalizedDevice, attempt: Attempt
) -> BiasFields | None:
    """展开一个尝试点；`FAIL`（无状态）返回 `None`，由 writer 只写 sweep/iteration 行。"""
    result: SolveResult | None = attempt.result
    if result is None or not result.converged or result.state is None:
        return None
    if result.physics is None:
        raise ValueError("converged SolveResult must carry PhysicsDiagnostics")

    state = result.state
    scales = device.scales
    mesh = device.mesh

    phi_V = scales.v0_v * state.u
    n_cm3 = state.n * scales.n0_cm3
    p_cm3 = state.p * scales.n0_cm3
    nd_cm3, na_cm3 = node_doping_cm3(config, mesh.x_nodes_cm)
    net_doping_cm3 = device.net_doping_nd * scales.n0_cm3
    charge_C_cm3 = q * (p_cm3 - n_cm3 + net_doping_cm3)
    srh_cm3_s = result.physics.recombination_nd * scales.r0_cm3_s

    residual = result.physics.residual_qf_nd.reshape(mesh.node_count, 3)
    node = NodeFields(
        phi_V=phi_V,
        n_cm3=n_cm3,
        p_cm3=p_cm3,
        nd_cm3=nd_cm3,
        na_cm3=na_cm3,
        net_doping_cm3=net_doping_cm3,
        charge_C_cm3=charge_C_cm3,
        srh_cm3_s=srh_cm3_s,
        f_poisson_nd=residual[:, 0],
        f_electron_nd=residual[:, 1],
        f_hole_nd=residual[:, 2],
    )

    fluxes = result.physics.fluxes
    du_dx = np.diff(state.u) / mesh.h_edges_cm
    jn = fluxes.jn_edges * scales.j0_a_cm2
    jp = fluxes.jp_edges * scales.j0_a_cm2
    edge = EdgeFields(
        position_cm=0.5 * (mesh.x_nodes_cm[:-1] + mesh.x_nodes_cm[1:]),
        E_V_cm=-scales.v0_v * du_dx,
        jn_A_cm2=jn,
        jp_A_cm2=jp,
        jtotal_A_cm2=jn + jp,
    )

    return BiasFields(
        bias_V=float(attempt.bias_V),
        node=node,
        edge=edge,
        iterations=int(result.iterations),
        runtime_s=float(attempt.runtime_s),
        reason=result.reason,
        iteration_history=tuple(result.iteration_history),
    )


def benchmark_fields(config: DeviceConfig, device: NormalizedDevice) -> dict[str, float]:
    """只依赖器件的基准量：接触浓度/接触势、内建电势、0 V 的质量作用与零电流。

    这些都是**精确**理论值，writer 会把它们填进 `theory_value`（`theory_kind=exact`）；
    有偏压的逐点 `phi/n/p` 没有通用闭式解，不在此列（§12.1 的理论/误差规则）。
    """
    left = device.net_doping_left_nd
    right = device.net_doping_right_nd
    from minitcad.core.contact import contact_concentrations_nd, contact_potential_nd

    n_l, p_l = contact_concentrations_nd(left)
    n_r, p_r = contact_concentrations_nd(right)
    a_l = contact_potential_nd(left)
    a_r = contact_potential_nd(right)
    return {
        "contact_n_left_cm3": n_l * device.scales.n0_cm3,
        "contact_p_left_cm3": p_l * device.scales.n0_cm3,
        "contact_n_right_cm3": n_r * device.scales.n0_cm3,
        "contact_p_right_cm3": p_r * device.scales.n0_cm3,
        "contact_potential_left_V": a_l * device.scales.v0_v,
        "contact_potential_right_V": a_r * device.scales.v0_v,
        "built_in_potential_V": (a_r - a_l) * device.scales.v0_v,
        "equilibrium_np_ratio_nd": 1.0,
        "equilibrium_zero_current_A_cm2": 0.0,
    }
