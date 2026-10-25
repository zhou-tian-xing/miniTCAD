"""W3 QF 残差集成 oracle：``evaluate_qf`` / ``residual_qf`` 的独立交叉核对。

A 的 QF 重构把内部未知量从 ``(u, n, p)`` 换成 ``(u, eta_n, eta_p)``，但物理方程
不应改变。这里用两条独立路径核对：

* 内部行：同一物理状态 ``(u, n, p)`` 分别代入旧 ``residual`` 与 ``residual_qf``，
  Poisson/SG/continuity 的内部行必须在舍入误差内一致；
* 接触行：把 ``eta_n/eta_p`` 设成 ``contact_targets_qf`` 的 Dirichlet 值后，
  QF 边界残差必须严格为 0；
* 精度：``n p - 1`` 在 ``np≈1`` 时发生灾难性抵消，QF 的 ``expm1(eta_n-eta_p)``
  仍能分辨小注入，作为精度提升的独立证据。
"""

from __future__ import annotations

import numpy as np
import pytest

from minitcad.core.contact import contact_targets, contact_targets_qf
from minitcad.core.device import (
    DeviceConfig,
    MeshSpec,
    NormalizationConfig,
    build_device,
)
from minitcad.core.material import SILICON
from minitcad.core.region import Region
from minitcad.core.state import QFState, State
from minitcad.discretization.residual import evaluate_qf, residual, residual_qf
from minitcad.physics.recombination import srh_recombination_nd, srh_recombination_qf

LENGTH_CM = 1.0e-3
JUNCTION_CM = 2.5e-4


def _device():
    return build_device(
        DeviceConfig(
            length_cm=LENGTH_CM,
            area_cm2=1.0e-4,
            temperature_k=300.0,
            material=SILICON,
            regions=(
                Region(0.0, JUNCTION_CM, "silicon", acceptor_cm3=1.0e16, donor_cm3=0.0),
                Region(JUNCTION_CM, LENGTH_CM, "silicon", acceptor_cm3=0.0, donor_cm3=2.0e15),
            ),
            mesh=MeshSpec(kind="uniform", nodes=41),
            normalization=NormalizationConfig(reference_diffusivity_cm2_s=10.0),
        )
    )


def _nonequilibrium_qf(device, bias: float = 0.1) -> QFState:
    """构造一个非平衡 QF 状态：u 平滑，eta 在两端之间线性变化。"""
    m = device.mesh.node_count
    targets = device.contact_targets(bias)
    u = np.linspace(targets[0], targets[3], m)
    eta_n = np.linspace(0.0, -0.35, m)
    eta_p = np.linspace(0.0, 0.22, m)
    return QFState(u=u, eta_n=eta_n, eta_p=eta_p)


def test_internal_rows_match_legacy_physical_residual():
    device = _device()
    qf = _nonequilibrium_qf(device)
    physical = qf.to_physical()
    old = residual(physical, device, 0.1)
    new = evaluate_qf(qf, device, 0.1)
    old_grid = old.reshape((device.mesh.node_count, 3))
    new_grid = new.f.reshape((device.mesh.node_count, 3))
    # Poisson 行：两条路径用同一组 n,p。
    np.testing.assert_allclose(
        new_grid[1:-1, 0], old_grid[1:-1, 0], rtol=1e-12, atol=1e-13
    )
    # continuity 行：物理 SG 与 QF SG 必须等价。
    np.testing.assert_allclose(
        new_grid[1:-1, 1:], old_grid[1:-1, 1:], rtol=1e-11, atol=1e-13
    )


def test_qf_boundary_rows_are_exactly_zero_at_targets():
    device = _device()
    bias = 0.15
    qf = _nonequilibrium_qf(device, bias)
    targets = contact_targets_qf(device, bias, ref_n=qf.ref_n, ref_p=qf.ref_p)
    u = np.array(qf.u, copy=True)
    eta_n = np.array(qf.eta_n, copy=True)
    eta_p = np.array(qf.eta_p, copy=True)
    u[0], u[-1] = targets[0], targets[3]
    eta_n[0], eta_n[-1] = targets[1], targets[4]
    eta_p[0], eta_p[-1] = targets[2], targets[5]
    qf = QFState(u=u, eta_n=eta_n, eta_p=eta_p)
    f = residual_qf(qf, device, bias).reshape(
        (device.mesh.node_count, 3)
    )
    assert f[0, 1] == 0.0 and f[0, 2] == 0.0
    assert f[-1, 1] == 0.0 and f[-1, 2] == 0.0
    assert f[0, 0] == 0.0 and f[-1, 0] == 0.0


def test_residual_qf_matches_evaluation_payload():
    device = _device()
    qf = _nonequilibrium_qf(device)
    evaluation = evaluate_qf(qf, device, 0.1)
    np.testing.assert_array_equal(
        residual_qf(qf, device, 0.1), evaluation.f
    )
    physical = qf.to_physical()
    np.testing.assert_allclose(
        evaluation.recombination_nd,
        srh_recombination_nd(physical.n, physical.p, device.tau_n_nd, device.tau_p_nd),
        rtol=1e-12,
        atol=0.0,
    )
    assert evaluation.fluxes.jn_edges.shape == (device.mesh.node_count - 1,)


def test_qf_srh_avoids_cancellation_near_equilibrium():
    """u=40、eta_n=1e-16：np-1 在双精度下舍入成 0，QF 的 expm1 仍能分辨。"""
    u = np.array([40.0, 40.0, 40.0])
    eta_n = np.array([0.0, 1.0e-16, 3.0e-16])
    eta_p = np.zeros(3)
    qf_rate = srh_recombination_qf(u, eta_n, eta_p, 1.0, 1.0)
    n = np.exp(u + eta_n)
    p = np.exp(-u - eta_p)
    physical_numerator = n * p - 1.0
    assert physical_numerator[1] == 0.0  # 物理路径已丢失 1e-16 的信号
    assert physical_numerator[2] == 0.0
    assert np.all(qf_rate[1:] > 0.0)
    expected = np.expm1(eta_n) / (1.0 * (n + 1.0) + 1.0 * (p + 1.0))
    np.testing.assert_allclose(qf_rate, expected, rtol=1e-14, atol=0.0)
    assert srh_recombination_nd(n[1], p[1], 1.0, 1.0) == 0.0


def test_qf_contact_targets_reconstruct_physical_contact():
    device = _device()
    bias = 0.2
    qf = contact_targets_qf(device, bias)
    physical = contact_targets(device, bias)
    assert np.exp(qf[0] + qf[1]) == pytest.approx(physical[1], rel=1e-14)
    assert np.exp(-qf[0] - qf[2]) == pytest.approx(physical[2], rel=1e-14)
    assert np.exp(qf[3] + qf[4]) == pytest.approx(physical[4], rel=1e-14)
    assert np.exp(-qf[3] - qf[5]) == pytest.approx(physical[5], rel=1e-14)


def test_equilibrium_state_has_zero_qf_current():
    """W3 交接 smoke：A 的平衡标量解在 QF 残差下应给出零电流与 np=1。"""
    from minitcad.simulation import solve_equilibrium_scalar

    device = _device()
    t = device.contact_targets(0.0)
    equilibrium = solve_equilibrium_scalar(
        device.mesh.x_nodes_nd, device.cbar_nd, t[0], t[3]
    )
    assert equilibrium.converged and equilibrium.state is not None
    state = equilibrium.state
    np.testing.assert_allclose(state.n * state.p, 1.0, rtol=0.0, atol=1e-12)
    qf = QFState.from_physical(state)
    evaluation = evaluate_qf(qf, device, 0.0)
    assert np.max(np.abs(evaluation.fluxes.jn_edges)) <= 1e-9
    assert np.max(np.abs(evaluation.fluxes.jp_edges)) <= 1e-9
    assert np.max(np.abs(evaluation.recombination_nd)) <= 1e-9
