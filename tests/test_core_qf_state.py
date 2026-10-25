"""W2/W3 QF 数据结构：``QFState``、``PhysicsDiagnostics`` 与 QF 接触目标。

A 在 2026-10-21 把求解器核心从 ``(u, n, p)`` 改成准费米变量 ``(u, eta_n, eta_p)``
以提高精度（见 ``C2_v1/src/minitcad/discretization/residual.py`` 头注释）。
``State`` 仍是公开的物理状态；``QFState`` 只是求解器内部表示，但它的只读、
3M 交错排序和接触边界必须与 W2 冻结口径一致（I-01 的 QF 扩展）。
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from minitcad.core.contact import contact_targets, contact_targets_qf
from minitcad.core.state import (
    Fluxes,
    PhysicsDiagnostics,
    QFState,
    SolveResult,
    State,
)


class _DeviceStub:
    """只提供接触目标所需的三个属性，避免依赖完整 DeviceConfig 构建。"""

    def __init__(self, c_left_nd: float, c_right_nd: float, vt_v: float) -> None:
        self.net_doping_left_nd = c_left_nd
        self.net_doping_right_nd = c_right_nd
        self.thermal_voltage_v = vt_v


def _state(m: int = 5) -> State:
    u = np.linspace(-1.0, 1.0, m)
    return State(u=u, n=np.exp(u), p=np.exp(-u))


def _qf(m: int = 5) -> QFState:
    return QFState.from_physical(_state(m), ref_n=0.3, ref_p=-0.2)


def test_qf_round_trip_recovers_physical_state():
    state = _state(7)
    qf = QFState.from_physical(state, ref_n=0.3, ref_p=-0.2)
    restored = qf.to_physical()
    np.testing.assert_array_equal(restored.u, state.u)
    np.testing.assert_allclose(restored.n, state.n, rtol=1e-14, atol=0.0)
    np.testing.assert_allclose(restored.p, state.p, rtol=1e-14, atol=0.0)
    assert qf.ref_n == 0.3 and qf.ref_p == -0.2


def test_pack_qf_is_interleaved_and_round_trips():
    qf = _qf(3)
    z = qf.pack_qf()
    assert z.shape == (9,)
    assert z[0] == qf.u[0] and z[1] == qf.eta_n[0] and z[2] == qf.eta_p[0]
    assert z[3] == qf.u[1] and z[4] == qf.eta_n[1] and z[5] == qf.eta_p[1]
    restored = QFState.unpack_qf(z, 3, ref_n=qf.ref_n, ref_p=qf.ref_p)
    np.testing.assert_array_equal(restored.pack_qf(), z)
    assert restored.ref_n == qf.ref_n and restored.ref_p == qf.ref_p


def test_qf_arrays_are_read_only_and_defensively_copied():
    u = np.zeros(3)
    eta_n = np.zeros(3)
    eta_p = np.zeros(3)
    qf = QFState(u=u, eta_n=eta_n, eta_p=eta_p)
    u[0] = 9.0
    assert qf.u[0] == 0.0
    with pytest.raises(ValueError):
        qf.eta_n[0] = 1.0
    with pytest.raises(ValueError):
        qf.eta_p[0] = 1.0


def test_qf_rejects_bad_shapes_and_nonfinite_references():
    with pytest.raises(ValueError):
        QFState(u=np.zeros(2), eta_n=np.zeros(2), eta_p=np.zeros(2))
    with pytest.raises(ValueError):
        QFState(u=np.zeros(3), eta_n=np.zeros(4), eta_p=np.zeros(3))
    with pytest.raises(ValueError):
        QFState(u=np.zeros(3), eta_n=np.zeros(3), eta_p=np.zeros(3), ref_n=math.inf)
    with pytest.raises(ValueError):
        QFState.unpack_qf(np.zeros(8), 3)


def test_qf_to_physical_rejects_overflow():
    qf = QFState(
        u=np.full(3, 1000.0),
        eta_n=np.zeros(3),
        eta_p=np.zeros(3),
        ref_n=1000.0,
    )
    with pytest.raises(FloatingPointError):
        qf.to_physical()


@pytest.mark.parametrize("bias", [-0.1, 0.0, 0.05, 0.3])
def test_contact_targets_qf_reconstructs_physical_targets(bias):
    device = _DeviceStub(c_left_nd=-1.0e6, c_right_nd=2.0e5, vt_v=0.025852)
    v = bias / device.thermal_voltage_v
    qf = contact_targets_qf(device, bias)
    physical = contact_targets(device, bias)
    u_l, eta_nl, eta_pl, u_r, eta_nr, eta_pr = qf
    assert u_l == pytest.approx(physical[0], abs=0.0)
    assert u_r == pytest.approx(physical[3], abs=0.0)
    # 局部平衡边界：eta_n = eta_p = -v（左）、0（右）。
    assert eta_nl == pytest.approx(-v, abs=1e-14)
    assert eta_pl == pytest.approx(-v, abs=1e-14)
    assert eta_nr == pytest.approx(0.0, abs=1e-14)
    assert eta_pr == pytest.approx(0.0, abs=1e-14)
    reconstructed = (
        math.exp(u_l + eta_nl),
        math.exp(-u_l - eta_pl),
        math.exp(u_r + eta_nr),
        math.exp(-u_r - eta_pr),
    )
    np.testing.assert_allclose(
        reconstructed, [physical[1], physical[2], physical[4], physical[5]],
        rtol=1e-14, atol=0.0,
    )


def test_contact_targets_qf_rejects_nonfinite_reference():
    device = _DeviceStub(c_left_nd=-1.0, c_right_nd=1.0, vt_v=0.025852)
    with pytest.raises(ValueError):
        contact_targets_qf(device, 0.0, ref_n=math.nan)
    with pytest.raises(ValueError):
        contact_targets_qf(device, math.inf)


def _diagnostics(m: int = 4) -> PhysicsDiagnostics:
    return PhysicsDiagnostics(
        fluxes=Fluxes(jn_edges=np.zeros(m - 1), jp_edges=np.zeros(m - 1)),
        recombination_nd=np.zeros(m),
        residual_qf_nd=np.zeros(3 * m),
    )


def test_physics_diagnostics_shapes_and_read_only():
    diagnostics = _diagnostics(4)
    assert diagnostics.recombination_nd.size == 4
    assert diagnostics.residual_qf_nd.size == 12
    with pytest.raises(ValueError):
        diagnostics.recombination_nd[0] = 1.0
    with pytest.raises(ValueError):
        PhysicsDiagnostics(
            fluxes=Fluxes(jn_edges=np.zeros(4), jp_edges=np.zeros(4)),
            recombination_nd=np.zeros(3),
            residual_qf_nd=np.zeros(9),
        )
    with pytest.raises(ValueError):
        PhysicsDiagnostics(
            fluxes=Fluxes(jn_edges=np.zeros(2), jp_edges=np.zeros(2)),
            recombination_nd=np.zeros(4),
            residual_qf_nd=np.zeros(12),
        )


def test_solve_result_physics_requires_converged_state():
    state = _state(4)
    diagnostics = _diagnostics(4)
    result = SolveResult(0.0, state, True, 3, 1e-9, 1e-7, "SUCCESS",
                         physics=diagnostics)
    assert result.physics is diagnostics
    with pytest.raises(ValueError):
        SolveResult(0.0, None, False, 3, 1e-9, 1e-7, "MAX_ITER", physics=diagnostics)
    with pytest.raises(ValueError):
        SolveResult(0.0, state, False, 3, 1e-9, 1e-7, "MAX_ITER", physics=diagnostics)
    with pytest.raises(ValueError):
        SolveResult(0.0, _state(5), True, 3, 1e-9, 1e-7, "SUCCESS",
                    physics=diagnostics)
