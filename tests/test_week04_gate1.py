"""W4-B1：T08（0 V 完整 PN）与 T09（总电流沿空间守恒）（第 11.2 节）。

两个测试都只消费 `api.run()` 的 `RunResult`，即与 GUI/CLI 完全相同的产物，
不直接调用求解器内部函数；期望值来自解析恒等式（`np = 1`、零电流、接触目标），
不复用 A 的实现产生期望值。
"""

from __future__ import annotations

import pathlib

import numpy as np
import pytest

from minitcad.api import run
from minitcad.io.config import load_config

BASELINE = pathlib.Path(__file__).resolve().parents[1] / "examples" / "pn_baseline_0_to_030.yaml"

NP_TOL = 1e-6
ZERO_CURRENT_TOL = 1e-6
CONTACT_TOL = 1e-10
CONSERVATION_TOL = 1e-3


@pytest.fixture(scope="module")
def baseline_run():
    return run(load_config(BASELINE))


def test_run_accepts_every_baseline_bias(baseline_run):
    assert baseline_run.accepted, baseline_run.failure_reason
    assert len(baseline_run.attempts) == 7
    assert len(baseline_run.successful) == 7
    assert [attempt.bias_V for attempt in baseline_run.attempts] == pytest.approx(
        [0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30]
    )


def test_t08_equilibrium_mass_action(baseline_run):
    """0 V：全器件 `np = n_i^2`，偏差 <= 1e-6（§11.2 T08）。"""
    equilibrium = baseline_run.successful[0]
    assert equilibrium.bias_V == 0.0
    state = equilibrium.result.state
    np_ratio = state.n * state.p
    assert np.max(np.abs(np_ratio - 1.0)) <= NP_TOL


def test_t08_equilibrium_has_zero_edge_current(baseline_run):
    equilibrium = baseline_run.successful[0]
    fluxes = equilibrium.result.physics.fluxes
    assert np.max(np.abs(fluxes.jn_edges)) <= ZERO_CURRENT_TOL
    assert np.max(np.abs(fluxes.jp_edges)) <= ZERO_CURRENT_TOL


def test_t08_boundary_neutrality_matches_contact_targets(baseline_run):
    """两端 `u, n, p` 必须落在 `contact_targets` 上，误差 <= 1e-10（§11.2 T03/T08）。"""
    device = baseline_run.device
    for attempt in baseline_run.successful:
        state = attempt.result.state
        targets = device.contact_targets(attempt.bias_V)
        assert abs(state.u[0] - targets[0]) <= CONTACT_TOL
        assert abs(state.u[-1] - targets[3]) <= CONTACT_TOL
        assert abs(np.log(state.n[0]) - np.log(targets[1])) <= CONTACT_TOL
        assert abs(-np.log(state.p[0]) + np.log(targets[2])) <= CONTACT_TOL
        assert abs(np.log(state.n[-1]) - np.log(targets[4])) <= CONTACT_TOL
        assert abs(-np.log(state.p[-1]) + np.log(targets[5])) <= CONTACT_TOL


def test_t09_total_current_is_conserved_along_the_device(baseline_run):
    """T09：`j_n + j_p` 沿所有边恒定，相对起伏 <= 1e-3（§11.2）。"""
    for attempt in baseline_run.successful:
        fluxes = attempt.result.physics.fluxes
        total = fluxes.jn_edges + fluxes.jp_edges
        # 0 V 的 J 恒为零；用 1e-9（约 4e-15 A/cm^2）的绝对下限吸收舍入噪声。
        scale = max(float(np.max(np.abs(total))), 1e-9)
        spread = float(np.ptp(total)) / scale
        assert spread <= CONSERVATION_TOL, (attempt.bias_V, spread)


def test_forward_current_is_exponential_with_ideality_near_one(baseline_run):
    """电学自洽性：正向电流随偏压指数上升，理想因子落在 [1.0, 1.5]。"""
    from minitcad.simulation.diagnostics import fields_for_attempt

    currents = [
        fields_for_attempt(baseline_run.config, baseline_run.device, attempt).terminal_current_density_A_cm2()
        for attempt in baseline_run.successful
        if attempt.bias_V >= 0.1
    ]
    biases = [attempt.bias_V for attempt in baseline_run.successful if attempt.bias_V >= 0.1]
    assert np.all(np.diff(currents) > 0.0)

    vt = baseline_run.device.thermal_voltage_v
    for index in range(1, len(currents)):
        ideality = (biases[index] - biases[index - 1]) / vt / np.log(currents[index] / currents[index - 1])
        assert 1.0 <= ideality <= 1.5, (biases[index], ideality)


def test_built_in_potential_matches_contact_potentials(baseline_run):
    """`V_bi = V_T (a(C_R) - a(C_L))`，与 CSV benchmark 行口径一致。"""
    from minitcad.simulation.diagnostics import benchmark_fields

    values = benchmark_fields(baseline_run.config, baseline_run.device)
    left = values["contact_potential_left_V"]
    right = values["contact_potential_right_V"]
    assert right > 0.0 > left
    assert values["built_in_potential_V"] == pytest.approx(right - left, rel=1e-12)
    assert 0.6 < values["built_in_potential_V"] < 1.0  # 1e16 / 2e15 硅 PN 的量级
