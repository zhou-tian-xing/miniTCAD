"""W5：Newton 的适配层契约，以及它与 Gummel 在同一状态下的一致性。

T11（双求解器正式对照）在 W6，这里先把**求解器级别**的契约钉住：

* 两个求解器的签名逐参数一致——`simulation/dc.py` 只按 `method` 选入口，靠的就是这一点；
* 在 0.30 V 上两者收敛到同一个解（`|Δu|`、终端电流、内部诊断量都对得上）；
* 失败时 `state is None`、`physics is None`，`iteration_history` 仍可读；
* 参数非法抛 `ValueError`，不靠返回值表达编程错误。

这条测试也是 Gate 2 的**前置证据**：如果两个求解器连同一个解都算不到，T11 就没法解释。
"""

from __future__ import annotations

import inspect
import pathlib

import numpy as np
import pytest

from minitcad.core.device import build_device
from minitcad.core.state import PhysicsDiagnostics, SolveResult
from minitcad.discretization.residual import evaluate_qf
from minitcad.io.config import load_config
from minitcad.simulation.sweep import equilibrium_start_state, run_fixed_sweep
from minitcad.solvers import solve_gummel, solve_newton
from minitcad.solvers.convergence import scaled_residual_qf

EXAMPLES = pathlib.Path(__file__).resolve().parents[1] / "examples"
BASELINE = EXAMPLES / "pn_baseline_0_to_030.yaml"


@pytest.fixture(scope="module")
def setup():
    config = load_config(BASELINE)
    device = build_device(config)
    equilibrium = equilibrium_start_state(device)
    outcome = run_fixed_sweep(device, config.sweep.biases(), equilibrium, config.solver)
    assert outcome.reached_V == pytest.approx(0.3)
    return config, device, equilibrium, outcome.last_converged_state


def test_solver_signatures_match():
    """两个入口的参数名与默认值必须完全相同。"""
    gummel = inspect.signature(solve_gummel)
    newton = inspect.signature(solve_newton)
    assert list(gummel.parameters) == list(newton.parameters)
    for name, parameter in gummel.parameters.items():
        other = newton.parameters[name]
        assert parameter.kind == other.kind, name
        if parameter.default is not inspect.Parameter.empty:
            assert parameter.default == other.default, name


def test_newton_converges_at_forward_bias(setup):
    _, device, equilibrium, _ = setup
    result = solve_newton(device, 0.05, equilibrium)
    assert isinstance(result, SolveResult)
    assert result.converged and result.state is not None
    assert result.reason == "SUCCESS"
    assert result.residual_norm <= 1e-8
    assert result.update_norm <= 1e-6
    assert isinstance(result.physics, PhysicsDiagnostics)


def test_newton_and_gummel_reach_the_same_solution(setup):
    """同一偏压、同一初值：两个求解器的解必须一致到远优于 T11 的判据。"""
    config, device, equilibrium, _ = setup
    # 两者都从同一个 0.25 V 解出发，避免“热启动路径不同”被误读成“解不同”。
    warm = run_fixed_sweep(
        device, (0.0, 0.05, 0.1, 0.15, 0.2, 0.25), equilibrium, config.solver
    ).last_converged_state
    assert warm is not None

    reference = solve_gummel(device, 0.3, warm)
    subject = solve_newton(device, 0.3, warm)
    assert reference.converged and subject.converged

    assert np.max(np.abs(reference.state.u - subject.state.u)) < 1e-8
    relative = np.max(
        np.abs(reference.state.n - subject.state.n)
        / np.maximum(reference.state.n, subject.state.n)
    )
    assert relative < 1e-8

    # 终端电流：Gate 2 要求非零电流相对差 ≤ 1%，这里实测在 1e-9 量级。
    j_gummel = reference.physics.fluxes.jn_edges[-1] + reference.physics.fluxes.jp_edges[-1]
    j_newton = subject.physics.fluxes.jn_edges[-1] + subject.physics.fluxes.jp_edges[-1]
    assert abs(j_gummel - j_newton) <= 1e-6 * max(abs(j_gummel), 1e-30)

    # 势峰值的差（T11 用 1e-3 V 作判据），这里换算成 V。
    peak = float(np.max(np.abs(reference.state.u - subject.state.u))) * device.thermal_voltage_v
    assert peak <= 1e-3


def test_diagnostics_match_an_independent_reevaluation(setup):
    """`diagnostics` 是给验证器用的内部量，必须与公开 `evaluate_qf` 逐位一致。"""
    _, device, equilibrium, _ = setup
    store: dict = {}
    result = solve_newton(device, 0.05, equilibrium, diagnostics=store)
    assert result.converged

    from minitcad.core.state import QFState

    qf = QFState(
        u=store["u"], eta_n=store["eta_n"], eta_p=store["eta_p"],
        ref_n=store["ref_n"], ref_p=store["ref_p"],
    )
    fresh = evaluate_qf(state=qf, device=device, bias_V=0.05)
    assert np.array_equal(store["residual"], fresh.f.reshape(qf.node_count, 3))
    assert np.array_equal(store["jn_edges"], fresh.fluxes.jn_edges)
    assert np.array_equal(store["rate"], fresh.recombination_nd)
    assert scaled_residual_qf(qf, device, 0.05, fresh) == pytest.approx(result.residual_norm)


def test_failure_returns_a_result_instead_of_raising(setup):
    """迭代上限调成 1 时必然失败；失败必须走返回值，且不带任何状态。"""
    _, device, equilibrium, _ = setup
    result = solve_newton(device, 0.05, equilibrium, max_iter=1)
    assert not result.converged
    assert result.state is None and result.physics is None
    assert result.reason in {"MAX_ITER", "LINE_SEARCH_FAILED", "LINEAR_SOLVE_FAILED"}
    assert result.iteration_history


def test_linear_solve_failure_is_reported(setup, monkeypatch):
    """线性系统解不动时判 `LINEAR_SOLVE_FAILED`，不静默返回一个假解。"""
    from minitcad.solvers import newton as newton_module

    def broken(*args, **kwargs):
        raise FloatingPointError("injected linear failure")

    monkeypatch.setattr(newton_module, "jacobian_qf", broken)
    _, device, equilibrium, _ = setup
    result = solve_newton(device, 0.05, equilibrium)
    assert not result.converged and result.state is None
    assert result.reason == "LINEAR_SOLVE_FAILED"
    assert any(record["status"] == "LINEAR_SOLVE_FAILED" for record in result.iteration_history)


@pytest.mark.parametrize(
    "overrides",
    [
        {"tol_residual": 0.0},
        {"tol_update": -1.0},
        {"max_iter": 0},
        {"max_iter": 2.5},
        {"alpha_min": 2.0},
        {"n_floor": 0.0},
        {"bias_V": float("nan")},
    ],
)
def test_invalid_parameters_are_rejected(setup, overrides):
    _, device, equilibrium, _ = setup
    kwargs = dict(overrides)
    bias_V = kwargs.pop("bias_V", 0.05)
    with pytest.raises(ValueError):
        solve_newton(device, bias_V, equilibrium, **kwargs)


def test_iteration_history_has_the_documented_keys(setup):
    _, device, equilibrium, _ = setup
    result = solve_newton(device, 0.05, equilibrium)
    assert result.converged
    keys = {"iteration", "alpha", "residual_norm", "update_norm", "status"}
    for record in result.iteration_history:
        assert keys & set(record), record
        assert record["iteration"] >= 0
    assert {record["status"] for record in result.iteration_history} <= {
        "ACCEPTED", "REJECTED", "LINEAR_SOLVE_FAILED", "INVALID_INITIAL_STATE",
    }
