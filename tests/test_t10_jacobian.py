"""T10：用独立的 FD oracle 验证 `jacobian_qf`（§14.3）。

判据：0 V 与 0.30 V 上随机方向 `d` 满足
`‖J_A d − J_FD d‖∞ / max(1, ‖J_FD d‖∞) ≤ 1e-4`。

取方向导数而不是逐元素误差：Jacobian 的元素跨 `1e12` 量级，逐元素比较会被小元素主导，
而 Newton 消费的是 `J d` 这个向量。
"""

from __future__ import annotations

import pathlib

import numpy as np
import pytest

from minitcad.core.contact import contact_targets_qf
from minitcad.core.device import build_device
from minitcad.core.state import QFState
from minitcad.io.config import load_config
from minitcad.simulation.sweep import equilibrium_start_state, run_fixed_sweep
from minitcad.solvers.jacobian import jacobian_qf
from minitcad.verification import fd_jacobian_oracle
from minitcad.verification.fd_jacobian_oracle import (
    CENTRAL,
    fd_jacobian,
    relative_difference,
)

EXAMPLES = pathlib.Path(__file__).resolve().parents[1] / "examples"
BASELINE = EXAMPLES / "pn_baseline_0_to_030.yaml"

# 判据本身写在 §14.3；这里不调松，也不按实测值“贴边”。
TOLERANCE = 1e-4
SEED = 20261110


@pytest.fixture(scope="module")
def setup():
    """按 Gummel 走到 0.30 V，拿到 0 V 与 0.30 V 两个真解作为差分基点。"""
    config = load_config(BASELINE)
    device = build_device(config)
    equilibrium = equilibrium_start_state(device)
    outcome = run_fixed_sweep(device, config.sweep.biases(), equilibrium, config.solver)
    assert outcome.reached_V == pytest.approx(0.3)
    return device, equilibrium, outcome.last_converged_state


def _qf_state(device, physical, bias_V: float) -> QFState:
    """用求解器同一组锚点（右端电子、左端空穴）构造 QF 状态。

    锚点不一致会让 `eta` 整体平移，与 `jacobian_qf` 内部的 `ref_n/ref_p` 对不上。
    """
    targets = device.contact_targets(bias_V)
    ref_n = float(np.log(targets[4]) - targets[3])
    ref_p = float(-np.log(targets[2]) - targets[0])
    return QFState.from_physical(physical, ref_n=ref_n, ref_p=ref_p)


@pytest.mark.parametrize("bias_V", [0.0, 0.3])
def test_jacobian_matches_fd_oracle(setup, bias_V):
    device, equilibrium, forward = setup
    physical = equilibrium if bias_V == 0.0 else forward
    qf = _qf_state(device, physical, bias_V)

    fd = fd_jacobian(qf, device, bias_V)
    # 出现单边说明有试探点越界，此时 1e-4 判据会比真实精度更松，要先看见。
    assert len(fd.one_sided_columns()) == 0
    assert all(status == CENTRAL for status in fd.statuses)

    analytic = jacobian_qf(state=qf, device=device, bias_V=bias_V)
    rng = np.random.default_rng(SEED)
    direction = rng.standard_normal(qf.pack_qf().size)

    reference = fd.matrix @ direction
    measured = analytic @ direction
    assert relative_difference(measured, reference) <= TOLERANCE


def test_jacobian_matches_fd_oracle_for_several_directions(setup):
    """单个方向的运气成分排除掉：5 个独立方向都要在判据内。"""
    device, _, forward = setup
    bias_V = 0.3
    qf = _qf_state(device, forward, bias_V)
    fd = fd_jacobian(qf, device, bias_V)
    analytic = jacobian_qf(state=qf, device=device, bias_V=bias_V)

    rng = np.random.default_rng(SEED + 1)
    errors = [
        relative_difference(analytic @ d, fd.matrix @ d)
        for d in (rng.standard_normal(qf.pack_qf().size) for _ in range(5))
    ]
    assert max(errors) <= TOLERANCE


def test_oracle_only_depends_on_the_residual():
    """oracle 不允许 import 生产实现：否则 T10 退化成“拿 A 的代码验 A 的代码”。"""
    source = pathlib.Path(fd_jacobian_oracle.__file__).read_text(encoding="utf-8")
    code_lines = [
        line for line in source.splitlines() if line.strip().startswith(("import ", "from "))
    ]
    assert code_lines, "oracle 应当有显式 import 行"
    assert not any("solvers" in line for line in code_lines), code_lines
    assert not any("jacobian_qf" in line for line in code_lines), code_lines


def test_one_sided_and_failed_columns_are_reported(monkeypatch):
    """越界时不静默混用差分：单侧失败标 `ONE_SIDED_PLUS`，两侧都失败要抛。

    直接替换 `_trial` 来模拟越界——真造一个只让单侧溢出的状态更脆。
    调用顺序是「基点，然后逐列 `+h`、`-h`」，所以第 5 次是第 1 列的 `-h` 侧，
    第 6、7 次是第 2 列的两侧。
    """
    device = build_device(load_config(BASELINE))
    qf = _qf_state(device, equilibrium_start_state(device), 0.0)
    real_trial = fd_jacobian_oracle._trial
    calls = {"count": 0}

    def fake_trial(state, device_, bias_V, z):
        calls["count"] += 1
        if calls["count"] == 5 or calls["count"] in (6, 7):
            return None
        return real_trial(state, device_, bias_V, z)

    monkeypatch.setattr(fd_jacobian_oracle, "_trial", fake_trial)
    with pytest.raises(fd_jacobian_oracle.FdOracleError):
        fd_jacobian(qf, device, 0.0)

    calls["count"] = 0
    result = fd_jacobian(qf, device, 0.0, allow_failed=True)
    assert result.statuses[1] == fd_jacobian_oracle.ONE_SIDED_PLUS
    assert result.statuses[2] == fd_jacobian_oracle.FAILED
    assert result.one_sided_columns() == (1,)
    assert result.failed_columns() == (2,)
    assert np.isnan(result.matrix[:, 2]).all()


def test_relative_difference_denominator_is_guarded():
    """零参考向量时用 `max(1, ...)` 兜底，避免 0/0。"""
    zeros = np.zeros(3)
    assert relative_difference(zeros, zeros) == 0.0
    assert relative_difference(np.array([1e-9, 0.0, 0.0]), zeros) == pytest.approx(1e-9)


def test_boundary_qf_targets_satisfy_local_equilibrium(setup):
    """`contact_targets_qf` 在左端给出 `eta_n = eta_p = -v`，右端为 0。

    锚点约定被改动时这里先失败，不用等 T11 的双求解器对照出问题。
    """
    device, _, _ = setup
    bias_V = 0.25
    v = bias_V / device.thermal_voltage_v
    targets = contact_targets_qf(device, bias_V, ref_n=0.0, ref_p=0.0)
    assert targets[1] == pytest.approx(-v)  # (u_L, eta_nL, eta_pL, u_R, eta_nR, eta_pR)
    assert targets[2] == pytest.approx(-v)
    assert targets[4] == pytest.approx(0.0)
    assert targets[5] == pytest.approx(0.0)
