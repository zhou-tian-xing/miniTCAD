"""W3-B3：RunResult 最小容器（第 10.2 节）。"""

from __future__ import annotations

import pytest

from minitcad.simulation.run_result import (
    DECISION_ACCEPTED,
    DECISION_REJECTED,
    Attempt,
    RunResult,
)


def test_attempt_validates_fields():
    with pytest.raises(ValueError):
        Attempt(bias_V=float("nan"), status="PASS")
    with pytest.raises(ValueError):
        Attempt(bias_V=0.0, status="OK")
    with pytest.raises(ValueError):
        Attempt(bias_V=0.0, status="PASS", runtime_s=-1.0)


def test_empty_result_is_rejected_by_default():
    result = RunResult(config=object(), device=object())
    assert result.decision == DECISION_REJECTED
    assert result.accepted is False
    assert result.attempts == () and result.successful == ()


def test_successful_must_be_a_subset_of_attempts():
    passed = Attempt(bias_V=0.0, status="PASS")
    failed = Attempt(bias_V=0.1, status="FAIL")
    result = RunResult(
        config=object(),
        device=object(),
        attempts=(passed, failed),
        successful=(passed,),
        decision=DECISION_ACCEPTED,
    )
    assert result.accepted is True
    with pytest.raises(ValueError):
        RunResult(
            config=object(),
            device=object(),
            attempts=(failed,),
            successful=(passed,),
        )
    with pytest.raises(ValueError):
        RunResult(
            config=object(),
            device=object(),
            attempts=(failed,),
            successful=(failed,),
            decision=DECISION_ACCEPTED,
        )


def test_attempt_keeps_qf_physics_diagnostics():
    """W3 QF 同步：已接受点的 SolveResult 必须保留 flux/rate/residual 诊断。"""
    import numpy as np

    from minitcad.core.state import (
        Fluxes,
        PhysicsDiagnostics,
        SolveResult,
        State,
    )

    state = State(u=np.zeros(4), n=np.ones(4), p=np.ones(4))
    diagnostics = PhysicsDiagnostics(
        fluxes=Fluxes(jn_edges=np.zeros(3), jp_edges=np.zeros(3)),
        recombination_nd=np.zeros(4),
        residual_qf_nd=np.zeros(12),
    )
    result = SolveResult(
        0.0, state, True, 2, 1e-9, 1e-7, "SUCCESS", physics=diagnostics
    )
    attempt = Attempt(bias_V=0.0, status="PASS", result=result)
    run = RunResult(
        config=object(),
        device=object(),
        attempts=(attempt,),
        successful=(attempt,),
        decision=DECISION_ACCEPTED,
    )
    assert run.successful[0].result.physics is diagnostics
    assert run.successful[0].result.physics.fluxes.jn_edges.size == 3
