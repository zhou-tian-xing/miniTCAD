"""W2-B3 / I-01：``State`` 只读性与 3M 交错打包（第 7.1、10.2 节）。"""

from __future__ import annotations

import numpy as np
import pytest

from minitcad.core.state import Fluxes, SolveResult, State


def _state(m: int = 4) -> State:
    return State(
        u=np.linspace(0.0, 1.0, m),
        n=np.ones(m) * 2.0,
        p=np.ones(m) * 0.5,
    )


def test_pack_unpack_round_trip_preserves_values():
    state = _state(5)
    restored = State.unpack(state.pack(), state.node_count)
    assert np.array_equal(restored.u, state.u)
    assert np.array_equal(restored.n, state.n)
    assert np.array_equal(restored.p, state.p)


def test_pack_uses_interleaved_ordering():
    state = _state(3)
    x = state.pack()
    assert x.size == 9
    assert x[0] == state.u[0] and x[1] == state.n[0] and x[2] == state.p[0]
    assert x[3] == state.u[1] and x[4] == state.n[1] and x[5] == state.p[1]


def test_state_arrays_are_read_only():
    state = _state()
    with pytest.raises(ValueError):
        state.n[0] = 3.0  # I-01：跨模块原地修改必须被阻止


def test_state_defensively_copies_inputs():
    u = np.zeros(4)
    n = np.ones(4)
    p = np.ones(4)
    state = State(u=u, n=n, p=p)
    u[0] = 99.0
    n[0] = 99.0
    assert state.u[0] == 0.0 and state.n[0] == 1.0


def test_nonpositive_concentration_is_rejected():
    with pytest.raises(ValueError):
        State(u=np.zeros(3), n=np.array([1.0, 0.0, 1.0]), p=np.ones(3))
    with pytest.raises(ValueError):
        State(u=np.zeros(3), n=np.ones(3), p=np.array([1.0, -1.0, 1.0]))


def test_mismatched_lengths_are_rejected():
    with pytest.raises(ValueError):
        State(u=np.zeros(3), n=np.ones(4), p=np.ones(3))


def test_unpack_rejects_wrong_length():
    with pytest.raises(ValueError):
        State.unpack(np.zeros(5), 2)


def test_fluxes_share_length_and_are_read_only():
    fluxes = Fluxes(jn_edges=np.zeros(4), jp_edges=np.zeros(4))
    with pytest.raises(ValueError):
        fluxes.jn_edges[0] = 1.0
    with pytest.raises(ValueError):
        Fluxes(jn_edges=np.zeros(4), jp_edges=np.zeros(3))


def test_solve_result_rejects_invalid_fields():
    with pytest.raises(ValueError):
        SolveResult(0.0, _state(), True, -1, 0.0, 0.0, "SUCCESS")
    with pytest.raises(ValueError):
        SolveResult(0.0, _state(), False, 3, 0.0, 0.0, "")


def test_solve_result_keeps_iteration_history():
    history = ({"step": 0, "alpha": 1.0, "residual_norm": 1.0, "update_norm": 1.0},)
    result = SolveResult(0.1, _state(), True, 1, 1e-9, 1e-7, "SUCCESS", history)
    assert result.iteration_history[0]["alpha"] == 1.0
