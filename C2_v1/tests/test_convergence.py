
"""W5 Day 2: shared QF convergence contract tests."""

import numpy as np
import pytest

from minitcad.core.contact import contact_targets_qf
from minitcad.core.state import State, QFState

from minitcad.discretization.residual import evaluate_qf

from minitcad.solvers.convergence import (
    scaled_residual_qf,
    scaled_update_physical,
)

from test_gummel import baseline


def test_physical_update_is_zero_for_identical_states():
    """Identical states must give zero physical update."""
    _, state = baseline()

    assert scaled_update_physical(
        state, state
    ) == 0.0


def test_physical_update_checks_densities_not_eta():
    """The update formula must operate on physical n/p."""
    old = State(
        u=[0, 0, 0],
        n=[1, 1, 1],
        p=[1, 1, 1],
    )

    new = State(
        u=[0, 0, 0],
        n=[1, 1.1, 1],
        p=[1, 1, 1],
    )

    assert scaled_update_physical(
        old, new
    ) == pytest.approx(0.1 / 1.1)


def test_n_floor_changes_denominator_not_concentration():
    """n_floor is a scale, never a density clipping rule."""
    old = State(
        u=[0, 0, 0],
        n=[1e-20] * 3,
        p=[1] * 3,
    )

    new = State(
        u=[0, 0, 0],
        n=[2e-20] * 3,
        p=[1] * 3,
    )

    norm = scaled_update_physical(
        old,
        new,
        n_floor=1e-15,
    )

    assert norm == pytest.approx(1e-5)

    # Concentration was never clipped to n_floor.
    assert new.n[0] == 2e-20


def test_residual_uses_qf_contact_targets_and_stable_fluxes():
    """Verify the boundary scale uses shifted QF targets."""
    device, state = baseline()

    bias = 0.05
    ref_n = 0.3
    ref_p = -0.4

    qf = QFState.from_physical(
        state,
        ref_n=ref_n,
        ref_p=ref_p,
    )

    ev = evaluate_qf(qf, device, bias)

    value = scaled_residual_qf(
        qf,
        device,
        bias,
        ev,
    )

    f = ev.f.reshape(-1, 3)

    targets = contact_targets_qf(
        device,
        bias,
        ref_n=ref_n,
        ref_p=ref_p,
    )

    left = np.max(
        np.abs(f[0])
        / np.maximum(1.0, np.abs(targets[:3]))
    )

    right = np.max(
        np.abs(f[-1])
        / np.maximum(1.0, np.abs(targets[3:]))
    )

    assert value >= left
    assert value >= right
    assert np.isfinite(value)


def test_bad_floor_is_rejected():
    """Invalid convergence settings must not be accepted."""
    _, state = baseline()

    with pytest.raises(ValueError):
        scaled_update_physical(
            state,
            state,
            n_floor=0.0,
        )


def test_zero_step_success_returns_stable_diagnostics():
    """An already converged state must retain full physics."""
    from minitcad.solvers.gummel import solve_gummel

    device, state = baseline()

    first = solve_gummel(
        device,
        0.05,
        state,
    )

    assert first.converged, first.reason

    repeat = solve_gummel(
        device,
        0.05,
        first.state,
    )

    assert repeat.converged
    assert repeat.iterations == 0
    assert repeat.update_norm == 0.0
    assert repeat.physics is not None

    m = device.mesh.node_count

    assert repeat.physics.residual_qf_nd.shape == (
        3 * m,
    )

    assert repeat.physics.fluxes.jn_edges.shape == (
        m - 1,
    )


def test_small_reverse_bias_has_valid_physics():
    """A mild reverse-bias case must also be verified."""
    from minitcad.solvers.gummel import solve_gummel

    device, state = baseline()

    result = solve_gummel(
        device,
        -0.005,
        state,
    )

    assert result.converged, result.reason
    assert result.physics is not None

    total = (
        result.physics.fluxes.jn_edges
        + result.physics.fluxes.jp_edges
    )

    variation = (
        np.ptp(total)
        / max(1.0, abs(np.median(total)))
    )

    assert variation <= 1e-3