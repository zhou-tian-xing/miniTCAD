
"""A-side W7 regression tests.

These tests support A's core acceptance.
B must perform independent Gate verification.
"""

import importlib

import numpy as np
import pytest

from test_stability import make_case
from minitcad.solvers import (
    solve_gummel,
    solve_newton,
)


def terminal_j(result, device):
    """Return terminal current density [A/cm^2]."""
    assert result.converged
    assert result.physics is not None

    flux = result.physics.fluxes

    return float(
        (
            flux.jn_edges[0]
            + flux.jp_edges[0]
        ) * device.scales.j0_a_cm2
    )


def check_physics(result, device):
    """Check the public production-solver contract."""
    assert result.converged, result.reason
    assert result.state is not None
    assert result.physics is not None

    assert result.residual_norm <= 1e-8
    assert result.update_norm <= 1e-6

    m = device.mesh.node_count
    physics = result.physics

    assert physics.fluxes.jn_edges.shape == (m - 1,)
    assert physics.fluxes.jp_edges.shape == (m - 1,)
    assert physics.recombination_nd.shape == (m,)
    assert physics.residual_qf_nd.shape == (3 * m,)

    assert np.isfinite(
        physics.residual_qf_nd
    ).all()

    j = (
        physics.fluxes.jn_edges
        + physics.fluxes.jp_edges
    )

    # Dimensionless current-conservation metric.
    nonconservation = (
        np.max(np.abs(j - np.median(j)))
        / max(1.0, abs(np.median(j)))
    )

    assert nonconservation <= 1e-3

    return j


@pytest.mark.parametrize(
    "doping",
    [1e12, 1e14, 1e16],
)
def test_t11_forward_sweep_same_start_and_tolerance(
    doping,
):
    """Compare independent Gummel/Newton continuations."""
    device, initial = make_case(doping=doping)

    g_state = initial
    n_state = initial

    for bias in (
        0.05, 0.10, 0.15,
        0.20, 0.25, 0.30
    ):
        g = solve_gummel(device, bias, g_state)
        n = solve_newton(device, bias, n_state)

        check_physics(g, device)
        check_physics(n, device)

        potential_error_v = (
            np.max(
                np.abs(g.state.u - n.state.u)
            ) * device.thermal_voltage_v
        )

        assert potential_error_v <= 1e-3

        jg = terminal_j(g, device)
        jn = terminal_j(n, device)

        # Relative comparison when current is nonzero.
        # Otherwise use a protected absolute comparison.
        if max(abs(jg), abs(jn)) > (
            device.scales.j0_a_cm2 * 1e-8
        ):
            assert (
                abs(jg - jn)
                / max(abs(jg), abs(jn))
                <= 0.01
            )
        else:
            assert (
                abs(jg - jn)
                / device.scales.j0_a_cm2
                <= 0.01
            )

        # Never continue from a failed result.
        g_state = g.state
        n_state = n.state


def test_nonuniform_121_nodes_two_solvers():
    """Full PN solve on the existing segmented mesh."""
    device, initial = make_case(
        doping=1e14,
        segmented=True,
    )

    assert device.mesh.node_count == 121

    ratio = (
        device.mesh.h_edges_cm.max()
        / device.mesh.h_edges_cm.min()
    )
    assert ratio > 2.0

    for solver in (
        solve_gummel,
        solve_newton,
    ):
        result = solver(
            device,
            0.05,
            initial,
        )
        check_physics(result, device)


def test_nonuniform_high_doping_small_step_continuation():
    """Regression for the difficult nonuniform case."""
    device, state = make_case(
        doping=1e16,
        segmented=True,
    )

    for bias in (
        0.025, 0.05, 0.075, 0.10,
        0.125, 0.15, 0.175, 0.20
    ):
        result = solve_newton(
            device,
            bias,
            state,
        )

        check_physics(result, device)

        # Only use a validated successful state.
        state = result.state


@pytest.mark.parametrize(
    "module_name,solver_name",
    [
        (
            "minitcad.solvers.gummel",
            "solve_gummel",
        ),
        (
            "minitcad.solvers.newton",
            "solve_newton",
        ),
    ],
)
def test_diagnostics_match_last_qf_evaluation(
    monkeypatch,
    module_name,
    solver_name,
):
    """Capture the actual final QF evaluation.

    This does not rebuild currents from physical
    concentrations through the legacy SG path.
    """
    device, initial = make_case(
        doping=1e14,
        nodes=101,
    )

    module = importlib.import_module(module_name)
    original = module.evaluate_qf

    snapshots = []

    def spy(*args, **kwargs):
        qf = (
            kwargs["state"]
            if "state" in kwargs
            else args[0]
        )

        evaluation = original(*args, **kwargs)
        snapshots.append((qf, evaluation))

        return evaluation

    monkeypatch.setattr(
        module,
        "evaluate_qf",
        spy,
    )

    result = getattr(
        module,
        solver_name,
    )(device, 0.05, initial)

    assert result.converged

    last_qf, last_ev = snapshots[-1]

    np.testing.assert_array_equal(
        result.state.pack(),
        last_qf.to_physical().pack(),
    )

    np.testing.assert_array_equal(
        result.physics.residual_qf_nd,
        last_ev.f,
    )

    np.testing.assert_array_equal(
        result.physics.fluxes.jn_edges,
        last_ev.fluxes.jn_edges,
    )

    np.testing.assert_array_equal(
        result.physics.fluxes.jp_edges,
        last_ev.fluxes.jp_edges,
    )

    np.testing.assert_array_equal(
        result.physics.recombination_nd,
        last_ev.recombination_nd,
    )


def test_zero_bias_public_handoff_from_equilibrium():
    """0 V production output must contain diagnostics."""
    device, eq = make_case(doping=1e14)

    for solver in (
        solve_gummel,
        solve_newton,
    ):
        result = solver(
            device,
            0.0,
            eq,
        )

        check_physics(result, device)
        assert result.iterations == 0


def test_failed_attempt_not_a_production_result():
    """A failed state must never be passed to B as success."""
    device, initial = make_case(doping=1e14)

    failed = solve_newton(
        device,
        0.05,
        initial,
        max_iter=1,
    )

    assert not failed.converged
    assert failed.state is None
    assert failed.physics is None

    assert failed.reason != "SUCCESS"
    assert failed.iteration_history