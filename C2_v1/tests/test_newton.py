
"""Week 5 Day 4 Newton smoke and interface tests."""
import numpy as np
import pytest

from minitcad.core.device import (
    DeviceConfig, MeshSpec, NormalizationConfig, build_device,
)
from minitcad.core.material import SILICON
from minitcad.core.region import Region
from minitcad.core.state import QFState
from minitcad.discretization.residual import evaluate_qf
from minitcad.simulation.equilibrium import solve_equilibrium_scalar
from minitcad.solvers import solve_newton


def example_case(nodes=61):
    device = build_device(DeviceConfig(
        length_cm=1e-3,
        area_cm2=1e-4,
        temperature_k=300.0,
        material=SILICON,
        regions=(
            Region(0, 2.5e-4, "silicon", 1e11, 0),
            Region(2.5e-4, 1e-3, "silicon", 0, 2e10),
        ),
        mesh=MeshSpec(kind="uniform", nodes=nodes),
        normalization=NormalizationConfig(
            reference_diffusivity_cm2_s=10.0
        ),
    ))
    t = device.contact_targets(0.0)
    eq = solve_equilibrium_scalar(
        device.mesh.x_nodes_nd, device.cbar_nd, t[0], t[3]
    )
    assert eq.converged
    return device, eq.state


def test_newton_equilibrium_zero_steps_has_physics():
    device, start = example_case()
    ans = solve_newton(device, 0.0, start)
    assert ans.converged, ans.reason
    assert ans.iterations == 0
    assert ans.state is not None
    assert ans.physics is not None
    assert ans.residual_norm <= 1e-8
    assert ans.update_norm == 0.0
    assert ans.physics.residual_qf_nd.shape == (3 * start.node_count,)
    assert ans.physics.fluxes.jn_edges.size == start.node_count - 1


@pytest.mark.parametrize("bias", [0.01, 0.05, -0.02])
def test_newton_runs_with_nonzero_bias(bias):
    device, start = example_case()
    ans = solve_newton(device, bias, start)
    assert ans.converged, (bias, ans.reason, ans.residual_norm)
    assert ans.iterations > 0
    assert ans.residual_norm <= 1e-8
    assert ans.update_norm <= 1e-6
    assert ans.state is not None and ans.physics is not None
    assert np.all(ans.state.n > 0)
    assert np.all(ans.state.p > 0)
    assert np.isfinite(ans.physics.fluxes.jn_edges).all()
    assert np.isfinite(ans.physics.fluxes.jp_edges).all()
    assert any(row["status"] == "ACCEPTED" for row in ans.iteration_history)
    t = device.contact_targets(bias)
    np.testing.assert_allclose(
        [ans.state.u[0], ans.state.n[0], ans.state.p[0],
         ans.state.u[-1], ans.state.n[-1], ans.state.p[-1]],
        t, rtol=1e-12, atol=1e-12,
    )
    ref_n = float(np.log(t[4]) - t[3])
    ref_p = float(-np.log(t[2]) - t[0])
    qf = QFState.from_physical(ans.state, ref_n=ref_n, ref_p=ref_p)
    again = evaluate_qf(qf, device, bias)
    np.testing.assert_allclose(
        ans.physics.residual_qf_nd, again.f, atol=1e-7, rtol=1e-7
    )
    np.testing.assert_allclose(
        ans.physics.fluxes.jn_edges, again.fluxes.jn_edges,
        atol=1e-7, rtol=1e-7,
    )
def test_newton_failure_returns_no_accepted_state():
    device, start = example_case()
    ans = solve_newton(device, 0.05, start, max_iter=1)
    assert not ans.converged
    assert ans.reason == "MAX_ITER"
    assert ans.state is None
    assert ans.physics is None
    assert len(ans.iteration_history) > 0
def test_newton_invalid_settings_are_rejected():
    device, start = example_case()
    with pytest.raises(ValueError):
        solve_newton(device, 0.05, start, max_iter=0)
    with pytest.raises(ValueError):
        solve_newton(device, 0.05, start, alpha_min=2.0)
    with pytest.raises(ValueError):
        solve_newton(device, np.nan, start)