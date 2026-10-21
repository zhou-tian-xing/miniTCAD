#测试平衡态的电压求解，vibecoding
import numpy as np
import pytest
from minitcad.core.state import State, SolveResult
from minitcad.discretization.poisson import poisson_residual_nd
from minitcad.simulation.equilibrium import solve_equilibrium_nd
def reference_device():
    q = 1.602176634e-19
    kb = 1.380649e-23
    eps = 11.7 * 8.8541878128e-14
    T = 300.0
    ni = 1e10
    VT = kb * T / q
    L0 = np.sqrt(eps * VT / (q * ni))
    x_cm = np.linspace(0, 1e-3, 201)
    xi = x_cm / L0
    cl = -1e16 / ni
    cr = 2e15 / ni
    x_cv_l = 0.5 * (
        x_cm[:-2] + x_cm[1:-1]
    )
    x_cv_r = 0.5 * (
        x_cm[1:-1] + x_cm[2:]
    )
    junction = 2.5e-4
    fraction_p = np.clip(
        (junction - x_cv_l) / (x_cv_r - x_cv_l),
        0.0,
        1.0,
    )
    cbar = fraction_p * cl + (1 - fraction_p) * cr
    u_left = np.arcsinh(cl / 2)
    u_right = np.arcsinh(cr / 2)
    return xi, cbar, u_left, u_right
def test_equilibrium_result():
    xi, cbar, ul, ur = reference_device()
    result = solve_equilibrium_nd(
        xi, cbar, ul, ur
    )
    assert isinstance(result, SolveResult)
    assert result.bias_V == 0.0
    assert result.converged
    assert result.reason == "SUCCESS"
    assert result.residual_norm <= 1e-8
    assert result.update_norm <= 1e-6
    state = result.state
    assert isinstance(state, State)
    assert state.node_count == 201
    assert not state.u.flags.writeable
    assert not state.n.flags.writeable
    assert not state.p.flags.writeable
    assert state.pack().shape == (603,)
    assert np.max(
        np.abs(state.n * state.p - 1)
    ) <= 1e-6
    f = poisson_residual_nd(
        xi=xi,
        u=state.u,
        n=state.n,
        p=state.p,
        cbar=cbar,
        u_left=ul,
        u_right=ur,
    )
    assert f.shape == (201,)
    assert np.isfinite(f).all()
    assert abs(f[0]) <= 1e-10
    assert abs(f[-1]) <= 1e-10
def test_failure_does_not_return_state():
    xi, cbar, ul, ur = reference_device()
    result = solve_equilibrium_nd(
        xi, cbar, ul, ur,
        max_iter=1,
    )
    assert not result.converged
    assert result.reason == "MAX_ITER"
    assert result.state is None
def test_invalid_mesh():
    xi = np.array([0.0, 0.5, 0.5, 1.0])
    cbar = np.array([-1.0, 1.0])
    with pytest.raises(
        ValueError,
        match="strictly increasing",
    ):
        solve_equilibrium_nd(
            xi, cbar, -1.0, 1.0
        )