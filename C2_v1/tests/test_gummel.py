"""Tests for QF-based Gummel solver."""

import numpy as np
import pytest

from minitcad import simulation, solvers

from minitcad.core.state import (
    State,
    QFState,
)

from minitcad.core.device import (
    build_device,
    DeviceConfig,
    MeshSpec,
    NormalizationConfig,
)

from minitcad.core.material import SILICON
from minitcad.core.region import Region


def baseline():
    device = build_device(
        DeviceConfig(
            length_cm=1e-3,
            area_cm2=1e-4,
            temperature_k=300.0,
            material=SILICON,
            regions=(
                Region(
                    0, 2.5e-4,
                    "silicon",
                    1e11, 0,
                ),
                Region(
                    2.5e-4, 1e-3,
                    "silicon",
                    0, 2e10,
                ),
            ),
            mesh=MeshSpec(
                kind="uniform",
                nodes=201,
            ),
            normalization=NormalizationConfig(
                reference_diffusivity_cm2_s=10.0
            ),
        )
    )

    t = device.contact_targets(0.0)

    eq = simulation.solve_equilibrium_scalar(
        device.mesh.x_nodes_nd,
        device.cbar_nd,
        t[0],
        t[3],
    )

    assert eq.converged
    return device, eq.state


def test_log_quasi_fermi_roundtrip():
    _, eq = baseline()

    qf = QFState.from_physical(
        eq,
        ref_n=0.3,
        ref_p=-0.2,
    )

    # Actual QF levels of an equilibrium state
    # should be approximately zero.
    assert np.max(
        np.abs(qf.eta_n + qf.ref_n)
    ) < 1e-12

    assert np.max(
        np.abs(qf.eta_p + qf.ref_p)
    ) < 1e-12

    reconstructed = qf.to_physical()

    np.testing.assert_allclose(
        reconstructed.n,
        eq.n,
        rtol=1e-14,
    )

    np.testing.assert_allclose(
        reconstructed.p,
        eq.p,
        rtol=1e-14,
    )


def test_gummel_zero_bias():
    device, eq = baseline()

    result = solvers.solve_gummel(
        device,
        0.0,
        eq,
    )

    assert result.converged
    assert result.residual_norm <= 1e-8
    assert result.update_norm <= 1e-6

    state = result.state

    assert isinstance(state, State)
    assert result.physics is not None

    # Equilibrium mass-action law
    assert np.max(
        np.abs(state.n * state.p - 1.0)
    ) < 1e-6

    # Use stable QF flux diagnostics
    flux = result.physics.fluxes

    assert max(
        np.max(np.abs(flux.jn_edges)),
        np.max(np.abs(flux.jp_edges)),
    ) < 1e-6


def test_gummel_forward_bias_conservation():
    device, eq = baseline()
    bias = 0.05

    result = solvers.solve_gummel(
        device,
        bias,
        eq,
    )

    assert result.converged, result.reason
    assert result.iterations > 0

    assert result.residual_norm <= 1e-8
    assert result.update_norm <= 1e-6

    state = result.state
    physics = result.physics

    assert physics is not None

    assert np.all(state.n > 0)
    assert np.all(state.p > 0)

    m = state.node_count

    assert physics.residual_qf_nd.shape == (
        3 * m,
    )

    # Six physical contact conditions
    targets = device.contact_targets(bias)

    actual_contacts = [
        state.u[0],
        state.n[0],
        state.p[0],
        state.u[-1],
        state.n[-1],
        state.p[-1],
    ]

    np.testing.assert_allclose(
        actual_contacts,
        targets,
        rtol=1e-12,
        atol=1e-12,
    )

    # Stable total edge current
    flux = physics.fluxes

    j_total = (
        flux.jn_edges
        + flux.jp_edges
    )

    variation = (
        np.ptp(j_total)
        / max(
            1.0,
            abs(np.median(j_total)),
        )
    )

    assert variation <= 1e-3

    # Independent discrete conservation identity:
    #
    # Fn + Fp = diff(jn + jp) / w
    #
    f = physics.residual_qf_nd.reshape(m, 3)

    np.testing.assert_allclose(
        f[1:-1, 1] + f[1:-1, 2],
        np.diff(j_total)
        / device.mesh.w_interior_nd,
        rtol=1e-11,
        atol=1e-10,
    )


def test_failure():
    device, eq = baseline()

    result = solvers.solve_gummel(
        device,
        0.05,
        eq,
        max_iter=1,
    )

    assert not result.converged
    assert result.state is None
    assert result.physics is None
    assert result.reason == "MAX_ITER"


def test_bad_input():
    device, eq = baseline()

    with pytest.raises(ValueError):
        solvers.solve_gummel(
            device,
            0.1,
            eq,
            max_iter=0,
        )

    with pytest.raises(ValueError):
        solvers.solve_gummel(
            device,
            0.1,
            eq,
            alpha_min=2.0,
        )

    with pytest.raises(ValueError):
        solvers.solve_gummel(
            device,
            0.1,
            eq,
            tol_residual=0,
        )

    with pytest.raises(ValueError):
        solvers.solve_gummel(
            device,
            float("nan"),
            eq,
        )
'''
#检验gummel迭代的结果
import pytest
import numpy as np
from minitcad import simulation,solvers,discretization
from minitcad.core.device import build_device,DeviceConfig,MeshSpec,NormalizationConfig
from minitcad.core.material import SILICON
from minitcad.core.region import Region
from minitcad.solvers.gummel import _quasi_fermi,_state
#基准的器件
def baseline():
    device = build_device(DeviceConfig(
        length_cm = 1e-3,
        area_cm2 = 1e-4,
        temperature_k = 300.0,
        material = SILICON,
        regions = (
            Region(0,2.5e-4,"silicon",1e11,0),
            Region(2.5e-4,1e-3,"silicon",0,2e10)
        ),
        mesh = MeshSpec(kind = "uniform",nodes = 201),
        normalization = NormalizationConfig(reference_diffusivity_cm2_s = 10.0)
    ))
    t = device.contact_targets(0)
    eq = simulation.solve_equilibrium_scalar(
        device.mesh.x_nodes_nd,device.cbar_nd,t[0],t[3])
    return device,eq.state
#准fermi变量测试
def test_log_quasi_fermi():
    _,eq = baseline()
    (fn,fp) = _quasi_fermi(eq)
    assert np.max(fn) < 1e-12
    assert np.max(fp) < 1e-12
    require = _state(eq.u,fn,fp)
    np.testing.assert_allclose(require.n,eq.n,rtol = 1e-14,atol = 0)
    np.testing.assert_allclose(require.p,eq.p,rtol = 1e-14,atol = 0)
#零偏压测试
def test_gummel_zero_bias():
    d,eq = baseline()
    result = solvers.solve_gummel(d,0.0,eq)
    assert result.converged
    assert result.residual_norm < 1e-8
    assert result.update_norm < 1e-6
    s = result.state
    assert np.max(abs(s.n * s.p - 1)) < 1e-6
    f = discretization.edge_flux(s,d)
    assert max(np.max(abs(f.jn_edges)),np.max(abs(f.jp_edges))) < 1e-6
#真实求解测试
def test_gummel_forward_bias_conservation():
    d,eq = baseline()
    bias = 0.05
    result = solvers.solve_gummel(d,bias,eq)
    assert result.converged
    assert result.iterations > 0
    assert result.residual_norm <= 1e-8
    assert result.update_norm <= 1e-6
    state = result.state
    assert (state.n > 0).all()
    assert (state.p > 0).all()
    assert np.max(
        np.abs(discretization.residual(state,d,bias))
    ) < 1e-6
    t = d.contact_targets(bias)
    np.testing.assert_allclose(
        [
            state.u[0],
            state.n[0],
            state.p[0],
            state.u[-1],
            state.n[-1],
            state.p[-1]
        ],
        t,
        rtol=1e-12,
        atol=1e-12
    )
    flow = discretization.edge_flux(state,d)
    j = flow.jn_edges + flow.jp_edges
    variation = (
        np.max(np.abs(j - np.median(j)))
        / max(1.0, abs(np.median(j)))
    )
    assert variation <= 1e-3
#求解失败测试
def test_failure():
    d,eq = baseline()
    r = solvers.solve_gummel(d,0.05,eq,max_iter = 1)
    assert not r.converged
    assert r.state is None
    assert r.reason == "MAX_ITER"
#异常输入测试
def test_bad_input():
    d,eq = baseline()
    with pytest.raises(ValueError):
        solvers.solve_gummel(d,0.1,eq,max_iter = 0)
    with pytest.raises(ValueError):
        solvers.solve_gummel(d,0.1,eq,alpha_min = 2.0)
'''