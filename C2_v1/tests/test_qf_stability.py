#测试准费米能级方法的稳健性，vibecoding
import numpy as np
import pytest
from minitcad.core.device import (
    build_device,
    DeviceConfig,
    MeshSpec,
    NormalizationConfig,
)
from minitcad.core.region import Region
from minitcad.core.material import SILICON
from minitcad import simulation, solvers
from minitcad.discretization.sg import edge_flux_qf
from minitcad.discretization.residual import residual_qf
from minitcad.physics.recombination import (
    srh_recombination_qf,
)
def standard_pn():
    """高掺杂 PN 器件。"""
    device = build_device(DeviceConfig(
        length_cm=1e-3,
        area_cm2=1e-4,
        temperature_k=300.0,
        material=SILICON,
        regions=(
            Region(
                0.0, 2.5e-4,
                "silicon", 1e16, 0.0
            ),
            Region(
                2.5e-4, 1e-3,
                "silicon", 0.0, 2e15
            ),
        ),
        mesh=MeshSpec(
            kind="uniform",
            nodes=201,
        ),
        normalization=NormalizationConfig(
            reference_diffusivity_cm2_s=10.0
        ),
    ))
    t = device.contact_targets(0.0)
    eq = simulation.solve_equilibrium_scalar(
        device.mesh.x_nodes_nd,
        device.cbar_nd,
        t[0],
        t[3],
    )
    assert eq.converged
    return device, eq.state
def test_constant_quasi_fermi_zero_current():
    """常数准费米量必须严格对应零 SG 电流。"""
    u = np.linspace(-12.0, 13.0, 101)
    h = np.full(100, 0.01)
    eta_n = np.zeros(101)
    eta_p = np.zeros(101)
    flux = edge_flux_qf(
        u,
        eta_n,
        eta_p,
        h,
        3.0,
        1.0,
        ref_n=0.0,
        ref_p=-1.934,
    )
    assert np.count_nonzero(flux.jn_edges) == 0
    assert np.count_nonzero(flux.jp_edges) == 0
def test_qf_srh_expm1():
    """准费米 SRH 与直接公式的一般状态一致。"""
    u = np.array([0.2, 0.3])
    eta_n = np.array([0.1, 0.4])
    eta_p = np.array([0.2, 0.05])
    actual = srh_recombination_qf(
        u, eta_n, eta_p,
        2.0, 4.0,
        ref_n=0.0,
        ref_p=-0.3,
    )
    n = np.exp(u + eta_n)
    p = np.exp(-u + 0.3 - eta_p)
    expected = (
        (n * p - 1.0)
        / (4.0 * (n + 1.0) + 2.0 * (p + 1.0))
    )
    np.testing.assert_allclose(
        actual,
        expected,
        rtol=1e-14,
        atol=1e-15,
    )
@pytest.mark.parametrize(
    "bias",
    [0.0, 0.01, 0.05, 0.10],
)
def test_standard_pn_forward_bias(bias):
    """标准高掺杂器件：双收敛、正性、接触、电流守恒。"""
    device, eq = standard_pn()
    out = {}
    result = solvers.solve_gummel(
        device,
        bias,
        eq,
        max_iter=75,
        diagnostics=out,
    )
    assert result.converged, (
        bias,
        result.reason,
        result.residual_norm,
    )
    assert result.residual_norm <= 1e-8
    assert result.update_norm <= 1e-6
    state = result.state
    assert np.all(state.n > 0)
    assert np.all(state.p > 0)
    # 重新计算稳定内部残差
    f, flux, rate = residual_qf(
        device,
        bias,
        out["u"],
        out["eta_n"],
        out["eta_p"],
        ref_n=out["ref_n"],
        ref_p=out["ref_p"],
    )
    np.testing.assert_allclose(
        f,
        out["residual"],
        rtol=0,
        atol=0,
    )
    # 检查总电流守恒
    j = flux.jn_edges + flux.jp_edges
    variation = (
        np.ptp(j)
        / max(1.0, abs(np.median(j)))
    )
    assert variation <= 1e-3
    # 检查六个接触条件
    t = device.contact_targets(bias)
    np.testing.assert_allclose(
        [
            state.u[0],
            state.n[0],
            state.p[0],
            state.u[-1],
            state.n[-1],
            state.p[-1],
        ],
        t,
        rtol=1e-12,
        atol=1e-12,
    )
def test_forward_sweep_state_warmstart():
    """检查由上一偏压解进行下一步初始化。"""
    device, state = standard_pn()
    for bias in np.arange(0.0, 0.101, 0.01):
        result = solvers.solve_gummel(
            device,
            float(bias),
            state,
            max_iter=75,
        )
        assert result.converged, (
            bias,
            result.reason,
            result.residual_norm,
        )
        state = result.state