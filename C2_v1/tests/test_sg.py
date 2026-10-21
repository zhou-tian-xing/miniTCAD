#Scharfetter-Gummel边通量测试
import pytest
import numpy as np
from types import SimpleNamespace
from minitcad import discretization
#归一化的器件
def device(h,m_n = 2.0,m_p = 3.0):
    return SimpleNamespace(
        mesh = SimpleNamespace(h_edges_nd = np.asarray(h,dtype = float)),
        m_n_nd = m_n,
        m_p_nd = m_p
    )
#电压为零的纯扩散测试
@pytest.mark.parametrize(
    "h,mn,mp,nl,nr,pl,pr",
    [
        (0.5, 2.0, 3.0, 1.0, 4.0, 5.0, 2.0),
        (0.1, 0.7, 1.6, 8.0, 2.0, 0.5, 2.5),
        (3.0, 2.4, 0.8, 7.0, 7.0, 1.0, 1.0),
    ]
)
def test_pure_diffusion(h,mn,mp,nl,nr,pl,pr):
    jn = discretization.electron_edge_flux(nl,nr,0.0,h,mn)
    jp = discretization.hole_edge_flux(pl,pr,0.0,h,mp)
    assert type(jn) == float and type(jp) == float
    assert jn == pytest.approx(mn * (nr - nl) / h,rel = 1e-12,abs = 0.0)
    assert jp == pytest.approx(mp * (pl - pr) / h,rel = 1e-12,abs = 0.0)
#Boltzman分布的稳定条件测试
@pytest.mark.parametrize(
    "du",
    [-20., -5., -1., -1e-12, 0., 1e-12, 1., 5., 20.]
)
def test_Boltzman_distribution(du):
    nl,pl = 2.0,0.5
    nr = nl * np.exp(du)
    pr = pl * np.exp(-du)
    jn = discretization.electron_edge_flux(nl,nr,du,0.3,1.7)
    jp = discretization.hole_edge_flux(pl,pr,du,0.3,0.8)
    assert abs(jn) <= 1e-11
    assert abs(jp) <= 1e-11
#浓度差为零的纯电流测试
@pytest.mark.parametrize(
    "du",
    [-1000.,1000.,1e-12,-1e-12]
)
def test_pure_drift(du):
    n,p,h,mn,mp = 2.0,5.0,0.5,1.5,0.7
    jn = discretization.electron_edge_flux(n, n, du, h, mn)
    jp = discretization.hole_edge_flux(p, p, du, h, mp)
    assert np.isfinite(jn).all() and np.isfinite(jp).all()
    assert jn == pytest.approx(-mn * n * du / h,rel = 1e-11,abs = 1e-11)
    assert jp == pytest.approx(-mp * p * du / h,rel = 1e-11,abs = 1e-11)
#非均匀边上两条通量，返回相同M-1的Fluxes，须和State一致
def test_edge_flux_shape_and_immutable():
    from minitcad.core.state import State,Fluxes
    state = State(u=[-0.5,0.0,0.8,1.0],n=[1,2,3,4],p=[5,4,3,2])
    obj = SimpleNamespace(mesh=SimpleNamespace(h_edges_nd=np.array([0.1,0.3,0.6])),
                          mu_n_nd=2.0,mu_p_nd=3.0)
    got = discretization.edge_flux(state,obj)
    assert isinstance(got,Fluxes)
    assert got.jn_edges.shape == (3,) and got.jp_edges.shape == (3,)
    with pytest.raises(ValueError):
        got.jn_edges[0] = 0
    with pytest.raises(ValueError):
        discretization.electron_edge_flux(-1,2,0,0.2,1)
    with pytest.raises(ValueError):
        discretization.hole_edge_flux(1,2,0,0.0,1)
    with pytest.raises(ValueError):
        discretization.hole_edge_flux(1,2,float('nan'),0.2,1)
@pytest.mark.parametrize("du",[-10.,-1.,-0.1,0.,0.1,1.,10.])
def test_swap_interval_reverses_both_currents(du):
    a = discretization.electron_edge_flux(2.,3.,du,0.7,1.2)
    b = discretization.electron_edge_flux(3.,2.,-du,0.7,1.2)
    c = discretization.hole_edge_flux(1.3,5.,du,0.7,3.)
    d = discretization.hole_edge_flux(5.,1.3,-du,0.7,3.)
    assert a == pytest.approx(-b,rel=2e-13,abs=1e-13)
    assert c == pytest.approx(-d,rel=2e-13,abs=1e-13)