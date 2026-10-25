"""T06/T07：Scharfetter-Gummel 边通量的独立 oracle（第 6.3、11.2 节）。

A 在 2026-10-21 把求解器核心改成准费米能级（QF）表示后，SG 的入口拆成：

* ``electron_edge_flux(n_L, n_R, du, h, m_n)`` / ``hole_edge_flux(p_L, p_R, du, h, m_p)``
  —— 物理量（n、p）形式的边通量，T07 纯扩散极限用它；
* ``edge_flux_qf(u, eta_n, eta_p, h, mu_n, mu_p, ref_n=, ref_p=)``
  —— 生产版 QF 边通量，T06 平衡零电流用它（常数 QF ⇒ 严格零电流）。

两条恒等式的期望值由代数构造，而不是调用被测函数：
  T06 结点热平衡 n_R=n_L e^{du}, p_R=p_L e^{-du} -> j_n=j_p=0；
      常数 eta_n/eta_p -> j_n=j_p=0（QF 版，不依赖 exp 相减）；
  T07 du=0 -> 纯扩散 j_n=m_n(n_R-n_L)/h, j_p=m_p(p_L-p_R)/h。
"""

from __future__ import annotations

import numpy as np
import pytest

sg = pytest.importorskip(
    "minitcad.discretization.sg",
    reason="A 的 W3 交付 sg.py 尚未落地（W3-A1）",
)


def electron_flux(n_l, n_r, du, h, m_n):
    function = getattr(sg, "electron_edge_flux", None)
    if function is None:
        pytest.skip("sg.py 缺 electron_edge_flux；A 冻结签名后在此适配")
    return function(n_l, n_r, du, h, m_n)


def hole_flux(p_l, p_r, du, h, m_p):
    function = getattr(sg, "hole_edge_flux", None)
    if function is None:
        pytest.skip("sg.py 缺 hole_edge_flux；A 冻结签名后在此适配")
    return function(p_l, p_r, du, h, m_p)


def edge_flux_qf(u, eta_n, eta_p, h, mu_n, mu_p, ref_n=0.0, ref_p=0.0):
    function = getattr(sg, "edge_flux_qf", None)
    if function is None:
        pytest.skip("sg.py 缺 edge_flux_qf；A 冻结签名后在此适配")
    return function(u, eta_n, eta_p, h, mu_n, mu_p, ref_n=ref_n, ref_p=ref_p)


@pytest.mark.parametrize("du", [0.0, 0.37, -1.2, 3.7, 12.0, -20.0])
@pytest.mark.parametrize("n_l,p_l", [(1.0, 1.0), (1.0e3, 2.5e-4), (1.0e-2, 4.0e2)])
def test_equilibrium_fluxes_vanish(du, n_l, p_l):
    h = 0.013
    n_r = n_l * np.exp(du)
    p_r = p_l * np.exp(-du)
    jn = electron_flux(n_l, n_r, du, h, 1.0)
    jp = hole_flux(p_l, p_r, du, h, 1.0)
    assert abs(jn) <= 1.0e-11 * max(abs(n_l), abs(n_r)) / h
    assert abs(jp) <= 1.0e-11 * max(abs(p_l), abs(p_r)) / h


@pytest.mark.parametrize("ref_n,ref_p", [(0.0, 0.0), (0.3, -1.934), (-4.2, 7.5)])
@pytest.mark.parametrize("du", [0.0, 0.37, -1.2, 3.7, 12.0, -20.0])
def test_constant_quasi_fermi_gives_zero_current(ref_n, ref_p, du):
    """QF 版 T06：eta_n、eta_p 在整条边上为常数时电流必须严格为零。"""
    h = 0.013
    u = np.array([0.0, 0.5 * du, du])
    eta_n = np.full(3, ref_n - 0.7)
    eta_p = np.full(3, ref_p + 0.4)
    flux = edge_flux_qf(u, eta_n, eta_p, np.full(2, h), 1.0, 1.0,
                        ref_n=ref_n, ref_p=ref_p)
    assert np.all(flux.jn_edges == 0.0)
    assert np.all(flux.jp_edges == 0.0)


@pytest.mark.parametrize("h", [0.013, 0.25, 1.5])
def test_pure_diffusion_limit_has_closed_form(h):
    mn, mp = 1.0, 1.0
    n_l, n_r, p_l, p_r = 1.0, 3.0, 5.0, 2.0
    jn = electron_flux(n_l, n_r, 0.0, h, mn)
    jp = hole_flux(p_l, p_r, 0.0, h, mp)
    assert abs(jn - mn * (n_r - n_l) / h) <= 1.0e-11 * max(1.0, abs(jn))
    assert abs(jp - mp * (p_l - p_r) / h) <= 1.0e-11 * max(1.0, abs(jp))


def test_diffusion_limit_scales_with_mobility():
    h, du = 0.25, 0.0
    n_l, n_r, p_l, p_r = 1.0, 2.0, 2.0, 1.0
    jn_1 = electron_flux(n_l, n_r, du, h, 1.0)
    jp_1 = hole_flux(p_l, p_r, du, h, 1.0)
    jn_2 = electron_flux(n_l, n_r, du, h, 2.0)
    jp_2 = hole_flux(p_l, p_r, du, h, 3.0)
    assert jn_2 == pytest.approx(2.0 * jn_1, rel=1e-12)
    assert jp_2 == pytest.approx(3.0 * jp_1, rel=1e-12)


def test_equilibrium_holds_for_a_whole_vector_of_du():
    h = 0.1
    du = np.array([-8.0, -1.0, 0.0, 0.5, 4.0])
    n_l, p_l = 10.0, 0.1
    n_r = n_l * np.exp(du)
    p_r = p_l * np.exp(-du)
    jn = electron_flux(n_l, n_r, du, np.full(du.size, h), 1.0)
    jp = hole_flux(p_l, p_r, du, np.full(du.size, h), 1.0)
    scale_n = np.maximum(n_l, n_r) / h
    scale_p = np.maximum(p_l, p_r) / h
    assert np.all(np.abs(jn) <= 1.0e-9 * scale_n)
    assert np.all(np.abs(jp) <= 1.0e-9 * scale_p)


def test_qf_and_physical_paths_agree_between_nodes():
    """同一物理状态写成 QF 变量后，QF 边通量应与 n/p 版边通量一致。"""
    h = 0.04
    u = np.array([1.7, 2.53, 2.9])
    du = u[1] - u[0]
    ref_n, ref_p = 0.3, -0.5
    eta_n = np.array([0.35, -0.42, 0.1])
    eta_p = np.array([-0.18, 0.61, 0.2])
    n = np.exp(u + ref_n + eta_n)
    p = np.exp(-u - ref_p - eta_p)
    jn_phys = electron_flux(n[0], n[1], du, h, 2.5)
    jp_phys = hole_flux(p[0], p[1], du, h, 0.7)
    flux = edge_flux_qf(
        u, eta_n, eta_p, np.full(2, h), 2.5, 0.7, ref_n=ref_n, ref_p=ref_p
    )
    assert flux.jn_edges[0] == pytest.approx(jn_phys, rel=1e-12, abs=1e-15)
    assert flux.jp_edges[0] == pytest.approx(jp_phys, rel=1e-12, abs=1e-15)
    assert abs(jn_phys) > 1e-3
    assert abs(jp_phys) > 1e-3
