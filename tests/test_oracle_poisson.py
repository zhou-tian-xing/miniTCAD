"""T05_POISSON_NONUNIFORM：非均匀网格 Poisson 的独立 oracle（第 17.2、11.2 节）。

u=xi^2 时 d2u/dxi2=2，配合 cbar=-2、n=p=1 的 FV 残差在每个内部结点都应为 0，
且与网格是否均匀无关。这同时钉死 w_i=(h_{i-1}+h_i)/2 的下标约定。
A 的 poisson.py 落地前整个模块 skip。
"""

from __future__ import annotations

import numpy as np
import pytest

poisson = pytest.importorskip(
    "minitcad.discretization.poisson",
    reason="A 的 W2 交付 poisson.py 尚未落地（W2-A3）",
)

GRIDS = [
    [0.0, 0.25, 0.75, 1.5],
    [0.0, 0.1, 0.35, 0.8, 1.3, 1.5],
    [0.0, 0.02, 0.05, 0.11, 0.4, 1.0, 1.5],
]


def residual(xi, u, n, p, cbar, u_left, u_right):
    function = getattr(poisson, "poisson_residual_nd", None)
    if function is None:
        pytest.skip("poisson.py 缺 poisson_residual_nd；A 冻结签名后在此适配")
    return function(xi, u, n, p, cbar, u_left, u_right)


@pytest.mark.parametrize("xi", GRIDS)
def test_quadratic_manufactured_solution_is_exact_on_any_grid(xi):
    xi = np.asarray(xi, dtype=float)
    u = xi**2
    n = np.ones_like(xi)
    p = np.ones_like(xi)
    cbar = np.full_like(xi, -2.0)  # -d2u/dxi2 = -2
    f = residual(xi, u, n, p, cbar, float(u[0]), float(u[-1]))
    assert f.shape == xi.shape
    assert np.max(np.abs(f[1:-1])) <= 1e-11
    assert abs(f[0]) <= 1e-15 and abs(f[-1]) <= 1e-15


@pytest.mark.parametrize("xi", GRIDS)
def test_wrong_denominator_is_detected(xi):
    """用单边 h 代替 (h_{i-1}+h_i)/2 必须留下可见残差，防止回归。"""
    xi = np.asarray(xi, dtype=float)
    u = xi**2
    h = np.diff(xi)
    gradient = np.diff(u) / h
    wrong = np.diff(gradient) / h[:-1] - 2.0  # 错误下标：h_i 而不是 CV 宽度
    assert np.max(np.abs(wrong)) > 1e-2


def test_linear_solution_and_dirichlet_rows():
    xi = np.array([0.0, 0.4, 0.9, 1.5])
    u = 0.3 + 2.0 * xi  # 常数梯度 -> 内部残差为 0
    n = p = np.ones_like(xi)
    cbar = np.zeros_like(xi)
    f = residual(xi, u, n, p, cbar, u[0] - 1e-3, u[-1] + 2e-3)
    assert np.max(np.abs(f[1:-1])) <= 1e-11
    assert f[0] == pytest.approx(1e-3, abs=1e-15)
    assert f[-1] == pytest.approx(-2e-3, abs=1e-15)
