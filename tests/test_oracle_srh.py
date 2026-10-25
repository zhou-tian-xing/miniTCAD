"""T02_SRH_SWAP：异寿命 SRH 的独立 oracle（第 2.3、11.2 节）。

用 tau_n != tau_p 才能检出“寿命乘错位置”：
n=10, p=2, tau_n=1, tau_p=4 -> (np-1)/(tau_p(n+1)+tau_n(p+1)) = 19/47。
乘反会得到 19/23。A 的 recombination.py 落地前整个模块 skip。
"""

from __future__ import annotations

from fractions import Fraction

import pytest

recombination = pytest.importorskip(
    "minitcad.physics.recombination",
    reason="A 的 W3 交付 recombination.py 尚未落地（W3-A2）",
)

CANDIDATE_NAMES = (
    "srh_rate_nd",
    "srh_recombination_nd",
    "recombination_rate_nd",
    "srh_nd",
    "srh_rate",
)


def rate(n: float, p: float, tau_n: float, tau_p: float) -> float:
    """按 (n, p, tau_n, tau_p) 的位置参数调用 A 的实现。"""
    for name in CANDIDATE_NAMES:
        function = getattr(recombination, name, None)
        if function is not None:
            return float(function(n, p, tau_n, tau_p))
    pytest.skip(
        "recombination.py 里没有已知的 SRH 入口 "
        f"{CANDIDATE_NAMES}；A 冻结签名后在此补一行适配"
    )


def test_asymmetric_lifetimes_pin_the_factor_order():
    expected = Fraction(19, 47)
    got = rate(n=10.0, p=2.0, tau_n=1.0, tau_p=4.0)
    assert abs(got - float(expected)) <= 1.0e-12
    # 乘反的实现会落在 19/23，必须能被这个测试区分出来
    assert abs(got - float(Fraction(19, 23))) > 1.0e-2


def test_equilibrium_state_gives_zero_rate():
    # np = 1 -> 分子为 0（质量作用定律）
    assert abs(rate(n=2.0, p=0.5, tau_n=1.0, tau_p=4.0)) <= 1.0e-14


def test_rate_sign_follows_injection():
    # 过剩载流子 np>1 -> 复合为正；抽取 np<1 -> 产生为负
    assert rate(n=3.0, p=3.0, tau_n=1.0, tau_p=4.0) > 0.0
    assert rate(n=0.3, p=0.3, tau_n=1.0, tau_p=4.0) < 0.0


def test_uses_dimensionless_normalization():
    # 无量纲式 r=(np-1)/(tau_p(n+1)+tau_n(p+1))，不是 (np-ni^2)/(...)
    got = rate(n=10.0, p=2.0, tau_n=2.0, tau_p=3.0)
    expected = Fraction(19, 3 * 11 + 2 * 3)
    assert abs(got - float(expected)) <= 1.0e-12
