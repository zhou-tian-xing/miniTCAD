"""T01_BERNOULLI：B(z) 的独立 oracle（第 6.4、6.5、11.2 节）。

期望值由 60 位 Decimal 直接算 z/(e^z-1)，不复用被测实现；A 的
discretization/bernoulli.py 落地前整个模块 skip。
最终门槛是 Gate 1 的 1e-11，这里先按第 6.5 节更严的 3e-14 量级自测。
"""

from __future__ import annotations

from decimal import Decimal, localcontext

import numpy as np
import pytest

bernoulli = pytest.importorskip(
    "minitcad.discretization.bernoulli",
    reason="A 的 W2 交付 bernoulli.py 尚未落地（W2-A1）",
).bernoulli

TOL = 3.0e-14
Z_VALUES = [
    0.0,
    1e-14, -1e-14, 1e-12, -1e-12,
    9.99e-4, -9.99e-4, 1.001e-3, -1.001e-3,
    1.0, -1.0, 49.9, -49.9, 50.1, -50.1,
    700.0, -700.0, 1000.0, -1000.0,
]


def reference_decimal(z: float) -> Decimal:
    """B(z)=z/(e^z-1) 的 60 位参照；大 |z| 走渐近支路避开 exp 溢出。"""
    with localcontext() as ctx:
        ctx.prec = 60
        zd = Decimal(z)
        if zd == 0:
            return Decimal(1)
        if zd > 60:  # z e^-z / (1 - e^-z)
            e = (-zd).exp()
            return zd * e / (1 - e)
        if zd < -60:  # -z / (1 - e^z)
            return -zd / (1 - zd.exp())
        return zd / (zd.exp() - 1)


@pytest.mark.parametrize("z", Z_VALUES)
def test_matches_high_precision_reference(z):
    want = float(reference_decimal(z))
    got = bernoulli(z)
    assert got >= 0.0
    assert abs(got - want) <= TOL * abs(want) + 1.0e-300


@pytest.mark.parametrize("z", Z_VALUES)
def test_reflection_identity(z):
    # B(-z) = B(z) + z，独立于参照实现（第 6.5 节）；
    # 等价于 B(-z) - B(z) - z = 0（旧版本这里写成了 +z，会误报 2z 的偏差）。
    left = bernoulli(-z) - bernoulli(z) - z
    assert abs(left) <= 1.0e-12 * max(1.0, abs(z))


def test_zero_is_one():
    assert bernoulli(0.0) == pytest.approx(1.0, rel=0.0, abs=1e-15)


def test_large_positive_underflows_to_zero():
    # 允许的正确下溢，不得靠 clipping 改物理值（第 6.4 节）
    assert bernoulli(1000.0) == 0.0


def test_scalar_returns_python_float_and_array_keeps_shape():
    assert isinstance(bernoulli(0.5), float)
    matrix = np.zeros((2, 3))
    out = bernoulli(matrix)
    assert isinstance(out, np.ndarray) and out.shape == (2, 3) and out.dtype == np.float64
    assert bernoulli(np.array([])).size == 0


def test_input_array_is_not_modified():
    z = np.array([0.5, 1.5, -2.0])
    before = z.copy()
    bernoulli(z)
    assert np.array_equal(z, before)


@pytest.mark.parametrize("bad", [float("nan"), float("inf"), float("-inf")])
def test_non_finite_raises(bad):
    with pytest.raises(ValueError):
        bernoulli(bad)
