#Bernoulli函数测试
#python -m pytest
import pytest
import numpy as np
from minitcad import discretization
def test_known_values():
    assert discretization.bernoulli(0) == 1.0
    assert discretization.bernoulli(1000) == 0.0
    assert discretization.bernoulli(-1000) == 1000.0
    assert discretization.bernoulli(750) > 0
    assert discretization.bernoulli(1.0) == pytest.approx(
        0.5819767068693265,rel = 2e-15,abs = 0.0
    )
    assert discretization.bernoulli(-1.0) == pytest.approx(
        1.5819767068693265,rel = 2e-15,abs = 0.0
    )
#B(-z) == B(z) + z检验
def test_reflection_identity():
    z = np.array([
        -1000,-700,-50,-1,-1e-12,
        0,1e-12,1,50,700,1000
    ],dtype  = float)
    np.testing.assert_allclose(
        discretization.bernoulli(-z),
        discretization.bernoulli(z) + z,
        atol = 3e-14,
        rtol = 3e-14
    )
#伯努利函数的可靠性
def test_shapes_types_and_no_mutation():
    assert type(discretization.bernoulli(1.0)) == float
    assert type(discretization.bernoulli(np.float64(1.0))) == float
    arr = np.array([
        [-1000,-1,0],
        [1,50,1000]
    ],dtype = float
    )
    before = arr.copy()
    got = discretization.bernoulli(arr)
    assert isinstance(got,np.ndarray)
    assert got.shape == before.shape
    assert got.dtype == np.float64
    np.testing.assert_array_equal(
        arr,
        before
    )
    assert discretization.bernoulli(np.empty((0,2))).shape == (0,2)
#非有限输入测试
@pytest.mark.parametrize("bad",[
    float("nan"),
    float("inf"),
    -float("inf"),
    [0.0,float("nan")],
    np.array([1.0,float("inf")])
])
def test_nonfinite_rejected(bad):
    with pytest.raises(ValueError):
        discretization.bernoulli(bad)
#输入类型错误测试
@pytest.mark.parametrize("bad",[
    "not a number",
    1 + 2j
])
def test_nonreal_or_unconvertible_rejected(bad):
    with pytest.raises((ValueError,TypeError)):
        discretization.bernoulli(bad)