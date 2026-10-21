#中能级SRH复合测试
import numpy as np
import pytest
from minitcad import physics
#数值正确性测试
def test_right_number():
    result = physics.srh_recombination_nd(n = 10.0,p = 2.0,tau_n_nd = 1.0,tau_p_nd = 4.0)
    assert isinstance(result,float)
    assert result == pytest.approx(19.0 / 47.0,rel = 0.0,abs = 1e-12)
#稳定条件测试
@pytest.mark.parametrize(
    "n,p",
    [(1.0,1.0),(2.0,0.5),(0.5,2.0)]
)
def test_equilibrium(n,p):
    equilibrium_result = physics.srh_recombination_nd(n,p,1.0,4.0)
    assert isinstance(equilibrium_result,float)
    assert equilibrium_result == pytest.approx(0.0,rel = 0.0,abs = 1e-12)
#tau的等比例缩放测试
def test_lifetime_scaling():
    base = physics.srh_recombination_nd(10.0,2.0,1.0,4.0)
    scaled = physics.srh_recombination_nd(10.0,2.0,7.0,28.0)
    assert scaled == pytest.approx(base / 7.0)
#输出类型以及不改变输出测试
def test_right_format():
    n = np.array([10.0,0.25,1.0])
    n_copy = n.copy()
    p = np.array([2.0,2.0,1.0])
    result = physics.srh_recombination_nd(n,p,1.0,4.0)
    assert isinstance(result,np.ndarray) and result.shape == (3,)
    np.testing.assert_allclose(result,[19.0 / 47.0,-1.0 / 16.0,0.0],atol = 1e-14,rtol = 0.0)
    np.testing.assert_array_equal(n,n_copy)
#空输入测试
def test_empty_input():
    result = physics.srh_recombination_nd(np.array([]),np.array([]),1.0,4.0)
    assert result.shape == (0,)
#异常输入测试
@pytest.mark.parametrize(("field","bad"),[
    ("n",[0.0,-1.0,float("nan"),float("inf")]),
    ("p",[0.0,-1.0,float("nan"),float("inf")]),
    ("tau_n_nd",[0.0,-1.0,float("nan"),float("inf")]),
    ("tau_p_nd",[0.0,-1.0,float("nan"),float("inf")])
])
def test_invalid_input(field,bad):
    args = dict(
        n = 10.0,
        p = 2.0,
        tau_n_nd = 1.0,
        tau_p_nd = 4.0
    )
    args[field] = bad
    with pytest.raises(ValueError):
        physics.srh_recombination_nd(**args)