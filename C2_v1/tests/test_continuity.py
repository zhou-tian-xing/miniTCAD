#连续性方程的测试
from types import SimpleNamespace
import numpy as np
import pytest
from minitcad.core.state import Fluxes
from minitcad import discretization
#构造辅助函数
def make_fluxes(jn, jp):
    return Fluxes(
        jn_edges = np.asarray(jn),
        jp_edges = np.asarray(jp)
    )
#非均匀网格精确值测试
def test_nonuniform_node():
    f = make_fluxes([2.0,5.0,11.0,20.0],[9.0,7.0,3.0,0.0])
    r = [100.0,1.0,2.0,3.0,-100.0]
    w = [1.5,2.5,3.5]
    (fn,fp) = discretization.continuity_residual_nd(f,r,w)
    assert fn.shape == fp.shape == (3,)
    np.testing.assert_allclose(
        fn,[1.0,0.4,-3.0 / 7.0],rtol = 0,atol = 1e-14
    )
    np.testing.assert_allclose(
        fp,[-1.0 / 3.0,0.4,15.0 / 7.0],rtol = 0,atol = 1e-14
    )
#边界条件忽略测试
def test_boundary_srh_values_do_not_enter_internal_equations():
    f = make_fluxes([1.0,2.0,4.0],[3.0,2.0,0.0])
    r_a = [0.0,0.25,-0.5,0.0]
    r_b = [1e12,0.25,-0.5,-1e12]
    w = [0.5,1.5]
    (fn_a,fp_a) = discretization.continuity_residual_nd(f,r_a,w)
    (fn_b,fp_b) = discretization.continuity_residual_nd(f,r_b,w)
    np.testing.assert_array_equal(fn_a,fn_b)
    np.testing.assert_array_equal(fp_a,fp_b)
#稳态解测试
def test_zero_srh_and_constant_individual_currents_yield_zero_residual():
    f = make_fluxes([2.0,2.0,2.0,2.0],[-3.0,-3.0,-3.0,-3.0])
    (fn,fp) = discretization.continuity_residual_nd(f,np.zeros(5),[0.5,1.0,2.0])
    np.testing.assert_array_equal(fn,[0.0,0.0,0.0])
    np.testing.assert_array_equal(fp,[0.0,0.0,0.0])
#符号正确测试
def test_srh_signs_are_opposite_in_electron_and_hole_equations():
    f = make_fluxes([0,0,0],[0,0,0])
    (fn,fp) = discretization.continuity_residual_nd(f,[0.0,0.25,-0.5,0.0], [0.5,1.5])
    np.testing.assert_array_equal(fn,[-0.25,0.5])
    np.testing.assert_array_equal(fp,[0.25,-0.5])
#代数守恒率：fn与fp相加消去复合率
def test_current_conservation():
    f = make_fluxes([0.2,2.3,-1.5,4.0],[3.5,-0.9,0.7,1.2])
    r = np.array([9.0,0.7,-0.3,1.1,-9.0])
    w = np.array([0.3,0.8,1.4])
    (fn,fp) = discretization.continuity_residual_nd(f,r,w)
    total_edges = f.jn_edges + f.jp_edges
    np.testing.assert_allclose(fn + fp,np.diff(total_edges) / w,rtol = 1e-13,atol = 1e-13)
#不改变原输入测试
def test_inputs_are_not_modified():
    jn = np.array([1.0,2.0,3.0])
    jp = np.array([4.0,5.0,6.0])
    r = np.array([9.0,0.2,0.3,8.0])
    w = np.array([0.5,0.6])
    arrays = (jn,jp,r,w)
    copies = tuple(arr.copy() for arr in arrays)
    discretization.continuity_residual_nd(make_fluxes(jn,jp),r,w)
    for actual,expected in zip(arrays,copies):
        np.testing.assert_array_equal(actual,expected)
#异常输入解模板
def invalid_fluxes(jn,jp):
    return SimpleNamespace(jn_edges = np.asarray(jn),jp_edges = np.asarray(jp))
#异常输入的数值测试
@pytest.mark.parametrize("fluxes,r,w",[
    (invalid_fluxes([1,np.nan],[2,1]),[0,0,0],[1]),
    (invalid_fluxes([1,np.inf],[2,1]),[0,0,0],[1]),
    (invalid_fluxes([1,2,],[2,np.inf]),[0,0,0],[1]),
    (invalid_fluxes([1,2,],[2,np.nan]),[0,0,0],[1]),
    (make_fluxes([1,2],[2,1]),[0.0,0.0,0.0],[True]),
    (make_fluxes([1,2],[2,1]),[0.0,0.0,1+2j],[1]),
    (make_fluxes([1,2],[2,1]),[0.0,0.0,0.0],["1"])
])
def test_invaild_input_numercial(fluxes,r,w):
    with pytest.raises(ValueError):
        discretization.continuity_residual_nd(fluxes,r,w)
#异常输入的形状测试
@pytest.mark.parametrize("fluxes,r,w",[
    (make_fluxes([1.0,2.0,3.0],[4.0,5.0,6.0]),[9.0,0.2,0.3,8.0],[0.5,0.6,0.0]),
    (make_fluxes([1.0,2.0,3.0],[4.0,5.0,6.0]),[9.0,0.2,0.3,8.0,7.5],[0.5,0.6]),
    (invalid_fluxes([1.0,2.0,3.0,0.0],[4.0,5.0,6.0]),[9.0,0.2,0.3,8.0],[0.5,0.6]),
    (invalid_fluxes([1.0,2.0,3.0],[4.0,5.0,6.0,0.0]),[9.0,0.2,0.3,8.0],[0.5,0.6])
])
def test_invaild_input_numercial(fluxes,r,w):
    with pytest.raises(ValueError):
        discretization.continuity_residual_nd(fluxes,r,w)