#无量纲 Poisson 与独立 SI 教学验证
#python -m pytest
import pytest
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve
from minitcad import discretization
from minitcad.core.constants import q as electron_charge
from minitcad.core.material import SILICON
from minitcad.physics.normalization import build_normalization_scales
from examples.poisson_si_oracle import poisson_residual_si_mms
#内部电荷等于0的直线构造解
@pytest.mark.parametrize(
    "xi",
    [
        np.linspace(0.0,1.0,9),
        np.array([0.0,0.07,0.26,0.61,1.0])
    ],
    ids = ["uniform","nonuniform"]
)
def test_linear_manufactured_solution_has_zero_residual(xi):
    m = xi.size
    u = 0.4 - 0.7 * xi
    p = np.ones(m)
    n = np.ones(m)
    Cbar = np.zeros(m)
    f = discretization.poisson_discretization(xi,u,n,p,Cbar,u[0],u[-1])
    np.testing.assert_allclose(f,0.0,rtol = 0.0,atol = 1e-11)
#内部电荷不为零的均匀网格抛物线解
def test_uniform_quadratic_solution():
    xi = np.linspace(0.0,1.0,11)
    u = 0.2 + 0.3 * xi - 0.6 * xi ** 2
    m = xi.size
    p = np.ones(m)
    n = np.ones(m)
    Cbar = np.full(m,1.2)
    f = discretization.poisson_discretization(xi,u,n,p,Cbar,u[0],u[-1])
    np.testing.assert_allclose(f,0.0,rtol = 0.0,atol = 1e-11)
#内部电荷不为零的非均匀网格抛物线解
@pytest.mark.parametrize("xi",[
    np.array([0.0,0.25,0.75,1.5]),
    np.array([0.0,0.1,0.4,1.0]),
    np.array([0.0,0.01,0.1,0.4,0.9,1.0])
])
def test_nonuniform_quadratic_solution(xi):
    u = 0.2 + 0.3 * xi - 0.6 * xi ** 2
    m = xi.size
    n = np.ones(m)
    p = np.ones(m)
    cbar = np.full(m,1.2)
    f = discretization.poisson_discretization(xi,u,n,p,cbar,u[0],u[-1])
    np.testing.assert_allclose(f,0.0,rtol = 0.0,atol = 1e-11)
#狄拉克边界条件
def test_dirichlet_rows_are_potential_differences():
    xi = np.array([0.0,0.25,0.75,1.5])
    u = 1.0 + 0.3 * xi
    m = xi.size
    p = np.ones(m)
    n = np.ones(m)
    Cbar = np.zeros(m)
    left_target = u[0] + 0.2
    right_target = u[-1] - 0.4
    f = discretization.poisson_discretization(xi,u,n,p,Cbar,left_target,right_target)
    np.testing.assert_allclose(
        f[[0,-1]],
        [-0.2,0.4],
        atol = 1e-12,
        rtol = 0.0
    )
    np.testing.assert_allclose(f[1:-1],0.0,atol = 1e-11,rtol = 0.0)
#电荷符号
@pytest.mark.parametrize("charge",[2.0,-2.0])
def test_charge_sign(charge):
    xi = np.array([0.0,0.3,0.8,1.0])
    m = xi.size
    f = discretization.poisson_discretization(xi,np.zeros(m),np.ones(m),
                                              np.ones(m),np.full(m,charge),0.0,0.0)
    np.testing.assert_allclose(f[1:-1],charge,rtol = 0,atol = 1e-12)
#SI单位与归一化换算测试,vibecoding
def test_si_and_dimensionless_residual_scaling():
    """Check SI/ND scaling for internal and boundary rows."""
    q_C = 1.602176634e-19
    kB_J_K = 1.380649e-23
    eps0_F_m = 8.8541878128e-12
    temperature_K = 300.0
    ni_cm3 = 1.0e10
    epsilon_F_m = 11.7 * eps0_F_m
    vt_V = kB_J_K * temperature_K / q_C
    epsilon_F_cm = epsilon_F_m / 100.0
    l0_cm = np.sqrt(
        epsilon_F_cm * vt_V / (q_C * ni_cm3)
    )
    ni_m3 = ni_cm3 * 1.0e6
    xi = np.array([0.0, 0.25, 0.75, 1.5])
    x_m = xi * l0_cm * 1.0e-2
    u = 0.5 + 0.25 * xi - 0.8 * xi**2
    n = np.array([1.1, 1.3, 2.1, 2.5])
    p = np.array([1.4, 1.1, 0.9, 0.6])
    cbar = np.array([0.0, -0.3, 0.8, 0.0])
    u_left_target = u[0] + 0.05
    u_right_target = u[-1] - 0.07
    f_nd = discretization.poisson_discretization(
        xi, u, n, p, cbar,
        u_left_target, u_right_target
    )
    rho_C_m3 = q_C * ni_m3 * (cbar + p - n)
    f_si = poisson_residual_si_mms(
        x_m,
        vt_V * u,
        rho_C_m3,
        epsilon_F_m,
        vt_V * u_left_target,
        vt_V * u_right_target,
    )
    # Internal rows: F_SI = -q * ni_SI * F_ND
    np.testing.assert_allclose(
        f_si[1:-1] / (-q_C * ni_m3),
        f_nd[1:-1],
        rtol=5e-13,
        atol=5e-13,
    )
    # Boundary rows: F_SI = VT * F_ND
    np.testing.assert_allclose(
        f_si[[0, -1]] / vt_V,
        f_nd[[0, -1]],
        rtol=5e-13,
        atol=5e-13,
    )
#异常网格测试
@pytest.mark.parametrize(
    "bad",
    [
        [0.0,1.0],
        [0.0,0.2,0.2,1.0],
        [0.0,0.7,0.4,1.0],
        [0.0,float("nan"),1.0],
        [[0.0,0.5,1.0]]
    ],
    ids = ["too short",
           "duplicate",
           "descending",
           "NaN",
           "not-1d"]
)
def test_rejects_invalid_mesh(bad):
    m = np.asarray(bad).size
    with pytest.raises(ValueError):
        discretization.poisson_discretization(
            bad,
            np.zeros(m),np.ones(m),np.ones(m),np.zeros(m),0.0,0.0
        )
#异常参数测试
@pytest.mark.parametrize(
    ("field","bad"),
    [
        ("u",[0.0,0.0]),
        ("n",[1.0,0.0,1.0,0.0]),
        ("p",[1.0,-1.0,1.0,1.0]),
        ("u",[1.0,float("inf"),1.0,1.0]),
        ("cbar",[0.0,float("nan"),0.0,0.0]),
        ("cbar",[[0.0,0.0,0.0,00.0]])
    ]
)
def test_reject_invaild_parameters(field,bad):
    args = dict(
        xi = np.array([0.0,0.25,0.75,1.0]),
        u = np.zeros(4),
        n = np.ones(4),
        p = np.ones(4),
        cbar = np.zeros(4),
        u_left = 0.0,
        u_right = 0.0
    )
    args[field] = bad
    with pytest.raises(ValueError):
        discretization.poisson_discretization(**args)
#异常边界条件测试
@pytest.mark.parametrize("side",["u_left","u_right"])
@pytest.mark.parametrize("bad",[float("nan"),float("inf"),-float("inf")])
def test_reject_invaild_boundary_potential(side,bad):
    args = dict(
        xi = np.array([0.0,0.25,0.75,1.0]),
        u = np.zeros(4),
        n = np.ones(4),
        p = np.ones(4),
        cbar = np.zeros(4),
        u_left = 0.0,
        u_right = 0.0
    )
    args[side] = bad
    with pytest.raises(ValueError):
        discretization.poisson_discretization(**args)
#正确形式的输出且不改变原数组
def test_output_shape_dtype_and_input_immutability():
    xi = np.array([0.0, 0.15, 0.6, 1.0])
    u = xi ** 2
    n = np.ones(4)
    p = np.ones(4)
    cbar = np.full(4, -2.0)
    originals = []
    for a in (xi,u,n,p,cbar):
        originals.append(a.copy())
    result = discretization.poisson_discretization(xi,u,n,p,cbar,u[0],u[-1])
    assert isinstance(result,np.ndarray)
    assert result.shape == (4,)
    assert result.dtype == np.float64
    assert np.isfinite(result).all()
    for actual,original in zip(
        (xi,u,n,p,cbar),originals
    ):
        np.testing.assert_array_equal(actual,original)
#正弦解测试函数
sine_lm = 1e-6
sine_epsilon_Fm = 11.7 * 8.8541878128e-12
def solve_sine_si(node_count:int):
    if (not isinstance(node_count,int)
        or node_count < 3):
        raise ValueError(
            "node_count must be integer bigger than 2"
        )
    x_m = np.linspace(0.0,sine_lm,node_count)
    h_m = x_m[1] - x_m[0]
    phi_exact_V = np.sin(
        np.pi * x_m / sine_lm
    )
    rho_exact_m3 = sine_epsilon_Fm * (np.pi / sine_lm) ** 2 * phi_exact_V
    interior_count = node_count - 2
    matrix = diags(
        diagonals = [
            np.full(interior_count - 1,- sine_epsilon_Fm / h_m ** 2),
            np.full(interior_count,2 * sine_epsilon_Fm / h_m ** 2),
            np.full(interior_count - 1,- sine_epsilon_Fm / h_m ** 2)
        ],
        offsets = [1,0,-1],
        format = "csc",
    )
    phi_numercial_V = np.zeros(node_count)
    phi_numercial_V[1:-1] += spsolve(matrix,rho_exact_m3[1:-1])
    return (x_m,phi_numercial_V,phi_exact_V,rho_exact_m3)
#测试函数的验证：节点和边界
@pytest.mark.parametrize(
    "node_count",
    [17,33,65,129]
)
def test_sine_solution_shape_and_boundary(node_count):
    x_m,phi_num,phi_exact,rho = solve_sine_si(node_count)
    assert x_m.shape == (node_count,)
    assert phi_num.shape == (node_count,)
    assert phi_exact.shape == (node_count,)
    assert rho.shape == (node_count,)
    np.testing.assert_allclose(
        phi_num[[0,-1]],
        [0.0,0.0],
        atol = 1e-14,
        rtol = 0.0
    )
    assert np.isfinite(phi_num).all()
#二阶收敛测试
def test_sine_solution_second_order():
    node_counts = [17,33,65,129]
    errors = []
    for N in node_counts:
        x_m,phi_num,phi_exact,rho = solve_sine_si(N)
        err = np.max(np.abs(phi_exact - phi_num))
        assert np.isfinite(err)
        if len(errors) != 0:
            ratio = errors[-1] / err
            order = np.log2(ratio)
            assert 3.8 <= ratio <= 4.2
            assert order == pytest.approx(2.0,abs = 0.003,rel = 0.0)
        errors.append(err)
    for i in range(len(errors) - 1):
        assert errors[i + 1] < errors[i]
#独立残差检查，vibecoding
@pytest.mark.parametrize(
    "node_count",
    [17, 33, 65, 129]
)
def test_sine_si_residual(node_count):
    # 获得数值解
    x_m, phi_num, phi_exact, rho = solve_sine_si(
        node_count
    )
    # 使用原代码中的独立 SI Poisson 验证器
    residual = poisson_residual_si_mms(
        x_m,
        phi_num,
        rho,
        sine_epsilon_Fm,
        0.0,
        0.0
    )
    # 源项尺度
    source_scale = max(
        1.0,
        np.max(np.abs(rho[1:-1]))
    )
    # 内部节点归一化残差
    error = (
        np.max(np.abs(residual[1:-1]))
        / source_scale
    )
    # 残差必须足够接近 0
    assert error < 1e-11
    # 两端 Dirichlet 残差必须为 0
    np.testing.assert_allclose(
        residual[[0, -1]],
        [0.0, 0.0],
        atol=1e-14,
        rtol=0.0
    )
#正弦解测试
@pytest.mark.parametrize(
    "node_count",
    [17,33,65,129]
)
def test_sine_solution(node_count):
    x_m, phi_num ,_, rho = solve_sine_si(
        node_count
    )
    scales = build_normalization_scales(
        SILICON,300,reference_diffusivity_cm2_s=10.0
    )
    l0_cm = scales.l0_cm
    v0_v = scales.v0_v
    n0_per_cm3 = scales.n0_cm3
    x_m_normalized = x_m * 1e2 / l0_cm 
    phi_num_normalized = phi_num / v0_v
    cbar = rho * 1e-6 / (n0_per_cm3 * electron_charge)
    source_scale = max(1,np.max(np.abs(cbar)))
    f = discretization.poisson_discretization(x_m_normalized,
                                              phi_num_normalized,
                                              np.ones(node_count),
                                              np.ones(node_count),
                                              cbar,
                                              phi_num_normalized[0],
                                              phi_num_normalized[-1])
    np.testing.assert_allclose(f[[0,-1]],0.0,atol = 1e-11,rtol = 0.0)
    relative_residual = np.max(np.abs(f[1:-1])) / source_scale
    assert relative_residual < 1e-10