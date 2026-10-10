# -*- coding: utf-8 -*-
"""
pn_junction.py 的测试
=====================

运行方式（不需要装 pytest）：

    python test_pn_junction.py

装了 pytest 的话也可以直接 `pytest test_pn_junction.py`，因为这些函数就是
普通的 test_* 函数。

这些测试的原则：**尽量和解析公式比，而不是和程序自己的输出比**。
程序画出"像二极管"的曲线并不说明它对；下面这些检查对的是教科书上的解析
结果（内建电势、质量作用定律、耗尽近似峰值电场、短二极管饱和电流公式）。
"""

import numpy as np

import pn_junction as pn          # 被测模块


# ============================================================================
# 小工具：把费时的求解结果缓存起来，多个测试共用
# ============================================================================

_CACHE = {}


def equilibrium():
    """平衡解，只算一次。"""
    if "eq" not in _CACHE:
        _CACHE["eq"] = pn.solve_equilibrium()
    return _CACHE["eq"]


def forward_sweep(v_max, step):
    """正偏扫描，只算一次（同一个参数）。"""
    key = ("sweep", v_max, step)
    if key not in _CACHE:
        _CACHE[key] = pn.run_sweep(v_max=v_max, step=step, verbose=False)
    return _CACHE[key]


def ideality_factors(voltages, currents):
    """相邻两点之间的理想因子 n = dV / (VT*ln(J2/J1))。

    从第 3 个点开始算：0 -> 第一个偏压那一段的 J(0) 是零，比值没有意义。
    """
    values = []
    for i in range(2, len(voltages)):
        ratio = currents[i] / currents[i - 1]
        values.append((voltages[i] - voltages[i - 1]) / (pn.V0 * np.log(ratio)))
    return values


# ============================================================================
# 1. 基础零件
# ============================================================================

def test_bernoulli_at_zero_is_one():
    """B(0) 取极限值 1（如果用朴素写法这里会算出 0/0 = nan）。"""
    assert pn.bernoulli(0.0) == 1.0
    assert pn.bernoulli(np.array([0.0, 0.0]))[0] == 1.0


def test_bernoulli_is_accurate_for_small_arguments():
    """小参数下必须准：B(z) ≈ 1 - z/2 + z^2/12。

    这一段专门盯住"用 exp(z)-1 还是 expm1(z)"这个坑：
    朴素写法在 z=1e-9 时只剩下 7 位有效数字。
    """
    z = 1.0e-9
    naive = z / (np.exp(z) - 1.0)            # 错误示范
    series = 1.0 - z / 2.0 + z * z / 12.0    # 解析展开
    assert abs(naive - series) > 1e-8        # 朴素写法确实掉了精度
    assert abs(pn.bernoulli(z) - series) < 1e-16   # 我们的写法没问题
    # 中等大小的参数也检查一下
    for value in (0.01, 0.5, 3.0, -2.0):
        expected = value / np.expm1(value)
        assert abs(pn.bernoulli(value) - expected) < 1e-15 * abs(expected) + 1e-15


def test_bernoulli_identity_used_by_sg_flux():
    """SG 通量变形的关键恒等式：exp(z)*B(z) = B(-z)。"""
    z = np.linspace(-50.0, 50.0, 201)
    assert np.allclose(np.exp(z) * pn.bernoulli(z), pn.bernoulli(-z), rtol=1e-11)


def test_mesh_and_control_volumes():
    """网格几何：节点递增、控制体宽度加起来等于器件长度。"""
    assert pn.N_NODES > 100
    assert np.all(np.diff(pn.NODE_X) > 0)
    assert abs(np.sum(pn.VOL) - pn.LENGTH) < 1e-18
    # 结附近应该比两侧细
    assert pn.H.min() < pn.H.max() / 3.0


def test_contact_potential_is_half_the_built_in_potential():
    """欧姆接触的电势 = ∓Vbi/2，Vbi 用解析公式核对。"""
    v_bi = pn.V0 * np.log(pn.ACCEPTOR * pn.DONOR / pn.NI ** 2)
    # P 区（受主）：掺杂是负的
    phi_p = pn.contact_phi(-pn.ACCEPTOR, 0.0)
    phi_n = pn.contact_phi(+pn.DONOR, 0.0)
    assert abs(phi_p - (-v_bi / 2)) < 1e-12
    assert abs(phi_n - (+v_bi / 2)) < 1e-12
    # 加偏压只是整体平移
    assert abs(pn.contact_phi(+pn.DONOR, 0.3) - (phi_n + 0.3)) < 1e-12


# ============================================================================
# 2. 初值：还没求解就应该满足两条物理性质
# ============================================================================

def test_initial_guess_has_zero_current_and_zero_continuity_residual():
    """电荷中性初值的准费米势是平的，所以电流精确为零。

    检查方式：连续性方程对应的那两段残差应该**恰好**是 0。
    这也是本文件里那个"expm1 变形"的必要性所在——如果用教科书形式，
    这里会得到一堆 ~1e-8 的舍入噪声，而不是 0。
    """
    X = pn.initial_guess(0.0, 0.0)
    n = pn.N_NODES
    F = pn.residual(X, 0.0, 0.0)
    assert np.all(F[n:2 * n] == 0.0)         # 电子连续性残差恰为 0
    assert np.all(F[2 * n:3 * n] == 0.0)     # 空穴连续性残差恰为 0
    assert abs(pn.terminal_current(X)) < 1e-25
    # 泊松残差只在结附近的过渡区不为零（两侧中性区已经满足）
    poisson = np.abs(F[0:n])
    assert np.sum(poisson > 1e-3 * poisson.max()) < n // 10


def test_initial_guess_satisfies_the_boundary_conditions():
    """初值在两端节点上必须正好等于边界条件。"""
    X = pn.initial_guess(0.0, 0.0)
    n = pn.N_NODES
    assert abs(X[pn.LEFT_INDEX] * pn.V0 - pn.contact_phi(pn.DOPING[pn.LEFT_INDEX], 0.0)) < 1e-12
    assert abs(X[pn.RIGHT_INDEX] * pn.V0 - pn.contact_phi(pn.DOPING[pn.RIGHT_INDEX], 0.0)) < 1e-12
    assert X[n + pn.LEFT_INDEX] == 0.0 and X[2 * n + pn.RIGHT_INDEX] == 0.0


# ============================================================================
# 3. 平衡解：对接解析结果
# ============================================================================

def test_equilibrium_converges():
    X, info = equilibrium()
    assert info["converged"], info.get("message", "")
    assert info["iterations"] <= 15


def test_equilibrium_obeys_mass_action_everywhere():
    """平衡时两个准费米势相等，所以 n*p = ni^2 逐点成立（机器精度）。"""
    X, info = equilibrium()
    assert pn.mass_action_error(X) < 1e-12


def test_equilibrium_carries_no_current():
    """平衡态电流为零，而且电流剖面在器件内恒定。"""
    X, info = equilibrium()
    assert abs(pn.terminal_current(X)) < 1e-20
    assert pn.current_spread(X) < 1e-18


def test_equilibrium_contact_potentials_and_peak_field():
    """接触电势 = ∓Vbi/2，峰值电场与耗尽近似一致（允许 10% 偏差）。

    数值上比解析值低约 4%，原因是耗尽近似假设耗尽区边界是突变的，
    而真实解在过渡区还有可动载流子分担一部分电荷。
    """
    X, info = equilibrium()
    phi = X[0:pn.N_NODES] * pn.V0
    v_bi = pn.V0 * np.log(pn.ACCEPTOR * pn.DONOR / pn.NI ** 2)
    assert abs(phi[0] - (-v_bi / 2)) < 1e-6
    assert abs(phi[-1] - (+v_bi / 2)) < 1e-6

    e_face = -np.diff(phi) / pn.H
    e_max = float(np.max(np.abs(e_face)))
    eps = pn.EPS_R * pn.EPS_VACUUM
    e_analytic = float(np.sqrt(pn.Q_CHARGE * pn.ACCEPTOR * v_bi / eps))
    ratio = e_max / e_analytic
    assert 0.90 < ratio < 1.02, "峰值电场 %.4e vs 解析 %.4e" % (e_max, e_analytic)

    # 电场应该在两端为零、峰值出现在结上
    assert abs(e_face[0]) < 1e-6 * e_max
    assert abs(e_face[-1]) < 1e-6 * e_max
    peak = pn.NODE_X[int(np.argmax(np.abs(e_face)))]
    assert abs(peak - pn.JUNCTION) < 2 * pn.H.min()


# ============================================================================
# 4. 正偏：I-V 曲线
# ============================================================================

def test_forward_current_increases_monotonically():
    voltages, currents, results = forward_sweep(0.5, 0.05)
    assert abs(currents[0]) < 1e-20                # 0 V 时电流为零
    assert np.all(currents[1:] > 0.0)              # 正偏电流为正
    assert np.all(np.diff(currents) > 0.0)         # 单调上升
    assert all(info["converged"] for (v, j, X, info) in results)


def test_ideality_factor_is_one_in_the_diffusion_range():
    """扩散电流占主导的区间，理想因子必须接近 1。

    这一条同时验证了扩散系数、几何尺寸和边界条件——因为 J ∝ exp(V/VT)
    这个斜率只由热电压决定，偏一点就会看出来。
    """
    voltages, currents, results = forward_sweep(0.5, 0.05)
    n_ideal = ideality_factors(voltages, currents)
    # n_ideal[0] 对应 0.05-0.10 V ... n_ideal[8] 对应 0.45-0.50 V
    high = n_ideal[6:9]                          # 0.35 V 以上
    assert all(1.0 <= value <= 1.03 for value in high), high
    # 低偏压时复合电流占主导，理想因子应该明显大于 1
    assert n_ideal[0] > 1.15, n_ideal[0]


def test_saturation_current_matches_the_short_diode_formula():
    """把寿命调长（复合可忽略）之后，正偏电流必须等于短二极管解析式：

        J = q*ni^2*( Dp/(Wn*ND) + Dn/(Wp*NA) ) * (exp(V/VT) - 1)

    这是整个测试里最强的一条：它把扩散系数、几何、边界条件钉在一起。
    """
    saved = (pn.TAU_N, pn.TAU_P, pn.TAU_N_T, pn.TAU_P_T)
    try:
        # 临时的长寿命（复合项几乎不起作用）
        pn.TAU_N = pn.TAU_P = 1.0e-3
        pn.TAU_N_T = pn.TAU_P_T = 1.0e-3 / pn.TAU0
        voltages, currents, results = pn.run_sweep(v_max=0.4, step=0.1, verbose=False)
    finally:
        pn.TAU_N, pn.TAU_P, pn.TAU_N_T, pn.TAU_P_T = saved

    w_p = pn.JUNCTION
    w_n = pn.LENGTH - pn.JUNCTION
    j_s_analytic = pn.Q_CHARGE * pn.NI ** 2 * (
        pn.MU_P * pn.V0 / (w_n * pn.DONOR) + pn.MU_N * pn.V0 / (w_p * pn.ACCEPTOR)
    )
    # 每个偏压反推出来的饱和电流都应该一样（说明斜率是 exp(V/VT)）
    for v, j in zip(voltages[1:], currents[1:]):
        j_s = j / np.expm1(v / pn.V0)
        assert 0.9 < j_s / j_s_analytic < 1.1, \
            "%.1f V: 反推 Js=%.4e, 解析 %.4e" % (v, j_s, j_s_analytic)


def test_current_profile_is_conserved_and_improves_with_current():
    """稳态下总电流沿器件恒定；电流越大，相对误差越小。"""
    voltages, currents, results = forward_sweep(0.5, 0.05)
    sample = []
    for (v, j, X, info) in results:
        if v < 0.15:                      # 电流太小时谈不上守恒
            continue
        spread = pn.current_spread(X) / abs(pn.total_current_profile(X)[0])
        sample.append((v, spread))
    assert len(sample) >= 5, sample
    assert all(value < 1e-2 for (v, value) in sample), sample
    # 高偏压（大电流）应该比低偏压准
    assert sample[-1][1] < sample[0][1], sample


def test_warm_start_preserves_the_carrier_densities():
    """延续法平移的是三个势，载流子浓度必须保持不变。"""
    X, info = equilibrium()
    X_shifted = pn.warm_start(X, 0.1)
    n = pn.N_NODES
    before = pn.exp_capped(X[0:n] - X[n:2 * n])
    after = pn.exp_capped(X_shifted[0:n] - X_shifted[n:2 * n])
    assert np.allclose(before, after, rtol=1e-12)
    # 两端边界条件也跟着平移了
    assert abs((X_shifted[n + pn.LEFT_INDEX] - X[n + pn.LEFT_INDEX]) * pn.V0) < 1e-12
    assert abs((X_shifted[n + pn.RIGHT_INDEX] - X[n + pn.RIGHT_INDEX]) * pn.V0 - 0.1) < 1e-12


# ============================================================================
# 5. 直接运行本文件时的小测试运行器
# ============================================================================

if __name__ == "__main__":
    import sys
    import traceback

    tests = [(name, func) for name, func in sorted(globals().items())
             if name.startswith("test_") and callable(func)]
    failed = 0
    for name, func in tests:
        try:
            func()
        except Exception:
            failed += 1
            print("失败  %s" % name)
            traceback.print_exc()
        else:
            print("通过  %s" % name)
    print()
    print("%d 个测试, %d 个失败" % (len(tests), failed))
    sys.exit(1 if failed else 0)
