# -*- coding: utf-8 -*-
"""
一维 PN 结漂移扩散模拟器（精简示例版）
======================================

这一个文件就完成了全部流程：

    建网格 -> 生成掺杂 -> 求平衡解 -> 逐步加偏压 -> 得到 I-V 曲线

只依赖 numpy 和 scipy（稀疏线性方程组交给 scipy 解）。


一、要解什么
------------
一维稳态漂移扩散方程组，三个方程耦合在一起：

    泊松方程:     -d/dx( eps * dphi/dx ) = q * (p - n + ND - NA)
    电子连续性:   dJn/dx =  q * R
    空穴连续性:   dJp/dx = -q * R

其中 R 是 SRH 复合率（电子和空穴通过缺陷能级复合）。
三个方程互相纠缠：电势决定载流子怎么流动，载流子又反过来决定电荷，
所以整体是一个非线性问题，只能迭代求解。


二、未知量为什么选"准费米势"，而不是 (phi, n, p)
-------------------------------------------------
每个节点上有三个未知量，但这里用的是三个势（单位都是伏特）：

    phi     静电势
    phi_n   电子准费米势      n = ni * exp( (phi   - phi_n) / VT )
    phi_p   空穴准费米势      p = ni * exp( (phi_p - phi  ) / VT )

好处有三个，后面写代码时会反复用到：
  1. 浓度永远由指数算出，一定是正的，Newton 迭代不会跑出负浓度；
  2. 三个未知量都在 0~30 这个量级（不像浓度跨 10^6 倍），矩阵条件数好；
  3. "平衡态"就是准费米势处处相等（phi_n = phi_p = 常数），所以初值和
     边界条件都特别好写。


三、怎么离散
------------
电流用 Scharfetter-Gummel (SG) 格式，这是半导体数值计算里的标准做法：
它把"漂移"和"扩散"两种机制合在一个式子里处理，在电场很强、浓度变化很陡
的地方也不会算出非物理的振荡。

数学上，教科书里 SG 通量写成：

    Jn[i+1/2] = (q*Dn/h) * [ n[i+1]*B(dphi) - n[i]*B(-dphi) ]
    B(z) = z / (exp(z) - 1)        （Bernoulli 函数）

但这份代码里用的是和它**数学上完全等价**的另一种写法（见 sg_flux 函数的
注释），原因写在那个函数里——这是整个文件里最值得注意的一处细节。


四、怎么求解
------------
把所有离散方程堆成一个向量 F(X) = 0（X 是全部未知量），用阻尼 Newton 法：

    1. 算残差 F(X)
    2. 数值差分算出 Jacobian 矩阵 J = dF/dX
    3. 解线性方程组 J*dX = -F，得到修正量 dX
    4. X = X + alpha*dX（alpha 由回溯线搜索确定，保证残差真的变小）
    5. 残差足够小就停，否则回到第 1 步

偏压扫描时，每一个已经收敛的解都作为下一个偏压的初值（这叫"延续法"），
因为直接从一个解跳到相差很远的另一个解，Newton 很容易发散。
"""

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve


# ============================================================================
# 1. 物理常量
# ============================================================================
# 全部使用"半导体常用单位"：长度 cm，电压 V，浓度 cm^-3，
# 迁移率 cm^2/(V*s)，介电常数 F/cm，电流密度 A/cm^2。

Q_CHARGE = 1.602176634e-19      # 电子电荷 [C]
K_BOLTZ = 1.380649e-23          # 玻尔兹曼常数 [J/K]
EPS_VACUUM = 8.8541878128e-14   # 真空介电常数 [F/cm]


# ============================================================================
# 2. 器件参数（想改器件就改这一段）
# ============================================================================

TEMPERATURE = 300.0             # 温度 [K]
LENGTH = 1.0e-3                 # 器件总长度 1e-3 cm = 10 um
JUNCTION = 0.5e-3               # PN 结的位置（器件正中间）
ACCEPTOR = 1.0e16               # P 区受主浓度 NA [cm^-3]
DONOR = 1.0e16                  # N 区施主浓度 ND [cm^-3]

# 硅的材料参数（这一版只支持单一材料，且都是常数）
NI = 1.0e10                     # 本征载流子浓度 [cm^-3]（300 K 的硅）
EPS_R = 11.7                    # 相对介电常数
MU_N = 1350.0                   # 电子迁移率 [cm^2/(V*s)]
MU_P = 480.0                    # 空穴迁移率 [cm^2/(V*s)]
TAU_N = 1.0e-6                  # 电子 SRH 寿命 [s]
TAU_P = 1.0e-6                  # 空穴 SRH 寿命 [s]

# 偏压扫描设置：从 0 V 扫到 0.8 V，每步 0.05 V。
# 注意方向：PN 结"正偏"是 P 区加正电压。这个器件 P 区在左边，
# 所以扫左接触；如果给右边（N 区）加正电压，那其实是反偏。
SWEEP_MAX = 0.8                 # 扫到多少伏
SWEEP_STEP = 0.05               # 步长 [V]


# ============================================================================
# 3. 归一化（把物理量换成无量纲量）
# ============================================================================
# 为什么需要：浓度是 1e16 量级、电势只有 0.1 V 量级、长度是 1e-4 量级，
# 直接把这些数扔进 Newton 迭代，矩阵的条件数会非常差。归一化之后所有量都
# 接近 1，迭代稳定得多。
#
# 规则（带波浪号的表示归一化量）：x~ = x/L0,  phi~ = phi/V0,  n~ = n/n0 ...

V0 = K_BOLTZ * TEMPERATURE / Q_CHARGE       # 热电压 VT = kT/q ≈ 0.02585 V
N0 = NI                                     # 浓度以本征浓度为单位
L0 = LENGTH                                 # 长度以器件长度为单位, x~ 在 [0,1]
MU0 = MU_N                                  # 迁移率以电子迁移率为单位

# 由上面几个基本尺度派生出来的量，用来把方程写成无量纲形式
J0 = Q_CHARGE * MU0 * V0 * N0 / L0          # 电流密度尺度 [A/cm^2]
TAU0 = L0 * L0 / (MU0 * V0)                 # 时间尺度 [s]（复合率尺度就是 n0/tau0）

# 泊松方程里空间电荷项的系数，它等于 (本征德拜长度 / 器件长度)^2。
# 这个数看着奇怪（掺杂归一化后是 1e6，乘上 1/LAMBDA2 ≈ 0.06 才和二阶导数
# 项量级相当），但这是这套标度本身的属性，不是写错了。
LAMBDA2 = EPS_R * EPS_VACUUM * V0 / (Q_CHARGE * N0 * L0 * L0)


# ============================================================================
# 4. 网格与掺杂
# ============================================================================

def build_mesh():
    """建网格：结附近加密，因为那里的电势和浓度变化最剧烈。

    分三段拼起来（左粗 / 结附近细 / 右粗），段与段之间共用的那个点
    只保留一次，所以最终节点数是 50 + 100 + 50 - 2 = 198。
    """
    x_coarse = 1.0e-4                        # 结两侧各 1 um 的加密窗口
    left = np.linspace(0.0, JUNCTION - x_coarse, 50)
    middle = np.linspace(JUNCTION - x_coarse, JUNCTION + x_coarse, 100)[1:]
    right = np.linspace(JUNCTION + x_coarse, LENGTH, 50)[1:]
    return np.concatenate([left, middle, right])


# 节点坐标。注意变量名没叫 X，因为 X 在下面专门表示"未知量向量"。
NODE_X = build_mesh()                # 节点坐标 [cm]
N_NODES = NODE_X.size                # 节点数

# 有限体积法要用到的几何量：
#   H   : 相邻节点之间的距离（"面"在两点中间），长度 N-1
#   VOL : 每个节点代表的控制体宽度，长度 N
H = np.diff(NODE_X)
VOL = np.empty(N_NODES)
VOL[0] = 0.5 * H[0]
VOL[-1] = 0.5 * H[-1]
VOL[1:-1] = 0.5 * (H[:-1] + H[1:])

# 归一化后的坐标、间距、控制体宽度
XT = NODE_X / L0
HT = H / L0
VOLT = VOL / L0

# 掺杂：在节点上取值。x < JUNCTION 是 P 区（受主），x >= JUNCTION 是 N 区。
# 净掺杂 DOPING = ND - NA 就是泊松方程右端项里的那一项。
NA_ARR = np.where(NODE_X < JUNCTION, ACCEPTOR, 0.0)
ND_ARR = np.where(NODE_X >= JUNCTION, DONOR, 0.0)
DOPING = ND_ARR - NA_ARR                    # [cm^-3]

# 归一化后的掺杂、迁移率、寿命
DOPING_T = DOPING / N0
MU_N_T = MU_N / MU0
MU_P_T = MU_P / MU0
TAU_N_T = TAU_N / TAU0                      # 注意：这是"归一化"寿命，不是秒
TAU_P_T = TAU_P / TAU0

# 两端的理想欧姆接触（节点下标）
LEFT_INDEX = 0
RIGHT_INDEX = N_NODES - 1


# ============================================================================
# 5. 两个基础函数：贝努利函数和带范围的指数
# ============================================================================

def bernoulli(z):
    """B(z) = z / (exp(z) - 1)，z = 0 处取极限值 1。

    写法上有一个坑：
        z / (np.exp(z) - 1.0)          <-- 不要这样写
    当 z 很小时（比如 1e-9），np.exp(z) 算出来是 1.000000001，减 1 之后
    前面的有效数字全被吃掉了，结果只剩 7 位精度。
    np.expm1(z) 是专门用来算 exp(z)-1 且不丢精度的函数，所以下面用它。
    另外 z -> +无穷 时它自然趋于 0；z -> -无穷 时 expm1(z) -> -1，
    结果趋于 -z，两个方向都不会溢出。
    """
    z = np.asarray(z, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        result = z / np.expm1(z)
    return np.where(z == 0.0, 1.0, result)


def exp_capped(z):
    """带范围的指数函数，防止迭代过程中指数溢出成 inf。

    正常状态下指数不会超过 ±50，这里限制到 ±300 只是为了在 Newton 偶尔
    迈出一大步时不至于算出 inf 或 nan（线搜索会拒绝这种步长）。
    """
    return np.exp(np.clip(z, -300.0, 300.0))


# ============================================================================
# 6. 边界条件（理想欧姆接触）
# ============================================================================
# 物理假设：接触处满足"电中性 + 平衡关系"，并且准费米势等于外加电压：
#
#     p - n = NA - ND ,   n*p = ni^2 ,   phi_n = phi_p = V_applied
#
# 把 n、p 用 phi 的表达式代进去，就能解出接触处的静电势：
#
#     phi = V_applied + VT * asinh( (ND - NA) / (2*ni) )
#
# 这个式子自动给出正确的内建电势：P 区接触是 -Vbi/2，N 区接触是 +Vbi/2。

def contact_phi(doping, v_applied):
    """欧姆接触处的静电势 [V]。doping = ND - NA [cm^-3]。"""
    z = doping / (2.0 * NI)
    return v_applied + V0 * np.arcsinh(z)


# ============================================================================
# 7. 残差函数 F(X)：整个程序的心脏
# ============================================================================
# 未知量向量 X 的排列方式（每段长度都是节点数 N）：
#
#     X = [ phi(0..N-1) | phi_n(0..N-1) | phi_p(0..N-1) ]
#
# F(X) 同样分三段：泊松残差、电子连续性残差、空穴连续性残差。
# 只要这个函数写对了，后面的 Newton、扫压都是通用套路；这个函数写错了，
# Newton 只会"高精度地求出错误方程的解"。

def residual(X, v_left, v_right):
    """给定归一化未知量 X 和两端电压，返回残差向量 F(X)。"""
    n = N_NODES
    phi = X[0:n]                 # 静电势（归一化，单位是 VT）
    phi_n = X[n:2 * n]           # 电子准费米势
    phi_p = X[2 * n:3 * n]       # 空穴准费米势

    # 由三个势算出载流子浓度（归一化，单位是 ni）。
    # 指数里的参数是两个势之差，只跟物理有关，和用什么单位无关。
    dens_n = exp_capped(phi - phi_n)
    dens_p = exp_capped(phi_p - phi)

    F = np.zeros(3 * n)
    interior = slice(1, n - 1)   # 内部节点（两端接触点由边界条件给出）

    # ------------------------------------------------------------------
    # 7.1 泊松方程
    # ------------------------------------------------------------------
    # 物理形式:      -d/dx( eps * dphi/dx ) = q * (p - n + ND - NA)
    # 有限体积形式:  相邻两个面上的电场之差 = 控制体内的净电荷
    #
    # 先算面上的电场（E = -dphi/dx，归一化后单位是 V0/L0）：
    e_face = -np.diff(phi) / HT

    charge = (dens_p - dens_n + DOPING_T)[interior]      # 净电荷（单位 q*ni）
    F[interior] = (e_face[1:] - e_face[:-1]) / VOLT[interior] - charge / LAMBDA2

    # ------------------------------------------------------------------
    # 7.2 连续性方程（电子和空穴）
    # ------------------------------------------------------------------
    # 物理形式:  dJn/dx = q*R  ,  dJp/dx = -q*R
    # 面通量用 SG 格式算，然后做相邻面之差。
    jn_face, jp_face = sg_flux(phi, phi_n, phi_p)

    # SRH 复合率（归一化形式）：
    #     R = (n*p/ni^2 - 1) / ( tau_p*(n/ni + 1) + tau_n*(p/ni + 1) ) / tau0
    # 这个式子对 (n, p, ni) 是一次齐次的、对 tau 是 -1 次齐次的，所以只要
    # 把浓度换成"以 ni 为单位"的归一化浓度（此时 ni 变成 1）、把寿命换成
    # "以 tau0 为单位"的归一化寿命，算出来的就直接是 R/R0。
    #
    # 这里有一个非常容易踩的坑：一边用物理浓度、一边用归一化寿命的话，
    # 复合率会差好几个数量级，而且表面上看起来还挺正常——因为复合项被
    # 悄悄关掉了，前向电流反而"更接近"理想二极管公式。
    #
    # 另外注意 tau 的配对：tau_p 乘的是电子浓度。在强 N 型材料里
    # (n >> ni >> p) 这给出 R ≈ p/tau_p，也就是说少数载流子空穴用它自己的
    # 寿命复合；如果两个 tau 写反了，这个物理性质就错了。
    #
    # 分子里的 n*p/ni^2 直接用准费米势之差来算（不写成 dens_n*dens_p）：
    #     n*p/ni^2 = exp(phi_p - phi_n)
    # 两种写法在数学上等价，但平衡时 phi_p 和 phi_n 是同一个数，这个写法
    # 给出精确的 1（于是 R 精确为 0），而 exp(a)*exp(-a) 在浮点上不一定
    # 等于 1，会留下一个 ~1e-16 的假复合率。
    driving_force = exp_capped(phi_p - phi_n) - 1.0
    rate = driving_force / (
        TAU_P_T * (dens_n + 1.0) + TAU_N_T * (dens_p + 1.0)
    )

    F[n + 1:2 * n - 1] = (jn_face[1:] - jn_face[:-1]) / VOLT[interior] - rate[interior]
    F[2 * n + 1:3 * n - 1] = (jp_face[1:] - jp_face[:-1]) / VOLT[interior] + rate[interior]

    # ------------------------------------------------------------------
    # 7.3 边界条件：两端接触点用 Dirichlet 条件
    # ------------------------------------------------------------------
    # 做法是把那一行的方程替换成 "X - 边界值 = 0"。
    # 接触处 phi_n = phi_p = 外加电压，phi 由电中性关系算出。
    F[LEFT_INDEX] = phi[LEFT_INDEX] - contact_phi(DOPING[LEFT_INDEX], v_left) / V0
    F[n + LEFT_INDEX] = phi_n[LEFT_INDEX] - v_left / V0
    F[2 * n + LEFT_INDEX] = phi_p[LEFT_INDEX] - v_left / V0

    F[RIGHT_INDEX] = phi[RIGHT_INDEX] - contact_phi(DOPING[RIGHT_INDEX], v_right) / V0
    F[2 * n - 1] = phi_n[RIGHT_INDEX] - v_right / V0
    F[3 * n - 1] = phi_p[RIGHT_INDEX] - v_right / V0

    return F


def sg_flux(phi, phi_n, phi_p):
    """算每个面上的电子、空穴电流（归一化）。

    教科书里的 SG 格式是：

        Jn[i+1/2] = (mu_n/h) * [ n[i+1]*B(dphi) - n[i]*B(-dphi) ]

    数学上它和下面这种写法完全等价（用恒等式 exp(x)*B(x) = B(-x)，
    再把 exp(a)-exp(b) 改写成 exp(b)*expm1(a-b)）：

        Jn[i+1/2] =  (mu_n/h) * n[i] * B(-dphi) * expm1(-dphi_n)
        Jp[i+1/2] = -(mu_p/h) * p[i] * B( dphi) * expm1( dphi_p)

    两种写法算出来的数在数学上一样，但在计算机上差别极大：
    平衡态下准费米势是"平的"，相邻节点之差只有 1e-15 量级。教科书那种
    写法是拿两个 1e8 量级的数相减，浮点舍入误差(~1e-8)比要算的电流还大，
    结果是一条纯粹由舍入噪声构成的"假电流"。
    下面这种写法在准费米势相等时给出**精确的 0**，而且能一眼看出物理意义：
    电流由准费米势的梯度驱动，Jn = -q*mu_n*n*grad(phi_n)。
    """
    dphi = np.diff(phi)                     # 面上两侧的静电势差
    dphi_n = np.diff(phi_n)                 # 面上两侧的电子准费米势差
    dphi_p = np.diff(phi_p)                 # 面上两侧的空穴准费米势差

    # 面左侧节点的浓度（归一化）
    n_left = exp_capped(phi[:-1] - phi_n[:-1])
    p_left = exp_capped(phi_p[:-1] - phi[:-1])

    jn = MU_N_T / HT * n_left * bernoulli(-dphi) * np.expm1(-dphi_n)
    # 空穴电流前面有负号：Jn = -q*mu_n*n*grad(phi_n)，Jp = -q*mu_p*p*grad(phi_p)，
    # 两者形式一样，但把 p 的指数关系代进去之后会多出一个负号。
    jp = -MU_P_T / HT * p_left * bernoulli(dphi) * np.expm1(dphi_p)
    return jn, jp


# ============================================================================
# 8. 初值：电荷中性猜测
# ============================================================================

def initial_guess(v_left, v_right):
    """构造一个电荷中性、电流为零的初值。

    思路：
      * 准费米势取"两端电压的线性插值"，这样两端正好满足边界条件；
      * 静电势按"每个节点都电中性"来定：phi = v(x) + VT*asinh(掺杂/(2ni))。
    这个初值在两侧中性区已经精确满足泊松方程，只有结附近的过渡区需要
    Newton 去修正；而且因为准费米势是平的，电流恰好为零。
    """
    v_applied = v_left + (v_right - v_left) * XT     # 外加电压的线性插值
    z = DOPING_T / 2.0
    phi = v_applied + np.arcsinh(z)                  # 归一化静电势
    return np.concatenate([phi, v_applied, v_applied])


# ============================================================================
# 9. 数值差分 Jacobian
# ============================================================================

def fd_jacobian(f, X):
    """用前向差分算 Jacobian 矩阵 J[i,j] = dF[i]/dX[j]。

    每个未知量扰动一次、算一遍残差，所以 3N 个未知量要算 3N 次残差。
    好处是永远不会和残差公式不一致（不需要另外维护一套解析导数）。
    这个矩阵真实的非零结构是块三对角的（一个节点只和左右邻居耦合），
    正式实现应该用解析导数 + 稀疏装配；这里为了代码简短，先组装成稠密
    矩阵再转成稀疏矩阵交给 scipy。
    """
    m = X.size
    f0 = f(X)
    J = np.empty((m, m))
    X_work = X.copy()
    for j in range(m):
        # 步长按变量大小缩放：X 里的量都是 O(1) 到 O(30)
        step = 1.0e-7 * max(abs(X[j]), 1.0)
        old = X_work[j]
        X_work[j] = old + step
        J[:, j] = (f(X_work) - f0) / step
        X_work[j] = old
    return csr_matrix(J)


# ============================================================================
# 10. 阻尼 Newton 法
# ============================================================================

def newton_solve(v_left, v_right, X_start, max_iterations=50, tol=1.0e-10):
    """求解 F(X) = 0，返回 (X, 信息字典)。"""
    X = X_start.astype(float).copy()
    F = residual(X, v_left, v_right)
    norm = float(np.max(np.abs(F)))
    norm0 = norm                      # 初始残差，用来做"相对收敛判据"
    history = [norm]

    # 收敛判据：残差降到初始值的 tol 倍以下。
    # 为什么不要求残差降到某个绝对小量？因为中性区的电流是两个大通量相减
    # 得到的，残差存在一个"浮点舍入地板"（连续性问题大约 1e-4 量级），
    # 只有相对判据才既有意义又不会永远不收敛。
    threshold = tol * norm0

    # 单步电位更新的上限（单位 VT）。残差对电位差是指数敏感的，整步 Newton
    # 从冷启动出发可能一步跨出去好几个数量级，光靠线搜索救不回来；
    # 把步长等比缩小（只缩放、不截断，保持下降方向）就稳了。
    max_update = 5.0

    for iteration in range(max_iterations):
        if norm <= threshold:
            return X, {"converged": True, "iterations": iteration,
                       "residual": norm, "history": history}

        # 解线性方程组 J*dX = -F
        J = fd_jacobian(lambda Z: residual(Z, v_left, v_right), X)
        dX = spsolve(J, -F)
        if not np.all(np.isfinite(dX)):
            return X, {"converged": False, "iterations": iteration,
                       "residual": norm, "history": history,
                       "message": "Jacobian 奇异"}

        # 步长限制：三个势场里移动最大的那个不超过 max_update
        biggest = 0.0
        for block in range(3):
            block_max = np.max(np.abs(dX[block * N_NODES:(block + 1) * N_NODES]))
            biggest = max(biggest, float(block_max))
        alpha = min(1.0, max_update / biggest) if biggest > 0.0 else 1.0

        # 回溯线搜索：只接受让残差变小的步长，最多折半 20 次
        accepted = False
        for _ in range(20):
            X_try = X + alpha * dX
            F_try = residual(X_try, v_left, v_right)
            norm_try = float(np.max(np.abs(F_try)))
            if np.isfinite(norm_try) and norm_try < norm:
                accepted = True
                break
            alpha *= 0.5
        if not accepted:
            return X, {"converged": False, "iterations": iteration,
                       "residual": norm, "history": history,
                       "message": "找不到能让残差下降的步长"}

        X, F, norm = X_try, F_try, norm_try
        history.append(norm)

    return X, {"converged": False, "iterations": max_iterations,
               "residual": norm, "history": history,
               "message": "达到最大迭代次数"}


# ============================================================================
# 11. 后处理：电流、以及"结果对不对"的检查
# ============================================================================

def total_current_profile(X):
    """每个面上的总电流密度（归一化）。稳态下它应该是常数。"""
    n = N_NODES
    jn, jp = sg_flux(X[0:n], X[n:2 * n], X[2 * n:3 * n])
    return jn + jp


def terminal_current(X):
    """端电流密度 [A/cm^2]（正号表示电流从左边流向右边）。

    取最左边那个面上的总通量。为什么不拿电势做差分算电流？因为只有通量
    形式才是离散方程里真正守恒的量。
    """
    return J0 * float(total_current_profile(X)[0])


def current_spread(X):
    """电流剖面在器件内的最大变化量，除以电流尺度 J0（无量纲）。

    稳态下 Jn+Jp 沿器件必须是常数，所以这个数越小越好。它既能检查残差是
    否正确，也能告诉我们"这么小的电流到底有几位可信"。
    """
    total = total_current_profile(X)
    return float(np.max(np.abs(total - total[0])))


def mass_action_error(X):
    """平衡态检查用：max |n*p/ni^2 - 1|。平衡时应该等于 0。"""
    n = N_NODES
    dens_n = exp_capped(X[0:n] - X[n:2 * n])
    dens_p = exp_capped(X[2 * n:3 * n] - X[0:n])
    return float(np.max(np.abs(dens_n * dens_p - 1.0)))


def solve_equilibrium():
    """求平衡解（两端都是 0 V）。返回 (X, 信息字典)。"""
    return newton_solve(0.0, 0.0, initial_guess(0.0, 0.0))


# ============================================================================
# 12. 偏压扫描（延续法）
# ============================================================================

def warm_start(X, delta_v):
    """把上一个偏压的收敛解平移到新偏压，作为 Newton 的初值。

    三个势场同时平移同一个量，载流子浓度不变（浓度只依赖 phi - phi_n 和
    phi_p - phi 这两个差），所以上一个解的结构被完整保留下来，Newton
    只需要重新分配多出来的那点电压。

    平移量取"从固定端到扫描端的线性斜坡"。经验数字（0.1 V 步长、从平衡
    出发、需要的迭代次数）：完全不平移 -> 50 次还不收敛；线性斜坡 -> 7 次。
    不平移的初值残差看起来很小（只有两个边界点不满足），但 Newton 之后
    要一点一点把准费米势挪到位，反而爬不动。
    """
    shift = delta_v / V0 * XT
    X_new = X.copy()
    X_new[0:N_NODES] += shift                 # phi
    X_new[N_NODES:2 * N_NODES] += shift       # phi_n
    X_new[2 * N_NODES:3 * N_NODES] += shift   # phi_p
    return X_new


def run_sweep(v_max=SWEEP_MAX, step=SWEEP_STEP, verbose=True):
    """从 0 V 扫到 v_max（左接触加正电压 = 正偏）。

    返回 (电压数组, 电流数组, 结果列表)，结果列表里每一项是
    (电压, 电流密度, 状态向量X, 信息字典)。
    """
    voltages = np.arange(0.0, v_max + step / 2.0, step)     # 含 0 和 v_max

    # 第一步：平衡解，整个扫描的起点
    X, info = solve_equilibrium()
    results = [(0.0, terminal_current(X), X, info)]
    if verbose:
        print("  %5.2f V  %s  迭代 %2d 次, 残差 %.2e" % (
            0.0, "收敛" if info["converged"] else "失败",
            info["iterations"], info["residual"]))

    # 后面每一步：用上一个解平移之后做初值
    for v in voltages[1:]:
        X = warm_start(X, step)
        X, info = newton_solve(v, 0.0, X)
        results.append((float(v), terminal_current(X), X, info))
        if verbose:
            print("  %5.2f V  %s  迭代 %2d 次, 残差 %.2e, J = %+.4e A/cm^2" % (
                v, "收敛" if info["converged"] else "失败",
                info["iterations"], info["residual"], results[-1][1]))
        if not info["converged"]:
            # 这一版不做自适应步长：某个点失败就停下，并记录失败原因
            if verbose:
                print("  在此停止: %s" % info.get("message", ""))
            break

    v_array = np.array([r[0] for r in results])
    j_array = np.array([r[1] for r in results])
    return v_array, j_array, results


# ============================================================================
# 13. 打印与自检
# ============================================================================

def describe_device():
    """打印器件参数和归一化尺度，便于核对。"""
    print("=" * 70)
    print("一维 PN 结漂移扩散模拟")
    print("=" * 70)
    print("器件长度     : %.3e cm (= %.1f um)" % (LENGTH, LENGTH * 1e4))
    print("结位置       : %.3e cm" % JUNCTION)
    print("掺杂         : NA = %.1e, ND = %.1e cm^-3" % (ACCEPTOR, DONOR))
    print("节点数       : %d  (结附近 %.0f nm 一个点, 两侧 %.0f nm)"
          % (N_NODES, H.min() * 1e7, H.max() * 1e7))
    print("温度         : %.1f K,  热电压 VT = %.6f V" % (TEMPERATURE, V0))
    print("材料         : 硅  eps_r=%.1f, ni=%.1e cm^-3, mu_n=%.0f, mu_p=%.0f"
          % (EPS_R, NI, MU_N, MU_P))
    print("归一化尺度   : L0=%.1e cm, V0=%.5f V, n0=%.1e cm^-3"
          % (L0, V0, N0))
    print("               J0=%.3e A/cm^2, tau0=%.2e s, lambda^2=%.3f"
          % (J0, TAU0, LAMBDA2))
    v_bi = V0 * np.log(ACCEPTOR * DONOR / NI ** 2)
    print("内建电势     : Vbi = VT*ln(NA*ND/ni^2) = %.6f V" % v_bi)
    print()


def check_equilibrium():
    """平衡解的各项检查（对照解析结果）。返回平衡状态。"""
    v_bi = V0 * np.log(ACCEPTOR * DONOR / NI ** 2)
    X, info = solve_equilibrium()

    print("-" * 70)
    print("平衡解检查")
    print("-" * 70)
    print("1) Newton 收敛          : %s (迭代 %d 次, 残差 %.2e)"
          % ("是" if info["converged"] else "否", info["iterations"], info["residual"]))
    # 平衡态下两个准费米势处处相等 -> n*p = ni^2 逐点成立
    print("2) 质量作用 n*p = ni^2  : 最大偏差 %.2e" % mass_action_error(X))
    # 平衡态没有电流
    print("3) 端电流应等于零      : |J| = %.2e A/cm^2" % abs(terminal_current(X)))
    print("4) 电流剖面(应恒定)    : 变化量 = %.2e (归一化)" % current_spread(X))

    # 接触电势应该正好是 ±Vbi/2
    phi = X[0:N_NODES] * V0
    print("5) 接触处静电势        : 左 %+.6f V, 右 %+.6f V" % (phi[0], phi[-1]))
    print("   解析值              : 左 %+.6f V, 右 %+.6f V" % (-v_bi / 2, v_bi / 2))

    # 峰值电场：对照耗尽近似 Emax = sqrt(q*N*Vbi/eps)
    e_face = -np.diff(phi) / H
    e_max = float(np.max(np.abs(e_face)))
    eps = EPS_R * EPS_VACUUM
    e_analytic = float(np.sqrt(Q_CHARGE * max(ACCEPTOR, DONOR) * v_bi / eps))
    print("6) 峰值电场            : %.4e V/cm" % e_max)
    print("   耗尽近似解析值      : %.4e V/cm  (偏差 %+.1f%%，过渡区可动载流子"
          % (e_analytic, 100.0 * (e_max / e_analytic - 1.0)))
    print("                          分担了部分电荷，所以数值略低，属正常)")
    peak_position = NODE_X[int(np.argmax(np.abs(e_face)))]
    print("7) 电场应集中在结附近  : 峰值位置 x = %.4e cm (结在 %.4e cm)"
          % (peak_position, JUNCTION))
    print()
    return X


def check_forward():
    """正偏检查：理想因子、饱和电流、电流守恒。"""
    print("-" * 70)
    print("正偏检查")
    print("-" * 70)
    voltages, currents, results = run_sweep(v_max=0.5, step=0.05, verbose=False)

    print("1) 电流单调上升        : %s"
          % ("是" if np.all(np.diff(currents) > 0.0) else "否"))

    # 理想因子 n = dV / (VT * ln(J2/J1))：扩散电流占主导时应该接近 1。
    # 从 i=2 开始：0 -> 0.05 V 那一段的 J(0) 是零，比值没有意义。
    ideal = []
    for i in range(2, len(voltages)):
        ratio = currents[i] / currents[i - 1]
        ideal.append((voltages[i] - voltages[i - 1]) / (V0 * np.log(ratio)))
    print("2) 理想因子 n = dV/(VT*ln(J2/J1))，扩散区应接近 1：")
    print("   0.05-0.15 V (复合电流占主导): %.3f, %.3f" % (ideal[0], ideal[1]))
    print("   0.20-0.30 V (过渡区)        : %.3f, %.3f, %.3f"
          % (ideal[3], ideal[4], ideal[5]))
    print("   0.35-0.50 V (扩散占主导)    : %.3f, %.3f, %.3f"
          % (ideal[6], ideal[7], ideal[8]))

    # 饱和电流：解析值（短二极管极限）J = q*ni^2*(Dp/(Wn*ND) + Dn/(Wp*NA))
    w_n = LENGTH - JUNCTION
    w_p = JUNCTION
    j_s_analytic = Q_CHARGE * NI ** 2 * (
        MU_P * V0 / (w_n * DONOR) + MU_N * V0 / (w_p * ACCEPTOR)
    )
    j_s_fit = currents[-1] / np.expm1(voltages[-1] / V0)
    print("3) 饱和电流 Js:")
    print("   解析式(短二极管)      : %.4e A/cm^2" % j_s_analytic)
    print("   由 %.2f V 的电流反推   : %.4e A/cm^2  (偏差 %+.1f%%)"
          % (voltages[-1], j_s_fit, 100.0 * (j_s_fit / j_s_analytic - 1.0)))

    print("4) 电流剖面守恒（稳态下 J 沿器件应为常数）:")
    for (v, j, X, info) in results:
        if v < 0.15:
            continue
        relative = current_spread(X) / abs(total_current_profile(X)[0])
        print("   %4.2f V  J = %+.4e A/cm^2   剖面相对变化 %.1e" % (v, j, relative))
    print()
    return voltages, currents, results


def main():
    describe_device()
    check_equilibrium()
    check_forward()

    print("=" * 70)
    print("偏压扫描（左接触加正电压 = 正偏）")
    print("=" * 70)
    voltages, currents, results = run_sweep()

    # 写 I-V 曲线，方便以后画图
    filename = "iv_curve.csv"
    with open(filename, "w", encoding="utf-8") as handle:
        handle.write("voltage_V,current_density_A_per_cm2,converged,iterations\n")
        for (v, j, X, info) in results:
            handle.write("%.6f,%.8e,%d,%d\n"
                         % (v, j, 1 if info["converged"] else 0, info["iterations"]))
    print()
    print("I-V 曲线已写入 %s" % filename)

if __name__ == "__main__":
    main()
