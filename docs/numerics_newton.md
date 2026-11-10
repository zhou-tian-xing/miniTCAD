# Newton 联立求解的数学说明（对应 `solvers/newton.py` 与 `solvers/jacobian.py`）

> 目的：从同一个 `F(z)` 推出 Newton 的每一步——Jacobian 的带宽从哪来、为什么按 9 种颜色
> 做差分、为什么先做行缩放、为什么线搜索一定会终止——并说明每一步在代码里的位置。
> 方程与离散格式沿用 `docs/numerics_gummel.md` §1–§2，不在此重复；本文只讲**解法**。
> 符号：`z = pack_qf() = [u_0, η_{n,0}, η_{p,0}, …, u_{M-1}, η_{n,M-1}, η_{p,M-1}]`，长度 `3M`。

## 1. 建模：把三个子问题换成一个联立方程

Gummel 解的是三个解耦后的子问题（`M-2` 个未知量各一次）。Newton 解的是同一个
`F(z) = 0`，只是不再拆分：

```text
F(z) = [F_u(z); F_n(z); F_p(z)] = 0,     z ∈ R^{3M}
```

其中 `F_u, F_n, F_p` 就是 `discretization/residual.py::evaluate_qf` 给出的三组行。
两端的 Dirichlet 行（`u_0, η_{n,0}, η_{p,0}` 与 `u_{M-1}, η_{n,M-1}, η_{p,M-1}`）也在这
`3M` 个方程里：它们形如 `z_j - target_j = 0`，对 `z_j` 的导数是 1，其余是 0。

**未知量与 Gummel 完全相同**，所以两者的解可以直接逐位比较；这也是 T11 能做“同一初值、
两个求解器、比同一个解”的前提。

## 2. 变形：Newton 步与行缩放

### 2.1 原始 Newton 步

```text
J(z_k) δ = -F(z_k),        z_{k+1} = z_k + α δ,      α ∈ (0, 1]
```

`J = ∂F/∂z` 是一个 `3M × 3M` 的矩阵。它按结点分块是**块三对角**的：第 `k` 行块只依赖
结点 `k-1, k, k+1` 的未知量。理由是 `F` 的每一项都只跨一条边：

| `F` 的项 | 依赖 | 因此 `J` 的第 `r` 行非零列 |
|---|---|---|
| `(g_{k+1/2} - g_{k-1/2})/w_k`，`g = du/dξ` | `u_{k-1}, u_k, u_{k+1}` | `3(k-1) … 3(k+1)+2` |
| `j_{n,k±1/2}` | `u_k, u_{k±1}, η_{n,k}, η_{n,k±1}` | 同上 |
| `r_k`（SRH） | `u_k, η_{n,k}, η_{p,k}` | 同上 |
| 接触行 `z_j - target_j` | 只有 `z_j` 自己 | 第 `j` 列 |

所以 `J` 的带宽是 9 行（3 个结点 × 3 个分量）。这一点决定了 3.1 的着色方案。

### 2.2 行缩放

Poisson 行的量级由 `∂F_{u,i}/∂u_i ≈ -n_i - p_i` 决定，在 `1e6`（约 `1e16 cm^-3`）量级；
连续性行由 `∂F_{n,i}/∂η_{n,i} ≈ -|j|/w` 决定，在 `1` 量级。直接解 `J δ = -F` 时，解向量
`δ` 的精度被大行主导：连续性方向需要的修正量约 `1e-3`，而 `1e-3/1e6 = 1e-9` 相对
`1e6` 的舍入（`ε·1e6 ≈ 2e-10`）只剩一位有效数字。表现就是“线性求解显示成功，
但外层残差不再下降”。

于是先给每行除以它的尺度 `s_r`（`solvers/convergence.py::residual_scales_qf`）：

```text
(D^{-1} J) δ = -D^{-1} F,     D = diag(s_1, …, s_{3M})
```

`D` 只缩放**方程**、不缩放未知量，所以 `δ` 仍是 `z` 的增量，不需要回代。
判据 `R = max_r |F_r| / s_r` 用的正是同一组 `s_r`——求解尺度与判据尺度必须来自同一个
函数，否则会出现“缩放后的系统解得很准，却按另一个尺度判成没收敛”。

代码位置：

| 数学量 | 代码 | 位置 |
|---|---|---|
| `s` | `residual_scales_qf(...)` | `solvers/convergence.py` |
| `D^{-1}J` | `diags(1.0 / scales, format="csc") @ jac` | `solve_newton` |
| `-D^{-1}F` | `-evaluation.f * inverse_scales` | 同上 |
| `δ` | `spsolve(jacobian_scaled, rhs_scaled)` | 同上 |

### 2.3 线性解的可信度

`spsolve` 遇到奇异矩阵只发 `MatrixRankWarning` 并返回垃圾。这里把该告警升为异常
（`warnings.simplefilter("error", MatrixRankWarning)`），再用**未缩放**的 `J`、`F`
核对后向误差：

```text
β = ‖J δ + F‖∞ / max(1, ‖F‖∞)
```

`β > 1e-6` 就判 `LINEAR_SOLVE_FAILED`。用未缩放系统算 `β` 是必要的：后向误差问的是
“`δ` 到底解没解对原方程”，用缩放后的系统会把 `D` 的量级也带进来。

## 3. 求解：Jacobian 装配、线搜索、收敛

### 3.1 九色差分

解析 Jacobian 要为 `expm1`、`B(z)` 的三种分支、SG 的两个方向各写一份导数，并且必须与
`discretization/` 保持逐字同步；一旦残差改了而 Jacobian 没跟上，Newton 的二次收敛会
悄悄退化成线性，很难定位。所以 `jacobian.py` 用差分，永远与残差同源。

差分要 `6M` 次残差求值（每列中心差分两次）。利用 §2.1 的带宽，可以按结点着色省掉大部分：

```text
颜色 c = k mod 3，分量 q ∈ {u, η_n, η_p}
→ 每一“色”是 3M/3 个列，同色结点相距 ≥ 3 个网格间距 > 带宽 1 个结点
```

同色列互不干扰，于是一次扰动整色得到该色所有列的差分，共 `3 × 3 = 9` 次求值
（3 个分量 × 3 种颜色）而不是 `6M` 次。取差分时只截取扰动结点附近的 9 行
（`3(k-1) … 3(k+1)+2`）：整向量差分的其余行理论上是 `0`，浮点上是 `1e-16` 量级的抵消
噪声，截掉既是去噪也是保持稀疏。

**步长**：`h_j = 1e-5 · max(1, |z_j|)`。中心差分的总误差是 `O(h²)` 截断加 `O(ε/h)` 舍入，
在 `h ≈ 1e-7` 处取极小；取 `1e-5` 偏向抑制舍入，因为 Newton 只需要方向足够准，
线搜索会用真实残差复核。T10 的独立 oracle 用 `1e-6`，故意与实现错开一个量级。

**越界**：整色求值失败（某个扰动点让 `exp` 溢出）时，该列退回单体中心差分；仍失败就把
步长减半重试，最多 12 次后抛 `FloatingPointError`，不静默填 0。

| 数学量 | 代码 | 位置 |
|---|---|---|
| `h_j` | `rel_step * np.maximum(1.0, np.abs(z))` | `jacobian_qf` |
| 着色 | `3 * np.arange(node_color, m, 3) + component` | 同上 |
| `(f_+ - f_-)/(2h)` | `(fp - fm) / (2.0 * steps[col])` | 同上 |
| 局部 9 行 | `start = 3*max(0, k-1)`，`end = 3*min(m, k+2)` | 同上 |

### 3.2 线搜索

Newton 方向只在解的邻域内保证下降。取 merit 函数

```text
φ(z) = ‖D^{-1} F(z)‖₂
```

`α` 从 1 开始折半，接受条件是

```text
R(z + αδ) ≤ tol_residual   或   φ(z + αδ) ≤ (1 - 1e-4 α) φ(z)
```

右式是 Armijo 条件的宽松版（只要求有下降，不要求下降满一个固定比例）。左式让“已经收敛”
的试探直接通过，避免在残差接近 `tol_residual` 时被 merit 的舍入噪声挡住。
`α` 折半到 `alpha_min` 仍不接受就记 `LINE_SEARCH_FAILED`。

**必定终止**：`α` 每次乘 `1/2`，从 1 到 `alpha_min = 1e-6` 最多 `⌈log₂(1/α_min)⌉ = 20` 次；
外层最多 `max_iter` 次。所以一次失败求解的残差求值次数有硬上界，不会挂住。

| 数学量 | 代码 | 位置 |
|---|---|---|
| `φ(z)` | `stable_l2_norm(evaluation.f / scales)` | `solve_newton` |
| Armijo | `trial_merit <= (1.0 - 1e-4 * alpha) * base_merit` | 同上 |
| `α ← α/2` | `alpha *= 0.5` | 同上 |

### 3.3 双判据

收敛要求 `R ≤ tol_residual` **且** `D ≤ tol_update`，与 Gummel 共用
`scaled_update_physical`：`D` 定义在物理量 `(u, n, p)` 上，不是 QF 变量上。
原因是 `η` 的零点由锚点 `ref_n/ref_p` 决定，换一个锚点就换一个数值，对 `η` 取相对变化
没有物理含义。

`ref_n = log(n_R) - u_R`（右端电子）、`ref_p = -log(p_L) - u_L`（左端空穴），与
`solve_gummel` 逐位相同。锚点只影响 `exp` 的取整，残差对它是平移不变的
（`test_oracle_qf_residual.py::test_internal_rows_are_invariant_to_the_qf_reference_shift`）。

### 3.4 与 Gummel 的计算量对比

| | Gummel | Newton |
|---|---|---|
| 每次外层迭代 | 1 次标量 Poisson 内层 + 2 次三对角解 | 1 次 Jacobian（9 次残差求值）+ 1 次 `3M×3M` 稀疏分解 |
| 每步复杂度 | `O(M)` | `O(M)`（带状），常数更大 |
| 收敛阶 | 线性（电流越大越慢） | 解的邻域内二次 |
| 依赖初值 | 弱（解耦后每步都逼近单调） | 强（阶跃过大时线搜索会退化，见 §5） |

## 4. 可验证的性质

| 性质 | 表达式 | 测试 |
|---|---|---|
| `J` 是真导数 | `‖J_A d − J_FD d‖∞ / max(1,‖J_FD d‖∞) = 1.65e-11` | `tests/test_t10_jacobian.py` |
| 带宽假设 | 9 色差分与稠密中心差分一致（上一条的推论） | 同上 |
| 两法同解 | `max\|Δu\| = 3.4e-11`（`8.7e-13 V`）、终端电流相对差 `7.5e-12` | `tests/test_newton_adapter.py` |
| 签名一致 | 参数名与默认值逐项相同 | `tests/test_newton_adapter.py::test_solver_signatures_match` |
| 失败不抛异常 | `state is None`、`physics is None`、`iteration_history` 可读 | `test_failure_returns_a_result_instead_of_raising` |
| 行缩放必要 | 未缩放时连续性方向的修正落到 `1e-16` 以下 | 见 §2.2；`tests/test_week04_gate1.py` 的 `max scaled residual 6.7e-10` |

## 5. 支持范围与失败模式

- **阶跃上限**：从 0.05 V 的解直接跳到 0.10 V 会 `LINE_SEARCH_FAILED`。原因是接触投影
  一瞬间把左端 `η_n` 抬 `Δv = 1.93`，而内部结点还是旧值，于是第一条边的
  `j_p ≈ μ_p p expm1(Δη)/h` 被算成 `1e12` 量级的虚假通量；该行的尺度也跟着变成 `1e12`，
  merit 失去分辨力，`α` 被压到 `1e-3`，60 次迭代也走不出这一层。
  实测可用的阶跃为 `ΔV ≤ 0.02 V`（`reports/week05_acceptance.md` 表 3）。
  这是**继续法**的问题，不是 Jacobian 的问题（T10 在该状态下仍通过），因此用
  `simulation/sweep.py` 的自适应步长解决，而不是放宽判据或改 Newton。
  附带后果：0.10–0.15 V 处在失败边界附近，自适应路径的**尝试序列**对 `bias_V` 末位
  敏感，只有“是否到达目标偏压”是可复现的。
- `|u| ≤ 700` 之外 `exp` 溢出：`QFState.to_physical` 与残差的有限性检查直接抛出，
  Newton 记 `INVALID_INITIAL_STATE` 或让线搜索拒绝该试探。
- 线性系统病态：后向误差 `β > 1e-6` 记 `LINEAR_SOLVE_FAILED`。
- 线搜索耗尽：记 `LINE_SEARCH_FAILED`；迭代超限记 `MAX_ITER`。三者都写进
  `SolveResult.reason`、CSV 的 `notes` 与 `convergence.md` 的偏压表。
