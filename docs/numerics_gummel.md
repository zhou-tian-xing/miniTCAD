# Gummel 与 QF 变量的数学说明（对应 `solvers/gummel.py`）

> 目的：让读者能从方程推出代码里的每个三对角矩阵，并解释为什么内部未知量取
> `(u, eta_n, eta_p)` 而不是 `(u, n, p)`。每个断言都给出可检验的依据：恒等式证明、
> 代码位置或测试名。符号与 `docs/math_contract.md`、`docs/numerical_conventions.md` 一致。

## 1. 建模：无量纲方程与离散格式

未知量是结点值 `u_i, n_i, p_i`（`i = 0…M-1`），边通量 `j_{n,i+1/2}`、`j_{p,i+1/2}`
（`i = 0…M-2`）。离散来自对控制体 `[ξ_i - h_{i-1}/2, ξ_i + h_i/2]` 的积分：

```text
Poisson      F_{u,i} = (g_{i+1/2} - g_{i-1/2})/w_i + C̄_i + p_i - n_i = 0,   g = du/dξ
电子连续性   F_{n,i} = (j_{n,i+1/2} - j_{n,i-1/2})/w_i - r_i = 0
空穴连续性   F_{p,i} = (j_{p,i+1/2} - j_{p,i-1/2})/w_i + r_i = 0
SG 通量      j_{n,i+1/2} = (μ_n/h_i)[n_{i+1}B(Δu_i) - n_i B(-Δu_i)]
             j_{p,i+1/2} = (μ_p/h_i)[p_i B(Δu_i) - p_{i+1}B(-Δu_i)]
SRH          r_i = (n_i p_i - 1) / (τ̂_p (n_i + 1) + τ̂_n (p_i + 1))
```

`B(z) = z/(e^z - 1)` 是 Bernoulli 函数，`B(0) = 1`。
这些式子的符号与 `discretization/poisson.py`、`sg.py`、`recombination.py`、
`residual.py` 的实现一一对应；`F_{u,i}` 的导数恒等式 `F_u = 0 ⇔ d²u/dξ² = -(C+p-n)`
可以直接把 `g = du/dξ` 代回验证：`(g_{i+1/2}-g_{i-1/2})/w_i` 是 `u''` 的二阶中心差分。

边界是 Dirichlet：`u_0, η_{n,0}, η_{p,0}` 与 `u_{M-1}, η_{n,M-1}, η_{p,M-1}` 由
`contact_targets_qf` 给出，不是未知量，所以每个子问题的未知量只有 `M-2` 个内部结点。

## 2. 变形：为什么取 `(u, η_n, η_p)`

定义（`ref_n, ref_p` 是一次求解内固定的平移量）：

```text
n_i = exp( u_i + ref_n + η_{n,i})        p_i = exp(-u_i - ref_p - η_{p,i})
```

**质量作用**：把两式相乘得

```text
n_i p_i = exp(ref_n - ref_p + η_{n,i} - η_{p,i})
⇒ n_i p_i - 1 = expm1( (ref_n - ref_p) + (η_{n,i} - η_{p,i}) )
```

这是恒等式，不是近似。用 `expm1` 计算的意义是：`np ≈ 1` 时右端是“小量的指数差”，
而左端是两个 `O(1)` 数相减。

**精度损失有多大**：`n_i p_i` 的浮点相对误差约为 `ε ≈ 2.2e-16`，因此 `n_i p_i - 1`
在 `|n_i p_i - 1| ≲ ε·n_i p_i` 时被完全舍入掉。`u = 350`、`η_n - η_p = 1e-16` 时
`ε·n ≈ 1e-13 ≫ 1e-16`，物理路径必然得到 `0`；QF 路径得到 `expm1(1e-16) ≈ 1e-16`。
这两种结果都由测试固定：
`tests/test_oracle_qf_residual.py::test_qf_srh_avoids_cancellation_near_equilibrium`（SRH）、
`tests/test_oracle_sg.py::test_qf_sg_resolves_a_splitting_that_the_physical_path_rounds_to_zero`（SG）。

**SG 通量的 QF 形式**：由 `n_{i+1}/n_i = exp(Δu_i + Δη_{n,i})` 与恒等式
`B(-z) = B(z) e^{z}`（把 `z/(e^z-1)` 的分子分母同乘 `e^{-z}` 即得）可推出

```text
j_{n,i+1/2} = (μ_n/h_i) n_i B(-Δu_i) expm1(Δη_{n,i})
```

即通量只依赖 `Δη_n`，与 `ref_n` 无关。空穴同理，方向相反。
`expm1(Δη)` 在 `Δη → 0` 时给出 `Δη`，这就是“同一物理状态下 QF 路径能分辨
`1e-16` 的 QF 分裂、物理路径不能”的原因。

**参考量不改变物理**：`ref_n, ref_p` 同时出现在 `n = exp(u + ref_n + η_n)` 与
`η` 的边界条件里，两者相消；`F(X)`、通量、SRH 都只依赖 `η` 的差与和 `u`。
测试 `test_internal_rows_are_invariant_to_the_qf_reference_shift` 用三组参考量验证了这一点，
所以 A 的求解器对 `ref_n` 取右端、`ref_p` 取左端不影响结果（见 `docs/a_b_interface_week02.md` §5）。

## 3. 求解：Gummel 解耦

一个外层迭代解三个子问题，每个子问题只解一个变量：

### 3.1 Poisson 子步（固定 `φ_n = ref_n + η_n`、`φ_p = ref_p + η_p`，解 `u`）

把 `n = e^{u+φ_n}`、`p = e^{-u-φ_p}` 代入 `F_u`，对 `u` 求导：

```text
∂F_{u,i}/∂u_i     = -(1/h_i + 1/h_{i-1})/w_i - n_i - p_i
∂F_{u,i}/∂u_{i±1} =  (1/w_i h_{i±1/2})
```

（`∂n_i/∂u_i = n_i`，`∂p_i/∂u_i = -p_i`。）Newton 步 `δ = -J_F^{-1} F`，代码里装配的是
`J_F` 本身、右端取 `-F`，两者等价：

| 数学量 | 代码 | 位置 |
|---|---|---|
| `diag_k = -(1/h_k + 1/h_{k+1})/w_{k+1} - n - p` | `diagonal` | `_poisson_step` |
| `∂F/∂u_{i-1}`、`∂F/∂u_{i+1}` | `left[1:]`、`right[:-1]` | 同上 |
| `δ = J^{-1}(-F)` | `_solve_tridiagonal(diagonal, left, right, -f[1:-1])` | 同上 |

内层用阻尼线搜索：`α` 从 1 折半，直到 scaled residual 下降；`α < 1e-6` 记
`POISSON_LINE_SEARCH_FAILED`，超过 35 次记 `POISSON_MAX_INNER`。

### 3.2 载流子子步（固定 `u` 与另一载流子的 QF，解 `η_n` 或 `η_p`）

电子侧把 `j` 对 `η_n` 求导：

```text
∂j_{n,i}/∂η_{n,i}     = -(μ_n/h_i) n_i B(-Δu_i)
∂j_{n,i}/∂η_{n,i+1}   =  (μ_n/h_i) n_{i+1} B(Δu_i)
```

代入 `F_{n,i}` 并冻结 SRH 分母（`τ̂_p(n+1) + τ̂_n(p+1)` 视为常数，见下），得到

```text
∂F_{n,i}/∂η_{n,i}    = [-(μ_n/h_i)n_i B(-Δu_i) - (μ_n/h_{i-1})n_i B(Δu_{i-1})]/w_i - n_i p_i/denom
∂F_{n,i}/∂η_{n,i-1}  =  (μ_n/h_{i-1}) n_{i-1} B(-Δu_{i-1})/w_i
∂F_{n,i}/∂η_{n,i+1}  =  (μ_n/h_i) n_{i+1} B(Δu_i)/w_i
```

代码装配的是 `M = -∂F/∂η`，解 `M δ = F`（等价 `δ = -J_F^{-1}F`），再 `η ← η + α δ`：

| 数学量 | 代码 | 位置 |
|---|---|---|
| `-∂F/∂η_i`（对角） | `diagonal` | `_qf_carrier_step` |
| `-∂F/∂η_{i-1}`（下次对角） | `left` | 同上 |
| `-∂F/∂η_{i+1}`（上次对角） | `right` | 同上 |
| `M δ = F` | `_solve_tridiagonal(diagonal, left, right, residual)` | 同上 |

**为什么每条边只连两个未知量**：`j_{n,i+1/2}` 只依赖 `η_{n,i}` 与 `η_{n,i+1}`，
`r_i` 只依赖本结点。所以第 `i` 行只含 `i-1, i, i+1` 三列，`M` 是三对角矩阵，
用 `scipy.sparse.diags` 装配后 `spsolve` 直接解，复杂度 `O(M)`。

**为什么冻结 SRH 分母**：分母里同时含 `n` 与 `p`，若对它也求导，本行的 `∂F/∂η`
会多出相邻结点的耦合项，三对角结构就不再成立。冻结分母相当于把这一步的
Gummel 线性化误差交给外层迭代消除；这不是精度损失，因为最终收敛判据用的是
**完整残差** `evaluate_qf`，而不是这个近似线性化。

### 3.3 外层线搜索与双判据

三个子步给 `star`，外层试 `X = old + α(star - old)`（`α` 从 1 折半）：

- 接受条件：scaled residual 下降，或已达到 `tol_residual`；
- 收敛判据（缺一不可）：`R ≤ tol_residual` **且** `D ≤ tol_update`
  （`solvers/convergence.py::is_converged`）。

`R` 与 `D` 的定义与尺度见 `docs/numerical_conventions.md` §5 与
`solvers/convergence.py` 的模块文档；要点是每行除以该行“物理上不可能更小”的量级
（`max(1, |C̄|, n, p)`、`max(1, |j|/w, |r|)`、`max(1, |target|)`），
否则在 `1e6` 量级的掺杂下，`1e-8` 的绝对残差没有意义。

## 4. 可验证的性质（每条都有测试）

| 性质 | 表达式 | 测试 |
|---|---|---|
| 常数 QF ⇒ 零电流 | `Δη = 0 ⇒ expm1(0) = 0` | `test_oracle_sg.py::test_constant_quasi_fermi_gives_zero_current` |
| 纯扩散极限 | `Δu = 0 ⇒ B(0)=1 ⇒ j_n = μ_n(n_R - n_L)/h` | `test_oracle_sg.py::test_pure_diffusion_limit_has_closed_form` |
| 平衡零电流 | `n_R = n_L e^{Δu} ⇒ j = 0` | `test_oracle_sg.py::test_equilibrium_fluxes_vanish` |
| SRH 异寿命 | `r = (np-1)/(τ̂_p(n+1)+τ̂_n(p+1))`，`τ_n≠τ_p` 时 `19/47` | `test_oracle_srh.py` |
| 参考量平移不变 | 三组 `ref_*` 内部行相同、接触行严格 0 | `test_oracle_qf_residual.py::test_internal_rows_are_invariant_to_the_qf_reference_shift` |
| 0 V 质量作用 | `max|np-1| = 2.2e-16` | `tests/test_week04_gate1.py::test_t08_equilibrium_mass_action` |

## 5. 支持范围与失败模式

- `|u| ≤ 700`（约 18 V）内 `exp` 不溢出；再往上 `Fluxes`/`State` 的有限性检查抛出
  `ValueError`，不会把 `inf`/`nan` 当电流交给外层（`test_oracle_sg.py::test_qf_sg_is_exact_up_to_u_700_and_fails_loudly_beyond`）。
- 子步失败以 `ArithmeticError` 抛给外层，外层把它记成 `SUBSTEP_FAILED: <原因>`
  并让 `SolveResult.state = None`；`sweep` 不会把失败状态当作下一偏压的起点。
- 外层达到 `max_iter` 记 `MAX_ITER`；线搜索折半到 `alpha_min` 仍不下降记
  `LINE_SEARCH_FAILED`。三者都会出现在 `comparison.csv` 的 `notes` 与
  `convergence.md` 的偏压表里。
