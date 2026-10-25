# 数学契约（math_contract.md）

> **Owner**：A（`目的与规划.md` §14.2 W1 的 W1-A1…W1-A6）。
> **状态**：占位。本文件由 A 在本周内填写并与 B 双签；在双签之前不得开始 `poisson.py`。
> **用途**：冻结符号、统一正方向、cm→无量纲尺度、接触条件与 SRH 归一化式，
> 使 A 的离散模块与 B 的器件/单位模块共用同一套记号。

## 待填写的内容（对应 W1-A1…W1-A6）

1. **方程符号**（第 2.2 节）：Poisson `-d(ε dφ/dx)/dx = q(p - n + C)`、`E = -dφ/dx`、
   `J_n`/`J_p` 的漂移项与扩散项符号、`dJ_n/dx = +qR`、`dJ_p/dx = -qR`。
2. **统一正方向与极性**（第 0 节决定 1、2）：左 P 右 N；`C = N_D - N_A`，P 区为负；
   `V_app = V_left - V_right`，`V_app > 0` 为正向。
3. **尺度公式与符号表**（第 3.1、3.3 节）：`V_T`、`L_0`、`D_0`、`t_0`、`E_0`、`J_0`、`R_0`、`μ_0`，
   以及 cm↔SI 双向换算表。
4. **接触条件**（第 4.1、4.2 节）：`a(C) = asinh(C/2)`、`u(0) = a(C_L) + v`、`u(ℓ) = a(C_R)`、
   `n = e^a`、`p = e^-a`。
5. **SRH**（第 2.3 节）：`r = (np - 1) / (τ̂_p (n + 1) + τ̂_n (p + 1))`，
   并声明单元测试必须使用 `τ_n ≠ τ_p`。
6. **数值压力场景与接口变量清单**（第 10.2、10.3 节）。

## 已由代码实现的对应关系（供 A 填写时引用）

| 契约内容 | 代码位置 |
|---|---|
| cm→无量纲尺度（含 `D_0` 二选一） | `src/minitcad/physics/normalization.py` |
| `V_T(T)` | `src/minitcad/physics/thermal_voltage.py` |
| 接触势与接触目标 | `src/minitcad/core/contact.py` |
| cm↔SI 换算因子与字段后缀契约 | `src/minitcad/io/units.py` |
| QF 状态与 QF 接触目标 | `src/minitcad/core/state.py`、`core/contact.py` |
| QF 唯一残差 | `src/minitcad/discretization/residual.py` |

## 附录：准费米能级（QF）内部表示（B 依 A 的 `b195ce5` 整理，待 A 确认）

A 在 2026-10-21 把求解器内部未知量从 `(u, n, p)` 改为 `(u, eta_n, eta_p)`。
公开的物理状态仍是 `State(u, n, p)`；QF 只存在于求解器内部。

$$
n_i=\exp(u_i+\text{ref}_n+\eta_{n,i}),\qquad
p_i=\exp(-u_i-\text{ref}_p-\eta_{p,i})
$$

- `ref_n`、`ref_p` 是每次求解固定的参考平移量（A 的 `solve_gummel` 以右接触锚定：
  `ref_n=ln n_R-u_R`、`ref_p=-ln p_R-u_R`）；B 的 `QFState` 默认为 0。
- 接触局部平衡时 `eta_n=eta_p=-v`（左端，`v=V_app/V_T`）、`0`（右端）。
- SG 边通量用 `expm1(diff(eta))` 避免 `n p-1` 的灾难性抵消：
  `j_{n,i+1/2}=-(m_n/h_i)n_{i+1}B(du_i)\operatorname{expm1}(-\Delta\eta_{n,i})`
  （`Δη>0` 分支；`Δη≤0` 分支对称），空穴同理。
- 中能级 SRH 的分子 `np-1` 改写为 `expm1((ref_n-ref_p)+(eta_n-eta_p))`。
- QF 版残差 `evaluate_qf` 的 6 行 Dirichlet 直接作用在 `eta_n/eta_p` 上；
  Poisson 行仍用 `n,p` 重建值。
- **注意**：`C2_v1` 中的旧 `residual(State,...)` 仅为历史/回归接口，W4 起
  求解器必须调用 `evaluate_qf`/`residual_qf`（唯一 `F(X)` 原则）。
