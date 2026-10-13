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
