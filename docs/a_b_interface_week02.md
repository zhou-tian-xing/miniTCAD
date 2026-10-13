# W2 接口会审记录（A/B 签字）

> **状态**：模板，待 W2 末 A/B 逐项会审后填写并签字。
> **依据**：`目的与规划.md` §14.7（W2 → W3 接口会审）、§10.3（A/B 交接标准）。
> **Owner**：B 提供本表，A 校核；冻结后 W3 起不得在无通知的情况下改字段名。

## 1. 版本与签字

| 项 | 值 |
|---|---|
| 接口版本 | （待填） |
| 批准人 A / B | （待填） / （待填） |
| 变更日期 | （待填） |
| 对应 git commit | （待填） |

## 2. 会审项

| 字段 | 正式生产版本（B 定义、A 校核） | 教学 SI 验证版本（A） | 核对方式 | 结论 |
|---|---|---|---|---|
| 网格 | `x_nodes_cm[M]`、`x_nodes_nd[M]`、`h_edges_nd[M-1]`、`w_interior_nd[M-2]` | `x_m[M]`，`h=diff(x)`（仅 SI 教学/适配） | `tests/test_core_mesh.py`：严格递增、单一尺度、往返 | 待填 |
| 势 | `State.u[M]` 无量纲 | `phi_V[M]`、`phi_L/R_V` | 核对 `phi = V_T * u` | 待填 |
| 载流子与掺杂 | 物理输入 `*_cm3`，内核 `n_nd,p_nd,Cbar_nd` | `n,p,N_D,N_A` 用 `m^-3` | A 核对符号，B 核对 cm↔SI 适配 | 待填 |
| 介电常数 | `epsilon_f_per_cm > 0`，仅用于 `L0_cm` | `epsilon_f_per_m > 0` | B 核对 F/cm↔F/m，A 负责 SI MMS | 待填 |
| 状态布局 | `X=[u0,n0,p0,...]` 长度 `3M` | 仅供布局练习 | `tests/test_core_state.py`：pack/unpack 往返、只读 | 待填 |
| Poisson 行 | 内部 `F_u^nd` 全部无量纲 | 内部 `C/m^3`，边界 `V` | A 提供，用 §3.6 的解析关系比较 | 待填 |
| 非法状态 | 非正 `n/p`、非有限数据拒绝 | 负 ε、非法 mesh/ρ/边界拒绝 | 双方交换反例，失败有可读原因 | 待填 |
| Week03 下游 | `bernoulli -> SG -> continuity -> F(X)` | SI MMS 验证器不参与耦合求解 | A 主集成，B 独立测试 | 待填 |

## 3. 已落地的契约（W2 交付证据）

| 契约 | 实现位置 | 证据 |
|---|---|---|
| `Mesh1D` + `uniform`/`segmented` | `src/minitcad/core/mesh.py` | `tests/test_core_mesh.py` |
| 分段掺杂与 `Cbar_i` 控制体平均 | `src/minitcad/core/region.py` | `tests/test_core_region.py` |
| cm-first 材料与温度策略 | `src/minitcad/core/material.py` | `tests/test_core_material.py` |
| 接触目标与方向校验 | `src/minitcad/core/contact.py` | `tests/test_core_contact.py`（T03/T04） |
| `State`/`Fluxes`/`SolveResult` | `src/minitcad/core/state.py` | `tests/test_core_state.py`（I-01） |
| `NormalizedDevice` 与 cm 基准尺度 | `src/minitcad/core/device.py` | `tests/test_core_device.py`（附录 C） |
| cm↔SI 双向适配与字段后缀契约 | `src/minitcad/io/units.py` | `tests/test_io_units.py` |
| YAML 校验与 cm/SI 同器件一致性 | `src/minitcad/io/config.py` | `tests/test_io_config.py` |

## 4. 待定开放项（会审时须给出结论）

1. `D_0` 参数化：默认 `reference_diffusivity_cm2_s`（示例 10.0）还是 `reference_mobility_cm2_per_Vs`（1000）？当前实现要求**二选一**，两者不得同时出现（第 3.5 节）。
2. `n_i(T)` 模型：W2 只接受与运行温度一致的参数；是否引入具名 `n_i(T)` 模型留待后续周。
3. 参考材料 `SILICON` 的迁移率取 1400/450（C2_v1 原型），示例 YAML 取 1350/480（第 13 节）；两张数值来源需在会审中统一口径。
