# W2 → W3 接口会审记录（A/B 签字）

> **状态**：B 已按 A 的 `b195ce5`（2026-10-23）准费米能级（QF）重构复核并同步；
> A 待签。冻结后 W3 起不得在无通知的情况下改字段名或 QF 变量定义。
> **依据**：`目的与规划.md` §14.7、§10.3；A 的 `C2_v1/src/minitcad/discretization/residual.py`
> 头注释与 `C2_v1/tests/test_qf_state.py`、`test_qf_stability.py`。
> **Owner**：B 维护本表，A 校核。

## 1. 版本与签字

| 项 | 值 |
|---|---|
| 接口版本 | `w2-w3/qf-v1`（在 2026-10-18 冻结版上增加 QF 内部表示） |
| 批准人 A / B | 待签 / B（2026-10-25 复核 `b195ce5`） |
| 变更日期 | 2026-10-25 |
| 对应 git commit | A：`b195ce5`；B：本工作区（根 `src/minitcad` 同步提交） |

## 2. 会审项（QF 版）

| 字段 | 正式生产版本（B 定义、A 校核） | 教学 SI 验证版本（A） | 核对方式 | 结论 |
|---|---|---|---|---|
| 网格 | `x_nodes_cm[M]`、`x_nodes_nd[M]`、`h_edges_nd[M-1]`、`w_interior_nd[M-2]` | `x_m[M]`，`h=diff(x)`（仅 SI 教学/适配） | `tests/test_core_mesh.py` | 已核对 |
| 公开势 | `State.u[M]` 无量纲 | `phi_V[M]`、`phi_L/R_V` | `phi = V_T * u` | 已核对 |
| **内部 QF 状态** | `QFState(u, eta_n, eta_p, ref_n, ref_p)`，`n=exp(u+ref_n+eta_n)`、`p=exp(-u-ref_p-eta_p)`，数组只读、M≥3、3M 交错 | 不进入教学 SI 验证器 | `tests/test_core_qf_state.py`：往返/只读/3M 布局 | 已核对（新增） |
| 载流子与掺杂 | 物理输入 `*_cm3`，内核 `n_nd,p_nd,Cbar_nd` | `n,p,N_D,N_A` 用 `m^-3` | A 核对符号，B 核对 cm↔SI 适配 | 已核对 |
| 介电常数 | `epsilon_f_per_cm > 0`，仅用于 `L0_cm` | `epsilon_f_per_m > 0` | B 核对 F/cm↔F/m | 已核对 |
| 公开状态布局 | `X=[u0,n0,p0,...]` 长度 `3M`（`State.pack()`） | 仅供布局练习 | `tests/test_core_state.py` | 已核对 |
| **QF 状态布局** | `Z=[u0,eta_n0,eta_p0,...]` 长度 `3M`（`QFState.pack_qf()`） | — | `tests/test_core_qf_state.py` | 已核对（新增） |
| **接触边界** | 物理：`(u_L,n_L,p_L,u_R,n_R,p_R)`；QF：`(u_L,eta_nL,eta_pL,u_R,eta_nR,eta_pR)`，左端 `eta_n=eta_p=-v`、右端 `0` | 物理边界 V | `tests/test_core_contact.py`、`test_core_qf_state.py` | 已核对（新增 QF） |
| **唯一残差** | `evaluate_qf(QFState, device, bias_V)` / `residual_qf(...)`，结点交错 `(Fu,Fn,Fp)`，6 行 Dirichlet 在 `eta_n/eta_p` 上 | 教学 `residual` 不参与耦合 | `tests/test_oracle_qf_residual.py`：内部行与旧 `residual` 一致、边界行严格为 0 | 已核对（新增） |
| **结果诊断** | `SolveResult.physics: PhysicsDiagnostics(fluxes, recombination_nd, residual_qf_nd)`，仅 `converged=True` 时出现 | — | `tests/test_core_qf_state.py`、`test_run_result.py` | 已核对（新增） |
| Poisson 行 | 内部 `F_u^nd` 全部无量纲 | 内部 `C/m^3`，边界 `V` | `tests/test_oracle_poisson.py` 非均匀二次 MMS | 已核对 |
| 非法状态 | 非正 `n/p`、非有限数据、QF 形状/参考值非法均拒绝 | 负 ε、非法 mesh/ρ/边界拒绝 | 双方交换反例 | 已核对 |

## 3. 已落地的契约（W2/W3 交付证据）

| 契约 | 实现位置 | 证据 |
|---|---|---|
| `Mesh1D` + `uniform`/`segmented` | `src/minitcad/core/mesh.py` | `tests/test_core_mesh.py` |
| 分段掺杂与 `Cbar_i` 控制体平均 | `src/minitcad/core/region.py` | `tests/test_core_region.py` |
| cm-first 材料与温度策略 | `src/minitcad/core/material.py` | `tests/test_core_material.py` |
| 接触目标与方向校验（物理 + QF） | `src/minitcad/core/contact.py` | `tests/test_core_contact.py`、`test_core_qf_state.py` |
| `State`/`QFState`/`Fluxes`/`PhysicsDiagnostics`/`SolveResult` | `src/minitcad/core/state.py` | `tests/test_core_state.py`、`test_core_qf_state.py` |
| `NormalizedDevice` 与 cm 基准尺度 | `src/minitcad/core/device.py` | `tests/test_core_device.py` |
| cm↔SI 双向适配与字段后缀契约 | `src/minitcad/io/units.py` | `tests/test_io_units.py` |
| YAML 校验与 cm/SI 同器件一致性 | `src/minitcad/io/config.py` | `tests/test_io_config.py` |
| Bernoulli（T01） | `src/minitcad/discretization/bernoulli.py`（A，已并入根包） | `tests/test_oracle_bernoulli.py` |
| SRH（T02，异寿命 `19/47`） | `src/minitcad/physics/recombination.py`（A） | `tests/test_oracle_srh.py` |
| SG 边通量（T06/T07 + QF 零电流） | `src/minitcad/discretization/sg.py`（A） | `tests/test_oracle_sg.py` |
| Poisson 非均匀（T05） | `src/minitcad/discretization/poisson.py`（A） | `tests/test_oracle_poisson.py` |
| QF 残差集成 | `src/minitcad/discretization/residual.py`（A） | `tests/test_oracle_qf_residual.py` |
| 平衡标量初值 | `src/minitcad/simulation/equilibrium.py`（A） | `tests/test_oracle_qf_residual.py` |
| `RunResult` / CLI / 22 列 CSV schema | `src/minitcad/simulation/run_result.py`、`cli.py`、`io/results_csv.py` | `tests/test_run_result.py`、`test_cli.py`、`test_results_csv.py` |

## 4. QF 重构对 W2/W3 的影响（B 的处理）

1. **公开 API 不变**：`State(u,n,p)` 仍是求解器输入/输出与 writer 的唯一物理状态；
   `contact_targets` 仍返回物理六元组。已有 W2 测试无需改写。
2. **新增内部表示**：`QFState` 与 `PhysicsDiagnostics` 进入 `core/state.py`；
   `SolveResult.physics` 为可选字段，默认 `None`，不破坏旧的位置参数调用。
3. **W3 的 A 模块并入根包**：`discretization/`、`physics/recombination.py`、
   `simulation/equilibrium.py` 从 `C2_v1/src/minitcad` 同步到 `src/minitcad`，
   保持单一包（§10.1“不可形成两套独立 Mesh/State/Residual”）。
   `C2_v1/` 保留为 A 的参考实现与对照，不再作为导入源。
4. **B 的 oracle 改为 QF 口径**：T06 增加“常数 `eta_n/eta_p` ⇒ 严格零电流”，
   并保留 T07 纯扩散与物理 SG 的交叉核对；新增 QF 残差与旧 `residual` 内部行的一致性测试。
5. **精度证据**：`tests/test_oracle_qf_residual.py::test_qf_srh_avoids_cancellation_near_equilibrium`
   在 `u=40`、`eta_n=1e-16` 时物理 `np-1` 已舍入为 0，而 QF 的 `expm1` 仍给出正速率。

## 5. 待 A 确认的开放项

1. **`ref_n/ref_p` 的冻结定义**：A 的 `solve_gummel` 以右接触为参考
   （`ref_n=log(n_R)-u_R`、`ref_p=-log(p_R)-u_R`）。B 的 `QFState` 默认 0，
   `contact_targets_qf` 显式接收 `ref_*`；两者必须写死在 W3-A7 的签名冻结中。
2. **`C2_v1` 遗留重复定义**：`C2_v1/.../solvers/gummel.py` 同时存在旧
   `solve_gummel`（约 419 行）与 QF 版（约 567 行），后者覆盖前者；
   合并到根包 W4 时必须只保留 QF 版。本次未把 `solvers/` 并入根包（属 W4）。
3. **A 的测试与实现不一致**：`C2_v1/tests/test_qf_stability.py` 仍按旧签名调用
   `residual_qf(device, bias, u, eta_n, eta_p, ref_n=..., ref_p=...)` 并期望三元组；
   当前 `residual_qf(QFState, device, bias)` 返回 `ndarray`。`C2_v1` 实测
   `4 failed, 164 passed`，需 A 更新测试或恢复兼容层。
4. **归一化模块分叉**：`C2_v1/.../physics/normalization.py` 把 `e0_v_per_cm`
   改名为 `e0_v_cm`、删除 `qr0_a_cm3`、增加 `normalize_cm/denormalize_cm`；
   根包仍用 B 的 W1 版本。需确认哪一版是 A 的正式交付（根包 W2 测试当前依赖 B 版）。
5. **`io/units.py` 分叉**：A 的 `C2_v1` 版本只有 `convert_units`，删除了
   `UnitError`、`require_suffixed_field` 与 cm/SI 后缀校验；这会违反 W2-B4。
   根包明确保留 B 版，**不得**把 A 的简化版合并进来；A 的求解器不依赖单位层。
6. **`device/mesh` 的 A 版校验增强**（`area_cm2` 可空、右端必须 0 V、
   `segmented` 必须含结区断点）：根包尚未合并；需 A 确认是否作为正式契约，
   其中“右端 0 V”与 `contact_targets` 的右端参考一致，建议 W4 前合并。

## 6. 签字

A：**待签**（请确认 §5 的 1–6 项，尤其是 `ref_n/ref_p` 与 `C2_v1` 遗留定义）

B：**已签**（2026-10-25，复核 commit `b195ce5`；根包 `pytest -q` 全绿）
