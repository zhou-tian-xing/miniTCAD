# W2（第 2 周）B 任务清单、验收要求与开始方案

> **本文用途**：把 `目的与规划.md` 中属于 **B（器件 / 工作流 / 独立验证 / 交付）** 的第 2 周任务与验收要求单独抽出，并结合仓库现有代码 `C2_v1/` 给出一个可直接开工的方案。
> **来源**：`目的与规划.md` §14.2（W2 周块）、§14.4（Week 02 五日执行表 B 列）、§14.6（Week 02 测试矩阵 I-01）、§14.7（W2→W3 接口会审）、§11.2（T03/T04）、§5.1/§5.3（网格）、§2.4（材料与温度）、§3.3（cm↔SI）、§10.1/§10.2（目录与接口契约）。
> **不覆盖**：A 在第 2 周的任务（`bernoulli.py`、`poisson.py` 及其测试）。本文只列出与 B 有接口依赖的部分。
> **阶段定位**：W1–W2 为「基础准备」。W2 结束时 B 必须交出可被 A 直接调用的 `Mesh1D / State / NormalizedDevice / 单位适配`，并冻结 Mesh/State/残差的单位契约。

---

## 1. B 在 W2 的定位与冻结口径

- B 负责的产物是 **A 的数值内核所需要的器件侧数据结构与单位边界**，不是物理公式本身。
- W2 末必须冻结（第 10.3 节交接标准、第 14.7 节）：

| 冻结对象 | 提供者 → 接收者 | 冻结后可用的证据 |
|---|---|---|
| `Mesh1D`（cm 网格 + 无量纲网格 + 控制体） | B → A | 严格递增、结点/断点、`M≥3`、cm↔nd 往返测试 |
| `NormalizedDevice`（含 cm 基准尺度） | B → A | 尺度由当前 `T` 计算，写入 `convergence.md` 可复现 |
| Dirichlet 接触目标 `contact_targets` | B → A | T03/T04 接触与极性测试 |
| `State`（`u,n,p` 无量纲） | B → A | 不可原地改写、3M 交错 pack/unpack 往返（I-01） |
| SI↔cm 双向适配 | B → A / 最终用户 | cm↔SI 双向往返在声明容差内复原 |

- 冻结之后，W3 起不得在无通知情况下改字段名（第 14.9 节 W2 末行）。

---

## 2. B 的 W2 任务清单（逐条要求与判据）

### W2-B1　`src/minitcad/core/mesh.py`：`Mesh1D` 与两个生成器

**必须完成的动作**：实现 `Mesh1D` 数据结构，并提供 `uniform(M)` 与 `segmented(breakpoints, per_segment_nodes)` 两个确定性生成器。

**具体要求（第 5.1、5.3 节）**

1. 采用 0-based 索引 `i=0,…,M-1`，`0=ξ_0<ξ_1<⋯<ξ_{M-1}=ℓ`，边间距 `h_i=ξ_{i+1}-ξ_i>0`（`i=0,…,M-2`）。
2. 内部结点控制体宽度 `w_i=(h_{i-1}+h_i)/2`（`1≤i≤M-2`）；注意 `w_i` 的索引必须与当前结点 `i` 两侧的 `h_{i-1}`、`h_i` 一致，不得错取下一条边的间距。
3. 字段必须同时提供 cm 与无量纲两套：`x_nodes_cm`（形状 `(M,)`）、`h_edges_cm`（`(M-1,)`）、`w_interior_cm`（`(M-2,)`），以及 `x_nodes_nd`、`h_edges_nd`、`w_interior_nd`。
4. `M≥3`；端点精确覆盖 `[0,L]`；**结区坐标必须包含在网格断点中**。
5. `Mesh1D` 为冻结结构：一旦构造，`mesh` 与 `device` 不允许被 Newton 迭代修改（第 5.3 节）。

**交付物**：`src/minitcad/core/mesh.py`。

**完成判据**：`x_nodes_cm` 严格递增；`h_edges_cm` 长度为 `M-1`、`w_interior_cm` 长度为 `M-2`；`segmented` 生成的结区断点确实落在结点上；cm 网格与 nd 网格可由单一尺度互相换算；`M<3`、重复/倒序结点直接报错。

**依据**：第 5.1、5.3 节。

### W2-B2　`src/minitcad/core/region.py` 与 `src/minitcad/core/material.py`

**必须完成的动作**：实现分段掺杂 `region.py`，以及 **cm-first** 的材料定义与单位校验 `material.py`。

**具体要求（第 2.1、2.4、10.2 节）**

1. 材料字段采用带单位后缀的命名：`epsilon_relative`、`intrinsic_density_cm3`、`mobility_n_cm2_per_Vs`、`mobility_p_cm2_per_Vs`、`lifetime_n_s`、`lifetime_p_s`、`bandgap_eV`；绝对介电常数 `ε_s=ε_r·ε_0`，其中 `ε_0≈8.8541878128e-14 F/cm`（不得把 SI 的 `F/m` 数值直接代入 cm 公式）。
2. 参数来源必须可追溯：至少输出 `material_parameter_policy`（取值 `explicit_at_T` 或已审核的 `named_model`）、`temperature_K` 以及提供的 `n_i, μ_n, μ_p`。
3. **若输入参数明示的标定温度与运行温度 `T` 不一致，校验阶段必须拒绝，不得自动修复**（第 2.4 节）。
4. `region.py` 表达任意分段 Si 掺杂 `C=N_D-N_A`，要求 `N_D,N_A≥0`；对以下情况必须在配置校验阶段报错：区域重叠、空段、负掺杂、未覆盖 `[0,L]`、网格未覆盖器件。
5. 分段掺杂在内部结点须由 **控制体内分段积分** 得到 `C̄_i=(1/w_i)∫_{CV_i}C(ξ)dξ`；不得把恰处阶跃点的一侧值当作整个控制体的掺杂（第 5.1 节）。

**交付物**：`src/minitcad/core/material.py`、`src/minitcad/core/region.py`。

**完成判据**：非法 region 组合全部显式失败；温度不一致参数被拒绝；`C̄_i` 可由控制体积分独立复算。

**依据**：第 2.1、2.4、5.1、10.2 节。

### W2-B3　`core/contact.py`、`core/state.py`、`core/device.py`

**必须完成的动作**：实现接触目标与方向提示、状态容器，以及配置对象与归一化器件。

**具体要求（第 4.2、7.1、10.2 节）**

1. `contact_targets(device, bias_V)` 返回六元组 `(u_L, n_L, p_L, u_R, n_R, p_R)`；`bias_V` 单位始终是 **V**，函数内部转成 `v=V_app/V_T`。
2. 接触边界公式：`u(0)=a(C_L)+v`、`u(ℓ)=a(C_R)`；`n=e^{a(C)}`、`p=e^{-a(C)}`、`a(C)=asinh(C/2)`。左 P 要求 `a(C_L)<0`，右 N 要求 `a(C_R)>0`；右端为 0 V 参考。
3. 右端外加电压为零**不表示** `u(ℓ)=0`：在所选本征参考下右接触仍有 `a(C_R)`。
4. `State`：`u`、`n`、`p` 均为形状 `(M,)` 的数组，`u` 无量纲、`n,p` 无量纲且严格为正；必须是冻结对象，**数组需防止跨模块原地修改**。
5. 状态布局采用结点交错排序 `X=[u_0,n_0,p_0,…,u_{M-1},n_{M-1},p_{M-1}]^T∈R^{3M}`，并提供等价的 pack/unpack。
6. `SolveResult` 至少含 `bias_V`、`state`、`converged`、`iterations`、`residual_norm`、`update_norm`、`reason`、`iteration_history`（第 10.2 节）。
7. `DeviceConfig` 收集 YAML/API 字段（`units.system`、`device`、`material`、`regions`、`mesh`、`normalization`、`contacts`、`solver`、`sweep`、`output`）；`NormalizedDevice` 保存归一化网格/掺杂/迁移率/寿命，以及 **cm 基准尺度** `L0_cm / ni_cm3 / D0_cm2_s / J0_A_cm2 / R0_cm3_s` 与温度等可复现元数据。

**交付物**：`src/minitcad/core/contact.py`、`src/minitcad/core/state.py`、`src/minitcad/core/device.py`。

**完成判据**：接触六元组满足 `np=1`、`p-n+C=0`；`State` 的原地写入尝试被阻止或报错；`NormalizedDevice` 的 cm 尺度与附录 C 的解析值一致。

**依据**：第 4.1、4.2、7.1、10.2 节。

### W2-B4　`src/minitcad/io/units.py`：SI↔cm 双向适配

**必须完成的动作**：实现单位制和双向换算，并拒绝无单位输入。

**具体要求（第 3.3 节）**

1. 配置须声明 `units.system: cm` 或 `units.system: si`；项目默认工程契约为 **cm**。
2. 必须保留 cm↔SI 双向换算，换算因子固定为：

| 物理量 | cm → SI | SI → cm |
|---|---|---|
| 长度 | `×1e-2` | `×1e2` |
| 面积 | `×1e-4` | `×1e4` |
| 载流子 / 掺杂浓度 | `×1e6` | `×1e-6` |
| 迁移率 | `×1e-4` | `×1e4` |
| 扩散系数 | `×1e-4` | `×1e4` |
| 电场 | `×1e2` | `×1e-2` |
| 电流密度 | `×1e4` | `×1e-4` |
| 体复合率 / 数密度变化率 | `×1e6` | `×1e-6` |
| 体电荷密度 | `×1e6` | `×1e-6` |
| 绝对介电常数 | `×1e2` | `×1e-2` |
| 电势 / 电压 / 电流 / 寿命 / 温度 | 数值不变 | 数值不变 |

3. `units.system: cm` 只能用 `_cm/_cm2/_cm3/_cm2_per_Vs/_cm2_s` 字段；`units.system: si` 只能用 `_m/_m2/_m3/_m2_per_Vs/_m2_s` 字段；**二者不得在同一 YAML 中混用**。
4. **无单位的 `density=1e16` 一律拒绝**，绝不按数值大小猜测单位。
5. SI 输入经输入适配层变换为等价 cm 参数后，进入同一条归一化路径；不另造求解内核。
6. 单元测试须覆盖 `cm→SI→cm`、`SI→cm→SI`、`cm→无量纲→cm` 与 `SI→适配 cm→无量纲→恢复 cm→SI`；另须以相同物理器件验证 cm 输入与 SI 输入产生一致的无量纲状态。

**交付物**：`src/minitcad/io/units.py` 及其往返测试。

**完成判据**：四类往返在声明容差内复原；无单位输入被拒绝；cm 与 SI 输入得到同一 `F(X)`。

**依据**：第 3.3 节。

### W2-B5　接触解析基准（供 T03/T04 使用）

**必须完成的动作**：提供接触的解析基准值，作为 A 的实现之外的独立对照。

**具体要求（第 11.2 节 T03/T04）**

1. 对多组 `C_L<0<C_R`，给出独立计算的理论值：`np=1`、`p-n+C=0`、`u_L=a(C_L)+v`；中性误差与 `|np-1|` 需 `≤1e-10`。
2. T04 极性检查：对比 `V=0` 与 `V>0`，势垒 `u_R-u_L` 应减少 `v`；接触电势偏差 `≤1e-10`（无量纲）；错误方向（左 N 右 P）必须被拒绝。

**交付物**：接触基准测试（`tests/` 内），且该基准不得通过调用 A 的实现来产生期望值。

**完成判据**：T03、T04 通过；PN 方向反转时抛出 `INVALID_PN_ORIENTATION` 并标出两端 `C` 的单位与实际值。

**依据**：第 4.3、11.2 节。

### W2-B6　会审并冻结单位契约，记入 `docs/a_b_interface_week02.md`

**必须完成的动作**：按第 14.7 节的接口会审表，与 A 逐项签署 Mesh/势/载流子与掺杂/介电常数/状态布局/Poisson 行/非法状态/Week03 下游，并把结果写入 `docs/a_b_interface_week02.md` 与 `docs/numerical_conventions.md`。

**具体要求（第 14.7 节）**

| 会审项 | 正式生产版本（B 定义、A 校核） | 教学 SI 验证版本 |
|---|---|---|
| 网格 | `x_nodes_cm[M]` 与 `x_nodes_nd[M]`、`h_edges_nd[M-1]`、`w_interior_nd[M-2]` | `x_m[M]` 单位 m，`h=diff(x)`（仅 SI 教学/适配） |
| 势 | `State.u[M]` 无量纲 | `phi_V[M]`、`phi_L/R_V` |
| 载流子与掺杂 | 物理输入 `n_cm3,p_cm3,ND_cm3,NA_cm3`，内核 `n_nd,p_nd,Cbar_nd` | `n,p,N_D,N_A` 用 `m^-3` |
| 介电常数 | `epsilon_F_cm>0` 用于 `L0_cm` | `epsilon_F_m>0` |
| 状态布局 | `X=[u0,n0,p0,…]` | 仅供布局练习 |
| 非法状态 | 非正 `n/p`、非有限数据拒绝 | 负 ε、非法 mesh/ρ/边界拒绝 |

**交付物**：`docs/a_b_interface_week02.md`（含接口版本、批准人、变更日期、A/B 签字）。

**完成判据**：双方交换反例，失败有可读原因；接口表 A/B 签字。

**依据**：第 10.3、14.7 节。

---

## 3. B 在 W2 负责的测试与验收项

Week 02 专项门槛（第 14.2 节）共九项，其中属于 B 的部分如下；任一项失败不得进入 W3 集成。

| 验收项 | Owner | 判据 | 依据 |
|---|---|---|---|
| `T03_CONTACT` | B | `np=1`、`p-n+C=0`、`u_L=a(C_L)+v`；相对中性误差与 `\|np-1\|` `≤1e-10` | §11.2 |
| `T04_POTENTIAL_POLARITY` | B | `V>0` 时 `u_R-u_L` 减少 `v`；接触电势偏差 `≤1e-10`；错误方向拒绝 | §11.2 |
| `I-01` A/B 交接 | B | `Mesh/State` 不可原地改写；3M 交错 pack/unpack 往返值/形状/复制语义正确 | §14.6 |
| cm↔nd 往返 | B | cm 网格/掺杂归一化后可复原 | §14.2 |
| cm↔SI 往返 | B | 四类往返在容差内复原；同器件 cm/SI 输入得同一无量纲状态 | §3.3 |
| 异常注入（配置/网格侧） | B | 重复/倒序结点、非正 `ε`、错误形状、region 重叠/缺口全部显式失败 | §5.4、§14.6 |

由 A 负责、但 B 需提供数据的项：`T01_BERNOULLI`、`T05_POISSON_NONUNIFORM`、非均匀二次 MMS、均匀正弦加密（`[3.8,4.2]`）。B 的职责是**不重用 A 的实现**去生成这些对照（第 14.4 节 Day 4）。

---

## 4. 现状盘点：仓库与 `C2_v1/`

### 4.1 仓库现状

| 目标产物（依据 W1/W2） | 当前是否存在 | 说明 |
|---|---|---|
| `src/minitcad/` 包骨架 | 否 | 尚无 `src/` |
| `pyproject.toml` | 否 | W1-B2 未落地 |
| `docs/math_contract.md` | 否 | W1-A1 未落地 |
| `examples/pn_left_p_right_n.yaml` | 否 | W1-B5 未落地 |
| `docs/a_b_interface_week02.md` | 否 | W2-B6 待建 |
| `convergence.md` 模板 | 否 | W1-B4 未落地 |
| `C2_v1/`（物理原型） | **是** | 13 项测试通过（已实测） |
| `ai_example/`（单文件示例） | 是 | 与本规范无直接接口关系 |

结论：**W1 的 B 交付物基本未落地，W2 不能直接从零开始编码**，需在 Step 0 补齐最小前置件。

### 4.2 `C2_v1/` 模块与目标位置映射

`C2_v1/` 已实现的部分与目标 `src/minitcad/` 的对应关系（Owner 按第 10.1 节）：

| `C2_v1/` 现有文件 | 内容 | 目标位置 / Owner | 处置建议 |
|---|---|---|---|
| `constants.py` | `q,kbj,kbv,epsilon_0` | `src/minitcad/core/constants.py` | 直接搬运；`epsilon_0` 已是 `F/cm`，符合第 2.4 节 |
| `materials/base.py` | `Material` 冻结 dataclass | `src/minitcad/core/material.py`（B） | **作为种子**，补单位后缀命名与温度一致性校验 |
| `materials/silicon.py` | 硅的 300 K 参考值 | `core/material.py`（B） | 搬运，字段重命名 |
| `materials/registry.py` | `get_material(name)` | `core/material.py`（B） | 搬运 |
| `PhysicsBase/thermal_voltage.py` | `V_T(T)`，非正温度报错 | `physics/thermal_voltage.py`（A） | 由 A 接管 |
| `PhysicsBase/IntrinsicCarrier.py` | `n_i` 的 Protocol 与两种实现 | `physics/`（A） | 由 A 接管 |
| `PhysicsBase/normalization.py` | `NormalizationScales` + `build_normalization_scales` | `physics/normalization.py`（A） | 由 A 接管 |
| `PhysicsBase/mobility.py` | `MobilityModel` 抽象基类 + 常量实现 | `physics/mobility.py`（A） | 由 A 接管 |
| `PhysicsBase/transport.py` | Einstein 关系 `D=μV_T` | `physics/transport.py`（A） | 由 A 接管 |
| `PhysicsBase/SRHrecombination.py` | SRH 复合率 | `physics/recombination.py`（A） | 由 A 接管 |
| `tests/test_PhysicsBase.py` | 13 项测试 | `tests/` | 拆分为 `test_thermal_voltage.py` 等，作为首批回归 |

**已实测结论**：`C2_v1` 的 `IntrinsicCarrier.py` 中 `ReferenceIntrinsicCarrierModel` 在温度不等于材料参考温度时直接抛错。这与第 2.4 节「标定温度与运行温度不一致时拒绝、不得自动修复」的策略一致，B 的 W2 可直接沿用该策略作为默认；相对地，`ConstantIntrinsicCarrierModel` 会忽略传入温度而恒返回参考 `n_i`，属于宽松变体，建议限制使用并在 `notes` 标注。

**`C2_v1` 相对 W2-B 范围的缺口**：无 `mesh.py`、无 `region.py`、无单位制声明与校验、无 `contact.py`、无 `state.py`、无 `device.py`、无 `io/units.py`、无控制体掺杂平均、无 `pyproject.toml`，且使用扁平导入（`from constants import q`）而非包导入。

### 4.3 需要先与 A 对齐的开放项

1. **参考扩散系数 `D_0` 的选择冲突（必须先定）**：`C2_v1/physics.md` 第 11 节与 `normalization.py` 采用 `μ_0=1000 cm²/(V·s)`、`D_0=μ_0·V_T`；而 `目的与规划.md` 第 13 节示例 YAML 用 `reference_diffusivity_cm2_s: 10.0`（固定 `D_0`）。第 3.5 节规定**两种基准只能选其一**，且报告必须记录生效的 `D_0` 与 `μ_0`。本项必须在 Step 1 前定案。
2. **字段命名映射**：`C2_v1` 用 `intrinsic_concentration_ref_per_cm3` / `electron_mobility_ref_cm2_per_vs`；规范用 `intrinsic_density_cm3` / `mobility_n_cm2_per_Vs`。需确定统一命名，避免同一量两种写法。
3. **`n_i(T)` 模型**：`C2_v1` 只有参考温度常量，没有温度模型。W2 阶段建议按第 2.4 节「用户在当前温度显式提供或经评审的具名模型」处理，**不得**用 300 K 数值按温度外推；是否存在具名 `n_i(T)` 模型留待后续周决定。
4. **包布局**（见 6.0 的决策 D1）。
5. **`State` 数组的只读保证方式**（`dataclass(frozen=True)` 只冻结属性不冻结 `ndarray` 内容）：需约定用 `setflags(write=False)` 或返回副本，并在 I-01 中测试。

---

## 5. 开始方案

### 5.0 决策 D1：包布局

**建议**：新建 `src/minitcad/` 作为**唯一正式包**，`C2_v1/` 保留为「捐赠者与回归基准」，本轮不删除、也不再新增功能。

理由：第 10.1 节要求「不可形成两套独立 `Mesh/State/Residual`」；第 10.2 节的签名契约与第 14.7 节的会审表都以 `src/minitcad/...` 为路径；后续 A 的 GUI 需要 `python -m minitcad.gui.app`，即必须有名为 `minitcad` 的可导入包。若改为「把 `C2_v1` 就地改造成等价布局」，仍需把目录名改成 `minitcad` 并支持包导入，成本更高。

### 5.1 Step 0　骨架与 W1 缺口补齐（半天）

- 建目录：`src/minitcad/{core,physics,discretization,solvers,simulation,verification,io,gui}`、`tests/`、`examples/`、`docs/`、`reports/`，各含 `__init__.py`。
- 写 `pyproject.toml`：`requires-python=">=3.10"`、NumPy、SciPy，`dev` extras 含 pytest，按 `src/` 布局配置 setuptools 发现路径。
- 建 `docs/math_contract.md` 占位（符号表与尺度公式，与 A 共同补全）。
- 建 `examples/pn_left_p_right_n.yaml`，顶部三行写明左 P 右 N、`V_app=V_left−V_right`、正向 `V_app>0`。
- 建 `convergence.md` 模板（十个固定章节）。

**完成判据**：`python -c "import minitcad"` 成功；`pytest` 能收集空测试目录且不报错。

### 5.2 Step 1　`constants.py` + `material.py`（含单位校验）

- 搬运 `C2_v1/constants.py` → `src/minitcad/core/constants.py`。
- 以 `C2_v1/materials/base.py` 为种子改写 `core/material.py`：字段改名到带单位后缀命名，`ε_s=ε_r·ε_0`（`ε_0` 取 `F/cm`）。
- 加入**单位制声明与温度一致性校验**：`material_parameter_policy`、`temperature_K`，标定温度与运行温度不一致即拒绝。
- 沿用 `ReferenceIntrinsicCarrierModel` 的「温度不一致抛错」策略作为默认。

**完成判据**：`materials.silicon` 的 `ε_s` 等于 `1.0359399740976e-12 F/cm`（附录 D.2）；温度不一致输入被拒绝。

### 5.3 Step 2　`mesh.py` + `region.py` + 控制体掺杂平均

- 实现 `Mesh1D` 与 `uniform(M)`、`segmented(breakpoints, per_segment_nodes)`；同时产出 cm 与 nd 两套字段。
- 实现区域分段掺杂，并实现 `C̄_i=(1/w_i)∫_{CV_i}C(ξ)dξ` 的控制体平均接口。
- 异常注入测试：重复/倒序结点、`M<3`、region 重叠/缺口/负掺杂。

**完成判据**：`uniform` 与 `junction-refined` 网格均可构造；非均匀网格的二次制造解测试所需 `h_left/h_right` 取自真实间距；cm↔nd 往返一致。

### 5.4 Step 3　`io/units.py` 双向适配

- 按第 3.3 节的换算因子表实现 cm↔SI，并对字段后缀做严格校验（不混用、不猜单位）。
- 写四类往返测试（`cm→SI→cm`、`SI→cm→SI`、`cm→nd→cm`、`SI→适配 cm→nd→恢复 cm→SI`）。

**完成判据**：四类往返在容差内复原；同一物理器件的 cm 输入与 SI 输入得到同一无量纲状态。

### 5.5 Step 4　`contact.py` + `state.py` + `device.py`

- `contact_targets` 返回 `(u_L,n_L,p_L,u_R,n_R,p_R)`；实现 `a(C)=asinh(C/2)`。
- `State` 冻结且数组只读（Step 0 决策 D1 的第 5 条开放项落地）。
- `DeviceConfig` / `NormalizedDevice`，后者保存 cm 基准尺度与温度元数据。

**完成判据**：接触六元组满足 `np=1`、`p−n+C=0`；`NormalizedDevice` 的 `L0_cm=4.08845452e-3 cm` 与附录 C 一致。

### 5.6 Step 5　测试收口与契约冻结

- 补齐 T03、T04、I-01、cm↔nd、cm↔SI、异常注入测试。
- 与 A 会审并填写 `docs/a_b_interface_week02.md`，双方签字。
- 写 `reports/week02_acceptance.md` 内部周报（按第 14.7 节模板）。

**完成判据**：第 3 节所列 B 的验收项全部通过；接口契约冻结；A 可用 `Mesh1D/State/NormalizedDevice` 驱动其 `poisson.py`。

### 5.7 建议节奏（对应第 14.4 节 Day1–Day5 B 列）

| 天 | 与上述 Step 的对应 | 当天可核对物 |
|---|---|---|
| Day 1 | Step 0 + Step 1 | 包可导入；材料单位校验与温度拒绝的测试 |
| Day 2 | Step 1 收尾 + 提供 SI 教学样本所需网格/电荷 | 接触电势验证样本；确认内部 `C/m³` 与边界 `V` 仅限 SI 教学 |
| Day 3 | Step 2 | `Mesh1D` 的 `x_nodes_nd/h_edges_nd/w_interior_nd` 与控制体掺杂平均 |
| Day 4 | Step 3 | 不重用 A 的 Poisson 实现，用解析正弦/独立线代交叉核对 |
| Day 5 | Step 4 + Step 5 | 四项必交文件复核、cm↔无量纲与 cm↔SI 双向测试、联合签核 |

---

## 6. 明确不做（W2 范围边界）

- 不实现 SG 电子/空穴、SRH 集成、continuity、耦合 `F(X)`、Gummel、Newton —— 这些属于 A 的 W3 起内容（第 14.4 节明示 W2 只承诺「稳定 Bernoulli + Poisson 离散及其独立证明」）。
- 不修改 `C2_v1` 中已通过测试的物理公式，只做搬运、改名与单位校验增强。
- 不引入跨材料通量模型、变量介电常数、场依赖迁移率（第 5.4 节「可选扩展」之外）。
- 不新增第四份验收文件（第 12 节：正式 solve 固定三文件）。

---

## 7. W2 完成定义（B 侧 DoD 清单）

- [ ] `src/minitcad/` 可导入，`pyproject.toml` 可 `pip install -e ".[dev]"`。
- [ ] `Mesh1D` 提供 `uniform/segmented` 与 cm、nd 两套字段，`M≥3`、严格递增、结区在断点内。
- [ ] `core/material.py` 为 cm-first，含 `_cm/_cm2/_cm3` 后缀命名与温度一致性校验。
- [ ] `region.py` 支持分段掺杂与控制体平均 `C̄_i`，非法配置全部报错。
- [ ] `io/units.py` 四类往返通过，无单位输入被拒绝。
- [ ] `contact_targets` 与 `State` 可被 A 直接调用；`NormalizedDevice` 含 cm 基准尺度。
- [ ] T03、T04、I-01、cm↔nd、cm↔SI、异常注入全部通过。
- [ ] `docs/a_b_interface_week02.md` 已冻结并双签。
- [ ] `reports/week02_acceptance.md` 记录实测命令、数值、依赖版本与未通过项。

---

## 附：一次可直接执行的开工顺序

```bash
# 1) 骨架
mkdir -p src/minitcad/{core,physics,discretization,solvers,simulation,verification,io,gui} tests examples docs reports

# 2) 从 C2_v1 搬运种子（先复制，再改名/加校验，不直接改动 C2_v1）
cp C2_v1/constants.py src/minitcad/core/constants.py
cp C2_v1/materials/base.py src/minitcad/core/material.py
cp C2_v1/materials/silicon.py src/minitcad/core/silicon.py

# 3) 验证包可导入
python -c "import minitcad; print(minitcad.__file__)"

# 4) 确认前置冲突已定案（D_0 基准二选一）后再写 mesh.py / units.py
```

> 顺序建议：**先定 `D_0` 基准 → 再写 `units.py` 与 `NormalizedDevice`**。
> 注意：长度基准 `L_0=sqrt(ε V_T/(q n_i))` 只依赖 `T`、`n_i`、`ε`，**不依赖 `D_0`**（第 3.1 节），
> 所以 `mesh.py` 与 cm↔nd 往返不受 `D_0` 选择影响，可先写。受 `D_0` 影响的是
> `t_0 = L_0²/D_0`、`J_0 = q D_0 n_i/L_0`、`R_0 = D_0 n_i/L_0²` 以及 `μ_0 = D_0/V_T`；
> 这些量进入 `NormalizedDevice` 的电流/时间/寿命尺度，须在写它们之前定案，
> 否则电流与寿命相关的往返测试会返工。
