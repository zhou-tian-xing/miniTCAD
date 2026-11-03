# 数值约定与单位契约（numerical_conventions.md）

> **Owner**：B（维护）；**校核**：A。
> **依据**：`目的与规划.md` §14.7（接口会审表）、§3.1/§3.3（尺度与 cm↔SI）、§5.1/§5.3（网格与控制体）、
> §7.1（状态布局与唯一残差）、§2.3（SRH）、§11.2（T01–T07 阈值）、§12.3（convergence.md）。
> **配套文件**：`docs/a_b_interface_week02.md`（会审签字与开放项）、`docs/math_contract.md`（A 的符号表）、
> `docs/convergence_template.md`（每次运行的收敛报告模板）。
> **接口版本**：`w2-w3/qf-v1`（在 2026-10-18 冻结版上增加准费米能级内部表示）。
> **变更日期**：2026-11-03。**批准**：B 已签；A 待签（见 `a_b_interface_week02.md` §6）。

本文件只记录**机器可检查的约定**：每条约定都给出实现位置与断言它的测试。
与 `docs/math_contract.md` 冲突时，以 A 的符号表为准，并在此文件登记差异。

---

## 1. 坐标系、极性与掺杂符号

| 约定 | 取值 | 实现 | 证据 |
|---|---|---|---|
| 坐标 | 一维，`x ∈ [0, L]`，结点 `i = 0…M-1`，0-based | `core/mesh.py` | `tests/test_core_mesh.py` |
| 器件方向 | 左端 `x=0` 为 P，右端 `x=L` 为 N | `core/contact.py` | `tests/test_core_contact.py` |
| 净掺杂 | `C = N_D - N_A`，P 区 `C < 0`，N 区 `C > 0`；`N_D, N_A ≥ 0` | `core/region.py` | `tests/test_core_region.py` |
| 偏压 | `V_app = V_left - V_right`；`V_app > 0` 为正向 | `core/contact.py`、`examples/*.yaml` 首部 | `tests/test_core_contact.py` |
| 右端参考 | 右端接地 `V_right = 0 V`；`DeviceConfig` 拒绝非零 `right_voltage_V` | `core/device.py` | `tests/test_core_device.py` |
| 接触电势 | `a(C) = asinh(C/2)`；`u(0)=a(C_L)+v`、`u(ℓ)=a(C_R)`，`v=V_app/V_T` | `core/contact.py` | T03/T04 |
| 右端电势 | `u(ℓ) = a(C_R) ≠ 0`：0 V 参考**不等于**无量纲电势为 0 | `core/contact.py` | `tests/test_core_contact.py` |

## 2. 单位：cm 输入 / 无量纲内核 / 可选 SI 输出

**唯一链路**：`cm 或 si 输入 →（适配层换算成 cm）→ 归一化（无量纲）→ 求解内核 → writer 按 `output.unit_system` 输出`。
内核不出现 `cm`、`F/m`、`m^-3` 等物理单位；SI 只在输入适配层与输出（writer）层出现（§3.3、§12.1）。

| 量 | 符号 | 定义 | 典型值（300 K，Si） |
|---|---|---|---|
| 热电压 | `V_T` | `k_B T / q` | `2.58519998e-2 V` |
| 本征浓度 | `n_i` | 材料参数（cm^-3） | `1.0e10 cm^-3` |
| 归一化长度 | `L_0` | `sqrt(ε_s V_T / (q n_i))` | `4.08845452e-3 cm` |
| 归一化扩散系数 | `D_0` | 配置直接给 `D_0`，或 `μ_0 V_T`（二选一） | `10 cm²/s` |
| 归一化时间 | `t_0` | `L_0² / D_0` | `1.67154603e-6 s` |
| 归一化电场 | `E_0` | `V_T / L_0` | `6.3234 V/cm` |
| 归一化电流密度 | `J_0` | `q D_0 n_i / L_0` | `3.91878307e-6 A/cm²` |
| 归一化复合率 | `R_0` | `D_0 n_i / L_0²` | `1.5327e19 cm^-3/s` |
| 归一化迁移率 | `μ_0` | `D_0 / V_T` | `386.82 cm²/(V·s)` |

- 绝对介电常数用 `ε_s = ε_r · ε_0`，`ε_0 = 8.8541878128e-14 F/cm`（`core/constants.py`）；
  **不得**把 `F/m` 的 8.854e-12 直接代进 cm 公式（§2.4）。
- 材料参数字段一律带单位后缀（`*_cm3`、`*_cm2_per_Vs`、`*_s`）；`units.system: cm` 与 `si` 不得在同一 YAML 混用，
  无单位字段（如裸 `density`）一律拒绝（`io/units.py`、`io/config.py`）。
- 材料标定温度与运行温度不一致时**拒绝**，不做外推（`core/material.py`，`MaterialTemperatureMismatch`）。

## 3. 网格与控制体

| 量 | 形状 | 约定 | 实现 |
|---|---|---|---|
| `x_nodes_cm/nd` | `(M,)` | 严格递增，`x_0=0`，`x_{M-1}=L`；`M ≥ 3` | `Mesh1D` |
| `h_edges_cm/nd` | `(M-1,)` | `h_i = x_{i+1} - x_i > 0`，必须等于 `diff(x)` | `Mesh1D` |
| `w_interior_cm/nd` | `(M-2,)` | `w_i = (h_{i-1} + h_i)/2`，**必须取结点两侧的边**，取错下标显式失败 | `Mesh1D` |
| cm↔nd | — | 全部字段由**同一个** `L_0` 换算；不成单一尺度即失败 | `Mesh1D` |
| 结区 | — | `segmented` 网格必须把每个结区坐标列为断点，且该坐标落在结点上；`uniform` 允许结区落在控制体内部 | `core/device.py::_require_junctions_on_nodes` |
| 结点数 | — | 必须是整数；浮点/布尔结点数被拒绝（不得静默截断） | `core/mesh.py`、`io/config.py` |
| 控制体掺杂 | `cbar_nd (M-2,)` | `C̄_i = (1/w_i)∫_{CV_i} C dξ`，跨结区的控制体取两段加权平均 | `region.py::cv_averaged_net_doping_cm3` |

## 4. 状态布局

| 状态 | 字段 | 布局 | 只读 |
|---|---|---|---|
| 公开物理状态 `State` | `u, n, p`（各 `(M,)`，`n,p > 0`） | `X = [u0,n0,p0, u1,n1,p1, …]`，长度 `3M` | `pack()/unpack()`，数组冻结 |
| 求解器内部 `QFState` | `u, eta_n, eta_p`（各 `(M,)`）+ 标量 `ref_n, ref_p` | `Z = [u0,eta_n0,eta_p0, …]`，长度 `3M` | `pack_qf()/unpack_qf()`，数组冻结 |

QF 关系（`docs/math_contract.md` 附录）：

```text
n_i = exp( u_i + ref_n + eta_{n,i})
p_i = exp(-u_i - ref_p - eta_{p,i})
```

- `ref_n/ref_p` 是每次求解固定的平移锚点；**物理残差与它无关**（`tests/test_oracle_qf_residual.py::
  test_internal_rows_are_invariant_to_the_qf_reference_shift`）。
- 接触局部平衡边界：`eta_n = eta_p = -v`（左）、`0`（右），与 `ref_*` 的取值无关
  （`contact_targets_qf`；`tests/test_core_qf_state.py`）。
- `State` 是求解器输入/输出与 writer 的唯一物理状态；QF 变量不对外暴露（§10.2）。

## 5. 残差与符号

唯一残差 `F(X) = 0`（`discretization/residual.py`）：结点交错 `(F_u, F_n, F_p)`，`F(X) ∈ R^{3M}`。
QF 生产入口是 `evaluate_qf(QFState, device, bias_V)` / `residual_qf(...)`；`residual(State, ...)` 只保留为历史/回归接口。

| 行 | 公式（无量纲） | 说明 |
|---|---|---|
| Poisson 内部 | `F_{u,i} = (g_{i+1/2} - g_{i-1/2})/w_i + C̄_i + p_i - n_i`，`g = du/dξ` | 负号已并入；`g_i = diff(u)/h` |
| Poisson 边界 | `F_{u,0} = u_0 - u_L`、`F_{u,M-1} = u_{M-1} - u_R` | 目标来自 `contact_targets` |
| 电子边通量 | `j_{n,i+1/2} = (μ_n/h_i)[n_{i+1}B(Δu_i) - n_i B(-Δu_i)]` | 正方向 = `+x`；`B` 为 Bernoulli |
| 空穴边通量 | `j_{p,i+1/2} = (μ_p/h_i)[p_i B(Δu_i) - p_{i+1}B(-Δu_i)]` | 同样的 `+x` 约定 |
| 电子连续性 | `F_{n,i} = (j_{n,i+1/2} - j_{n,i-1/2})/w_i - r_i` | `r` 为 SRH 净复合 |
| 空穴连续性 | `F_{p,i} = (j_{p,i+1/2} - j_{p,i-1/2})/w_i + r_i` | 空穴符号相反 |
| SRH | `r = (np - 1)/(τ̂_p (n+1) + τ̂_n (p+1))` | 测试必须用 `τ_n ≠ τ_p`（T02 用 `19/47`） |
| QF 边界行 | 6 行 Dirichlet 直接作用在 `eta_n/eta_p` 与 `u` 上 | 内部行仍由 `n,p` 重建值构成 |

QF 精度改写（A 的 `b195ce5`）：SG 用 `expm1(Δη)`、SRH 分子用 `expm1((ref_n-ref_p)+(η_n-η_p))`，
避开 `np-1` 在 `np≈1` 处的灾难性抵消。可复现证据（B 侧独立核对）：

- SRH：`u=40`、`η_n=1e-16` 时物理路径的 `np-1` 已为 0，QF 仍给出正速率
  （`test_qf_srh_avoids_cancellation_near_equilibrium`）。
- SG：`u=350`、`Δη=1e-16` 时物理路径 `(n_{i+1}-n_i)/h = 0.0`，QF 给出
  `5.0e137`，与 `n·expm1(Δη)/h` 相对误差 `<=1e-12`
  （`test_qf_sg_resolves_a_splitting_that_the_physical_path_rounds_to_zero`）。

**支持范围**：`|u| <= 700`（约 18 V）内 QF 边通量精确（常数 QF 严格零电流）；
超出后 `exp` 溢出，`Fluxes`/`State` 的有限性检查让整条链路显式抛出 `ValueError`，
不会把 `inf`/`nan` 当电流交给求解器
（`test_qf_sg_is_exact_up_to_u_700_and_fails_loudly_beyond`）。SRH 在同一区间内返回有限值；
更极端的 `u` 属于非物理工作点，不作为支持范围声明。

## 6. 冻结的容差与验收门槛

| 项 | 判据 | 证据 |
|---|---|---|
| 接触中性 / 质量作用 | `abs(np-1) <= 1e-10`、`abs(p-n+C) <= 1e-10*max(1,abs(C))` | T03 |
| 极性 / 势垒 | `V>0` 时 `u_R-u_L` 减少 `v`，偏差 `≤ 1e-10` | T04 |
| Bernoulli | 相对误差 `≤ 3e-14`；反射恒等式 `B(-z) = B(z)+z` | T01 |
| SRH 异寿命 | `19/47`，误差 `≤ 1e-12` | T02 |
| 非均匀 Poisson（二次制造解） | 内部残差 `≤ 1e-11` | T05 |
| SG 平衡零电流 | `j_n = j_p = 0`（QF 口径严格 `0.0`） | T06 |
| SG 纯扩散 | `j_n = m_n(n_R-n_L)/h`、`j_p = m_p(p_L-p_R)/h` | T07 |
| 均匀正弦加密 | 相邻误差比落在 `[3.8, 4.2]` | §14.6 P-03 |
| cm↔nd、cm↔SI 往返 | 在声明容差内复原；同器件 cm/SI 输入得同一无量纲状态 | §3.3 |

`Q_floor` 与各量的零阈值按单位分别存入 `convergence.md`；零理论值不做无说明的相对误差除法（§12.1）。

## 7. 非法状态与显式失败（不得静默修复）

| 分类 | 拒绝条件 | 位置 / 错误 |
|---|---|---|
| 网格 | `M < 3`、重复/倒序结点、`h` 与 `diff(x)` 不符、`w` 取错、cm/nd 不成单一尺度、非整数结点数 | `ValueError`（`core/mesh.py`） |
| 分段网格 | 结区坐标不是断点、结区不落在结点上 | `ValueError` 且消息含 `junction` 与坐标（`core/device.py`） |
| 掺杂 | 区域重叠、缺口、未覆盖 `[0,L]`、负 `N_D/N_A` | `ValueError`（`core/region.py`） |
| 材料 | 非正参数、标定温度 ≠ 运行温度 | `MaterialTemperatureMismatch` |
| 接触 | 反向 PN（左 N 右 P） | `INVALID_PN_ORIENTATION`，消息含两端 `C` 值 |
| 接触浓度 | 非有限 `C`；接触浓度超出 float64 可表示范围 | `ValueError`（`core/contact.py`） |
| 右端偏压 | `right_voltage_V ≠ 0` | `ValueError`（`core/device.py`） |
| 单位 | 混用 cm/SI 后缀、无单位字段、非数字 | `UnitError`（`io/units.py`、`io/config.py`） |
| 状态 | 非正 `n/p`、非有限值、`QFState` 形状/参考值非法、非收敛却带 `physics` | `ValueError`/`FloatingPointError` |

## 8. 变更记录

| 日期 | 变更 | 影响 | 提交 |
|---|---|---|---|
| 2026-10-25 | 并入 A 的 QF 残差/SG/SRH/平衡初值，新增 `QFState`/`contact_targets_qf`/`PhysicsDiagnostics` | 内部表示变化，公开 API 不变 | `98ff013` |
| 2026-11-03 | 收紧 W2 校验：`segmented` 结区必须落结点、结点数必须整数、`w` 一致性检查、右端必须 0 V、接触浓度越界拒绝；新增 `docs/numerical_conventions.md`、`docs/convergence_template.md` | 仅拒绝此前被静默接受的非法配置；合法配置结果不变（235 → 253 passed，新增 18 项测试） | 本次提交 |

## 9. 证据索引

```bash
python -m pytest -q                        # 全量
python -m pytest tests/test_oracle_bernoulli.py tests/test_oracle_srh.py \
    tests/test_oracle_poisson.py tests/test_oracle_sg.py \
    tests/test_oracle_qf_residual.py -q    # T01/T02/T05/T06/T07 + QF 集成
PYTHONPATH=src python -m examples.verify_poisson   # SI MMS 与正弦加密误差比
PYTHONPATH=src python -m minitcad.cli validate examples/pn_left_p_right_n.yaml
```
