# Mini-TCAD：一维左 P 右 N 结的数学模型、数值方法与双人 8 周实施验收规范（VS Code 版）

> **项目性质**：**8 周 / 2 人** / Python / 1D Si Poisson–Drift–Diffusion（DD）科研计算软件；**A 在 W1–W7 完成全部物理、离散与数值求解核心任务，W8 独立负责 PySide6 GUI 开发并与 B 接口联调；B 在 W1–W8 完成器件工作流、独立验收、三文件交付与发布任务**。  
> **文档用途**：项目实施、代码审查和验收的统一规范；本文件是**可独立执行的 8 周实施基线**，完整给出数学推导、软件契约、测试阈值、每周任务、Gate、Week 02 教学验证和 W8 GUI/Windows 发布要求；不代表软件已经实现或通过验收。  
> **适用范围**：均匀或非均匀一维网格、非简并 Boltzmann 统计、完全电离、恒温、Si 同质结、常数迁移率、Einstein 关系、SRH 复合、理想欧姆接触、稳态 DC。  
> **核心目标**：不仅能够画 PN 曲线，还必须能用独立的离散恒等式、收敛数据、理论/基准对照证明数值结果可信。第 8 周 GUI 只是同一计算核心的桌面入口，**不得形成第二套物理计算与验收口径**。
>
> **现有源码参考**：随附 `MiniTCAD_1D_PN_PySide6_WindowsEXE_Source_Complete.zip` 已包含 `gui/app.py`、`gui/main_window.py`、`gui/form_model.py`、`gui/worker.py`、`gui/plots.py`、Windows 构建脚本与相关手册，可作为第 8 周开发和接口设计的可复用基线；**是否符合本版 Gate，仍须重新执行验收**。压缩包中的高掺杂压力样例不应被默认视作通过。

> **8 周执行主线**：W1 数学/契约 → W2 Bernoulli/Poisson → W3 SG/SRH/continuity/唯一残差 → W4 Equilibrium/Gummel（Gate 1） → W5 Newton/Jacobian/continuation → W6 双求解器/自适应扫压（Gate 2） → W7 非均匀网格、三文件、复现与核心冻结 → W8 GUI/Windows 发布（Gate 3）。

***

## 0. 必须遵守的决定（先读）

1. **方向固定：从左到右为 P → N**。坐标从左接触 `x=0` 向右接触 `x=L` 增大；左端净掺杂 `C=N_D-N_A<0`，右端 `C>0`。
2. **偏压极性固定：`V_app=V_P(left)-V_N(right)`**；右端为 0 V 参考，`V_app>0` 为**正向**，`V_app<0` 为反向。接触边界必须与左 P、右 N 的方向保持一致。
3. **中能级 SRH 复合公式固定为**：`R=(np-1)/[τ_p(n+1)+τ_n(p+1)]`（所有量已按本文约定归一化）。`τ_p` 乘 `(n+1)`，`τ_n` 乘 `(p+1)`；必须用不相等的电子/空穴寿命进行单元测试。
4. **非均匀 Poisson 分母必须是 `h_{i-1}+h_i`**；相邻间距与当前控制体积必须采用一致的节点索引。
5. **每次 `solve` 固定且仅须生成三件顶层验收文件**：`comparison.csv`（完整数值、有效理论值、误差、迭代序列），`comparison.png`（图形化的数值/理论/误差对照），`convergence.md`（每偏压迭代、容差及失败诊断）。内部缓存与自动单元测试产物不属于交付输出。
6. **负责人 B 必须在用户界面显著提醒方向和极性**：YAML 示例开头、CLI 启动摘要、Python API 文档、PNG 标题/图例和收敛报告均明确写“左 P、右 N；左端施加 `V_app`，右端接地；正向为 `V_app>0`”。PN 示例输入方向错误时直接报错，不能静默反转。
7. **理论值不是万能的**：对于不存在通用闭式解的有偏置 Poisson–DD–SRH 剖面，`theory_value` 留空并标记 `none`；独立细网格/外部求解器数据放入 `reference_value`，不得谎称“精确理论值”。精确零值的**相对误差无定义**，此时报告绝对误差。
8. **物理单位必须显式声明，核心统一无量纲；项目默认使用 cm 工程单位，SI 通过适配层支持**。默认配置/API/GUI 使用 `cm`、`cm^-3`、`cm^2/(V·s)`、`V/cm`、`A/cm^2`、`s`（电压仍为 `V`）；直接以同一套 cm 尺度归一化进入内核，不强制经过 SI。SI 输入经显式适配层转换，若交付需要 SI 文件，仅由输出 writer 转换物理量、列名和图轴，**不得修改无量纲内核、残差、求解器或阈值**。
9. **温度不是默认固定为 300 K**：`T>0` 必须显式提供，`V_T(T)`、`L_0(T)`、`D_0(T)` 等在当前温度计算；`n_i(T)` 必须由明确的材料模型或标明适用温度的用户输入给出。300 K 只用于示例和回归基准，不可拿其 `n_i`/迁移率数值自动套到其他温度。
10. **第 2 周的可测交付固定为**：第 2 周完成 `bernoulli.py`、`poisson.py`、`test_bernoulli.py`、`test_poisson.py` 四项主交付；第 3 周负责 SG、SRH、continuity 和统一 `F(X)` 的集成。
11. **教学 SI Poisson 与生产无量纲 Poisson 明确隔离**：前者可返回内部 `C/m³`、边界 `V` 的混合单位残差，仅用于 Week 02 制造解/教学检查；生产 `residual()` 只能返回第 5、7 节的统一无量纲 residual。两者通过解析缩放与同组制造解交叉校验，不能直接混用范数。
12. **范围与变更控制**：方向/电压、SRH、Poisson/SG、无量纲内核和固定三文件契约为实施基线；教学示例、工程单位换算和开发环境设置不能改变正式物理模型。**W8 PySide6 GUI 是本版明确新增的必做交付**，除此之外不扩大物理模型范围。
13. **八周交付边界**：A 的归一化、Poisson/SG/SRH 离散、continuity/唯一残差、Equilibrium、Gummel、Newton、稀疏 Jacobian、阻尼与正性、强偏置、非均匀网格、性能、派生量及物理终审，必须在 **W7 末**完成并形成可复现核心版本；**W8 A 的主要工作量必须用于 GUI 及与 B 的接口联调**，不得以“核心尚未完成”为由挪走 GUI 周期。
14. **B 的八周交付边界**：B 必须在 W7 末交付并冻结 GUI 所需的配置校验、`run()`/`RunResult`、进度/取消、固定三文件 writer、结果解读契约；W8 负责接口联调、独立端到端验收、Windows 打包和发布，不得把核心 CLI/输出首次开发拖到 W8。
15. **里程碑与不能放松的证据**：Gate 1＝W4（可信平衡态与非零偏压）；Gate 2＝W6（双求解器、Jacobian、自适应扫压）；Gate 3＝W8（完整科学验收、GUI、三文件、可复现发布）。**8 周工期不得通过删除物理校验、失败记录、测试项或放宽验收阈值来实现**。

***

## 1. 项目边界、输入与预期行为

### 1.1 输入与输出

输入：区域及 `N_D/N_A`、Si 材料参数、长度/截面积、温度、网格、欧姆接触、求解器及 DC 扫压范围。计算主链：

```text
YAML / Python API / CLI / PySide6 GUI（A，W8）
         ↓  (B: 显式单位校验；默认 cm → 无量纲；可选 SI 输入适配为 cm；GUI 不自行变换单位)
Device + Mesh + Contact + State
         ↓  (A: Poisson / SG / SRH / 唯一 Residual)
      F(X; V) = 0
         ↓  (A: Gummel 或 Fully-coupled Newton)
Equilibrium → Bias continuation (B，支持 progress / should_stop)
         ↓  (A/B: 独立验证；统一 RunResult)
comparison.csv + comparison.png + convergence.md  (B: 默认 cm；可选 SI 输出适配)
         ↓  (A：GUI 只读取 RunResult / 已验收文件进行展示)
```

求解器物理状态是结点变量 `u_i,n_i,p_i`；网格边上的通量是 `j_n,i+1/2,j_p,i+1/2`。给出电势、电子/空穴、掺杂、空间电荷、SRH、边电场、电子/空穴/总电流密度、终端 `I–V`、按需输出相对能带、两套求解器的收敛记录。

### 1.2 不做与必须做

| 必须独立实现 | 使用成熟库 | 本期不做 |
|---|---|---|
| 无量纲化、欧姆接触、Poisson FV、SG、稳定 Bernoulli、SRH、共享 residual、Gummel、Newton 编排、Jacobian 生成、稀疏组装、线搜索、继续法、验收检查；**W8 PySide6 GUI 与同一计算 API 对接** | NumPy、SciPy 的 `spsolve`/稀疏矩阵、PyYAML、Matplotlib、pytest；**GUI 用 PySide6，Windows 发布用 PyInstaller** | 2D/3D、非 Si 异质结、Fermi–Dirac、复杂温度耦合、Schottky、隧穿、雪崩、量子修正、自写 sparse LU、GUI 内另造物理求解器 |

Python API、CLI、YAML 与 **W8 GUI** 统一调用同一个仿真核心；不存在四套互不一致的实现。GUI 必须显著展示左 P 右 N、正向偏压定义，调用 B 的校验与 `run()`，不能自行计算 PN 方程。第一版的“多 region”主要表示**同一种 Si 的分段掺杂**，并不隐含跨材料通量模型。

***

## 2. 物理量、符号和统一正方向

### 2.1 默认 cm 工程单位与可选 SI 对照符号表

| 符号 | 物理意义 | 项目默认 cm 工程单位 | SI 适配单位 |
|---|---|---|---|
| `x`, `L`, `A` | 位置、器件长度、截面积 | cm、cm、cm²；`x: P(left) → N(right)` | m、m、m² |
| `φ`, `E` | 静电势、电场 | V、V/cm；`E=-∂φ/∂x` | V、V/m |
| `n`, `p`, `N_D`, `N_A`, `n_i` | 电子、空穴、施主、受主及本征浓度 | cm⁻³；`N_D,N_A≥0` | m⁻³ |
| `C=N_D-N_A` | **净施主掺杂** | cm⁻³；P 区负、N 区正 | m⁻³ |
| `J_n,J_p,J=J_n+J_p` | **常规电流方向**的电流密度 | A/cm²；向右为正 | A/m² |
| `μ_n,μ_p`, `D_n,D_p` | 迁移率、扩散系数 | cm²/(V·s)、cm²/s | m²/(V·s)、m²/s |
| `R` | 净 SRH 复合（正为复合） | cm⁻³s⁻¹；本期 `G=0` | m⁻³s⁻¹ |
| `τ_n,τ_p` | SRH 电子、空穴寿命 | s，禁止交换含义 | s |
| `V_app`, `I` | 端电压、终端电流 | V（正向为正）、A | V、A |
| `q,k_B,T,ε` | 元电荷、Boltzmann 常量、温度、Si 绝对介电常数 | C、J/K、K、F/cm | C、J/K、K、F/m |

**接口规则**：此表只规定带单位的输入/输出；核心 `State`、通量和 `residual()` 一律无量纲。配置须通过 `units.system: cm` 或 `units.system: si` 声明输入单位制；默认项目工程契约为 cm。SI 由输入适配层及可选输出 writer 支持，不定义第二套求解内核。

物理假设：稳态、非简并、完全电离、均匀温度、`ε=ε_Si` 常数，`D_{n,p}=μ_{n,p} V_T`。不启用光生载流子（`G=0`）；若未来允许 `G≠0`，必须同步修改两条连续性方程，不能只改一条。

### 2.2 Poisson–漂移扩散–连续性方程

令 `C=N_D−N_A`，物理单位下有：

$$
-\frac{d}{dx}\!\left(\epsilon\frac{d\phi}{dx}\right)
=q\left(p-n+C\right), \qquad E=-\frac{d\phi}{dx}.
$$

$$
J_n=q\mu_n nE+qD_n\frac{dn}{dx}
=q\mu_n\left(V_T\frac{dn}{dx}-n\frac{d\phi}{dx}\right),
$$

$$
J_p=q\mu_p pE-qD_p\frac{dp}{dx}
=-q\mu_p\left(V_T\frac{dp}{dx}+p\frac{d\phi}{dx}\right).
$$

稳态连续性、`R>0` 代表复合时：

$$
\frac{dJ_n}{dx}=qR,\qquad \frac{dJ_p}{dx}=-qR
\quad\Longrightarrow\quad \frac{d(J_n+J_p)}{dx}=0.
$$

**符号核验**：电子项梯度 `dn/dx` 前是正号，空穴项 `dp/dx` 前是负号；电子连续性导数为 `+qR`，空穴为 `-qR`。总电流常数由相加**直接推得**，不是一个可选的可视化检查。

### 2.3 中能级 SRH 复合

一般形式：

$$
R_{\mathrm{SRH}}=\frac{np-n_i^2}{\tau_p(n+n_1)+\tau_n(p+p_1)}.
$$

本期只考虑**中能级陷阱** `n_1=p_1=n_i`，所以：

$$
\boxed{R_{\mathrm{SRH}}=\frac{np-n_i^2}
{\tau_p(n+n_i)+\tau_n(p+n_i)}}.
$$

接下来的无量纲式是：

$$
\boxed{r=\frac{np-1}{\widehat\tau_p(n+1)+\widehat\tau_n(p+1)}}.
$$

此处第二式的 `n,p` 已是无量纲浓度，`\widehat\tau` 是无量纲寿命。**寿命系数必须与对应载流子项匹配。**必须使用 `τ_n≠τ_p` 的测试数据，否则这个错误完全无法检出。

### 2.4 温度、`n_i(T)` 与迁移率参数来源

正式模型是给定器件温度 `T`、材料参数和掺杂后求解，并非默认“所有 Si 固定在 300 K”。定义：

$$
V_T(T)=\frac{k_BT}{q},\qquad
\varepsilon_{\mathrm{Si}}(T)=\varepsilon_{r,\mathrm{Si}}(T)\varepsilon_0,\qquad
D_n(T)=\mu_n(T)V_T(T),\quad D_p(T)=\mu_p(T)V_T(T).
$$

- `q=1.602176634×10^{-19} C`、`k_B=1.380649×10^{-23} J/K` 为基本常数；真空介电常数 `ε_0≈8.8541878128×10^{-14} F/cm`（即 SI 的 `8.8541878128×10^{-12} F/m`）。默认计算 `ε_Si=ε_r ε_0` 时必须采用 `F/cm`，不可将 `F/m` 数值直接代入 cm 公式；参考 `ε_r≈11.7` 仅在本期简化模型中取常量。
- `n_i=n_i(T)` 由用户在**当前温度**显式提供，或通过经评审的具名材料模型计算。300 K 参考参数不足以确定一般 `n_i(T)` 的模型常数；代码**不得**仅用 300 K 的参考数值按温度盲目外推。
- 第一版迁移率 `μ_n,μ_p` 在每次运行中均为常数，但允许对不同温度输入不同的常数；这和“有场依赖/空间依赖的迁移率模型”是两回事，后者仍在本期范围外。
- `material` 的参数来源需可追溯：至少输出 `material_parameter_policy=explicit_at_T` 或已审核的 `named_model`、`temperature_K` 和提供的 `n_i,μ_n,μ_p`。若输入参数明示的标定温度与运行 `T` 不一致，校验阶段拒绝或要求用户提供新参数，**不得自动修复**。
- 源漏附近 `10^19–10^20 cm^-3` 的高掺杂量级仅为工艺参考，并**不**表示本期非简并、完全电离、恒定迁移率假设对任意高掺杂情况都有效。测试基准需声明模型的使用假设。

***

## 3. cm 工程单位到无量纲的完整推导（SI 作为适配）

### 3.1 选取基准尺度

$$
V_T=\frac{k_BT}{q},\quad
L_0=\sqrt{\frac{\epsilon V_T}{q n_i}},\quad
D_0>0\ \text{（配置中明确的参考扩散系数）},\quad
t_0=\frac{L_0^2}{D_0}.
$$

$$
E_0=\frac{V_T}{L_0},\quad
J_0=\frac{qD_0n_i}{L_0},\quad
R_0=\frac{D_0n_i}{L_0^2},\quad
\mu_0=\frac{D_0}{V_T}.
$$

基准尺度采用**相互一致的 cm 单位**：`ε[F/cm]`、`n_i[cm^-3]`、`D_0[cm^2/s]`，由此得到 `L_0[cm]`、`t_0[s]`、`E_0[V/cm]`、`J_0[A/cm^2]`、`R_0[cm^-3/s]`、`μ_0[cm^2/(V·s)]`。绝对介电常数 `ε=ε_r ε_0`；不要把相对介电常数或未经适配的 SI 数值直接代进 `L_0`。所有尺度在同一运行设置下计算一次，并写入 `convergence.md` 供复现。

| cm 工程物理量 | 内核无量纲变量 | 恢复为 cm 工程量 |
|---|---|---|
| `x` | `ξ=x/L_0` | `x=L_0 ξ` |
| `φ`, `V_app` | `u=φ/V_T`，`v=V_app/V_T` | `φ=V_T u` |
| `n,p,N_D,N_A,C` | 各自除以 `n_i`，`Ĉ=(N_D−N_A)/n_i` | 乘 `n_i` |
| `E` | `e=E/E_0` | `E=(V_T/L_0)e` |
| `J_n,J_p` | `j_{n,p}=J_{n,p}/J_0` | `J=J_0 j` |
| `R` | `r=R/R_0` | `R=R_0 r` |
| `μ_n,μ_p` | `m_{n,p}=μ_{n,p}/μ_0=μ_{n,p}V_T/D_0` | `μ=m μ_0` |
| `D_n,D_p` | `d_{n,p}=D_{n,p}/D_0` | `D=D_0 d` |
| `τ_n,τ_p` | `\widehat τ_{n,p}=τ_{n,p}/t_0` | `τ=\widehat τ\,t_0` |

**记号约定**：从本节以下，若未标注具体物理单位，`ξ,u,n,p,C,j_n,j_p,r,m_n,m_p,\widehat τ` 全部表示**无量纲量**，不得把它们直接传给文件输出模块。符号表中带单位的同名物理量须通过 `n_i`、`L_0` 等尺度恢复，不能误把内核值作为 cm 或 SI 值。

### 3.2 代入 Poisson、连续性和输运方程

由 `L_0²=εV_T/(qn_i)`：

$$
\boxed{\frac{d^2u}{d\xi^2}=-(C+p-n)}.
$$

由 `D=μV_T` 与 `J_0=qD_0n_i/L_0`：

$$
\boxed{j_n=m_n\left(\frac{dn}{d\xi}-n\frac{du}{d\xi}\right)},
\qquad
\boxed{j_p=-m_p\left(\frac{dp}{d\xi}+p\frac{du}{d\xi}\right)}.
$$

由 `R_0=D_0n_i/L_0²`：

$$
\boxed{\frac{dj_n}{d\xi}=r,\qquad
\frac{dj_p}{d\xi}=-r},\qquad
r=\frac{np-1}{\widehat\tau_p(n+1)+\widehat\tau_n(p+1)}.
$$

**符号约定**：部分文献采用 `N` 表示净掺杂，但含义可能是净施主或净受主；本规范统一使用 `C=(N_D-N_A)/n_i` 表示无量纲净施主掺杂。

### 3.3 默认 cm → 无量纲；SI 仅通过显式适配层

**不再强制 `cm → SI → 无量纲`**。生产的唯一求解状态保持 `NormalizedDevice + State(u,n,p)`；归一化尺度和工程量先统一在 cm 单位中定义。

```text
YAML / Python API / GUI（声明 units.system: cm；项目默认工程契约）
    → B: Config/物理单位校验
    → cm DeviceConfig、cm Mesh/Material
    → 直接按 L0[cm]、ni[cm^-3]、D0[cm^2/s] 归一化
    → NormalizedDevice + State(u,n,p) → 唯一 F(X)、Gummel/Newton
    → 按 cm 尺度反变换 → 默认 cm CSV / PNG / API 派生量

可选 SI 输入（显式 units.system: si）
    → 输入适配层：SI → cm → 同一条归一化路径（不另造内核）

可选 SI 交付（output.unit_system: si）
    → 求解完成后 writer 将物理输出 cm → SI
    → 同一组三文件、SI 列名和图轴；内核及 RunResult 的无量纲状态不变
```

**cm → 无量纲的直接公式（尺度均取自同一温度与单位制）：**

| 物理量 | 默认 cm 输入或恢复值 | 内核无量纲量 |
|---|---|---|
| 位置/长度 | `x_cm`, `L0_cm` | `ξ=x_cm/L0_cm` |
| 电势/端电压 | `phi_V`, `bias_V`, `VT_V` | `u=phi_V/VT_V`, `v=bias_V/VT_V` |
| 载流子、掺杂、本征浓度 | `n_cm3,p_cm3,ND_cm3,NA_cm3,ni_cm3` | `n_nd=n_cm3/ni_cm3`，`p_nd=p_cm3/ni_cm3`，`Cbar_nd=(ND_cm3−NA_cm3)/ni_cm3`（掺杂先作 CV 平均） |
| 迁移率 | `mobility_cm2_per_Vs` | `m=μ_cm2_per_Vs×VT_V/D0_cm2_s` |
| 扩散系数 | `D_cm2_s` | `d=D_cm2_s/D0_cm2_s` |
| 电场 | `E_V_cm` | `e=E_V_cm×L0_cm/VT_V` |
| 电流密度 | `J_A_cm2` | `j=J_A_cm2/J0_A_cm2` |
| 体复合率 | `R_cm3_s`（值的物理单位为 cm⁻³s⁻¹） | `r=R_cm3_s/R0_cm3_s` |
| 电子/空穴寿命 | `tau_s` | `tau_nd=tau_s/t0_s` |
| 介电常数与面积 | `epsilon_F_cm`、`area_cm2` | `ε` 只用于形成 `L0_cm`；终端 `I_A=area_cm2×Jterminal_A_cm2` |

**必须保留 cm ↔ SI 的双向换算与往返测试**，但它们属于边界适配/验收，不构成默认内核的中间步骤：

| 物理量 | cm → SI | SI → cm |
|---|---|---|
| 长度 | `x_m=1e-2*x_cm` | `x_cm=1e2*x_m` |
| 面积 | `A_m2=1e-4*A_cm2` | `A_cm2=1e4*A_m2` |
| 载流子、掺杂 | `N_m3=1e6*N_cm3` | `N_cm3=1e-6*N_m3` |
| 迁移率 | `mu_m2_Vs=1e-4*mu_cm2_Vs` | `mu_cm2_Vs=1e4*mu_m2_Vs` |
| 扩散系数 | `D_m2_s=1e-4*D_cm2_s` | `D_cm2_s=1e4*D_m2_s` |
| 电场 | `E_V_m=1e2*E_V_cm` | `E_V_cm=1e-2*E_V_m` |
| 电流密度 | `J_A_m2=1e4*J_A_cm2` | `J_A_cm2=1e-4*J_A_m2` |
| 体复合率/数密度变化率 | `R_m3_s=1e6*R_cm3_s` | `R_cm3_s=1e-6*R_m3_s` |
| 体电荷密度 | `rho_C_m3=1e6*rho_C_cm3` | `rho_C_cm3=1e-6*rho_C_m3` |
| 绝对介电常数 | `eps_F_m=1e2*eps_F_cm` | `eps_F_cm=1e-2*eps_F_m` |
| 电势、电压、电流、寿命、温度 | 数值不变（V、A、s、K） | 数值不变（V、A、s、K） |

单元测试必须同时覆盖 `cm→SI→cm`、`SI→cm→SI`、`cm→无量纲→cm` 及 `SI→适配 cm→无量纲→恢复 cm→SI`，在声明的浮点容差内恢复原值；另以相同物理器件验证 cm 输入和 SI 输入产生一致的无量纲状态与 `F(X)`。对于 `output.unit_system: si`，测试**仅输出层**物理量、位置列名、`quantity`/`unit`、理论/参考值、绝对误差及图轴的换算；无量纲残差、相对误差、迭代结果不随输出单位改变。无单位的 `density=1e16` 一律拒绝或要求显式单位，绝不按数值大小猜测。

独立 SI 教学 Poisson 的 `rho_C_m3` 单位仍为 `C/m³`，不是 `m^-3`；该教学验证器不参与生产归一化路径。

### 3.4 两种 Debye 长度必须区别

参考长度 `L_0` 以当前温度下的本征浓度 `n_i(T)` 为浓度基准：

$$
L_0(T)=L_{D,n_i}(T)=\sqrt{\frac{\varepsilon_{\mathrm{Si}}V_T(T)}{q\,n_i(T)}}.
$$

局部屏蔽长度的数量级可用当前实际浓度 `N` 定义：

$$
L_{D,\mathrm{local}}(N,T)
=\sqrt{\frac{\varepsilon_{\mathrm{Si}}V_T(T)}{q\,N}}.
$$

两者**不是** PN 耗尽层宽度。以 300 K、`ε_r=11.7`、`n_i=10^10 cm^-3` 作为演示参数时，`L_0≈40.88 μm`；对于 `N=10^16 cm^-3`，局部长度约 `40.88 nm`，相差三个数量级。真实结附近能否分辨屏蔽区，要看**物理网格 `h_cm` 与同为 cm 的局部物理尺度之比**，不能因为无量纲总长度看起来很小就误以为网格足够细。

### 3.5 `D_0` 为参考选择，不是“固定的 Si 材料常数”

参考扩散系数 `D_0>0` 必须通过配置明确确定。另一种可选的参数化方法是：

$$
\mu_0=1000\ \mathrm{cm^2/(V\,s)}
=0.1\ \mathrm{m^2/(V\,s)},\qquad
D_0(T)=\mu_0 V_T(T).
$$

这只是**备选参数化**：使用时明确 `reference_mobility_cm2_per_Vs=1000.0`，据当前 `T` 计算 `D_0(T)=1000 V_T(T)`，单位 `cm²/s`；本规范第 13 节示例选择 `reference_diffusivity_cm2_s=10.0`，代表当次运行固定给定的 `D_0`，不可悄悄改成 `1000 V_T`。`1000 cm²/(V·s)=0.1 m²/(V·s)`、`10 cm²/s=1.0e-3 m²/s` 只用于 SI 适配测试。两种基准只能选其一；报告记录生效的 `D_0[cm²/s]` 和 `μ_0[cm²/(V·s)]`。

### 3.6 SI 教学 Poisson 与正式无量纲 Poisson 的精确关系

为核对独立 SI Poisson 验证器与第 5 节正式无量纲残差的符号和尺度，定义（常数 `ε=ε_Si`）：

$$
F_{P,i}^{SI}
=-\frac{[\epsilon(\phi_{i+1}-\phi_i)/h_i^{SI}]
       -[\epsilon(\phi_i-\phi_{i-1})/h_{i-1}^{SI}]}
       {(h_{i-1}^{SI}+h_i^{SI})/2}
-q\,n_{\mathrm{int}}(T)(\bar C_i+p_i-n_{e,i}),
$$

此式用 `n_int(T)` 特指**本征浓度**，`n_{e,i}` 特指**第 i 个节点的无量纲电子浓度**，避免两者同写成 `n_i` 的歧义；**SI 教学 oracle** 使用 `ni_material_m3`，生产 cm 配置使用 `ni_material_cm3`，内核节点使用 `n_node_nd`。教学侧由 SI 适配获得相同物理设备的 `x_m`、`ε_F_m`，代入 `x=L_0ξ`、`φ=V_Tu` 后：

$$
\boxed{F_{P,i}^{SI}=-q\,n_{\mathrm{int}}(T)\,F_{u,i}^{nd}},
\quad i=1,\dots,M-2.
$$

边界行则有 `F_{P,0}^{SI}=V_T F_{u,0}^{nd}`、`F_{P,M-1}^{SI}=V_T F_{u,M-1}^{nd}`。因此**内部行要除以 `-q n_i`，边界行要除以 `V_T`**，不能用同一个系数盲目缩放所有行。

此等价关系只在两套程序使用相同的节点、控制体积、绝对介电常数、体电荷、Dirichlet 和常数 Si 假设时成立。若在 cm 物理单位下构造独立 Poisson 验证行，相应关系是 `F_P^cm=−q n_i[cm^-3] F_u^nd`（内部单位为 `C/cm³`），边界仍为 `V_T F_u^nd`；它与 SI 教学行经体积单位适配后等价。若教学样例启用节点变化的 `epsilon` 调和平均，必须作为**额外的 SI 离散实验**单独报告，而不是假装仍属于本期恒定 `ε_Si` 的生产模型。

***

## 4. 左 P 右 N 欧姆接触与偏压条件

### 4.1 由局部电中性 + 质量作用定律求接触浓度

理想欧姆接触处满足：

$$
C+p-n=0,\qquad np=1
\quad\Rightarrow\quad n-p=C.
$$

将 `p=1/n` 代入：`n²−Cn−1=0`，物理解为：

$$
n_{\mathrm{ct}}(C)=\frac{C+\sqrt{C^2+4}}2,\quad
p_{\mathrm{ct}}(C)=\frac{-C+\sqrt{C^2+4}}2=\frac1{n_{\mathrm{ct}}(C)}.
$$

**稳定实现**（在高度 P 掺杂时避免相减消失）：

$$
a(C)=\operatorname{asinh}(C/2),\quad
n_{\mathrm{ct}}=e^{a(C)},\quad p_{\mathrm{ct}}=e^{-a(C)},\quad u_{\mathrm{ct}}=a(C).
$$

这是选取热平衡共同准费米能级为零时的规范；整体电势允许加常数，但**本项目固定右端参考和上述热平衡规范**。

### 4.2 接触 Dirichlet 条件及正向定义

令 `C_L<0`（P）、`C_R>0`（N），且右金属端接地：

$$
\boxed{u(0)=a(C_L)+v,\qquad u(\ell)=a(C_R)},
\quad v=\frac{V_{\mathrm{app}}}{V_T},\quad \ell=L/L_0.
$$

$$
\boxed{n(0)=e^{a(C_L)},\ p(0)=e^{-a(C_L)},\ 
 n(\ell)=e^{a(C_R)},\ p(\ell)=e^{-a(C_R)}}.
$$

这些条件表达“接触浓度按当地掺杂保持局部平衡，左侧金属整体上移 `V_app`”。**右金属端外加电压为零并不表示右侧半导体内部的静电势必须设成 `u(\ell)=0`**；在所选本征能级参考下，右接触仍有掺杂引起的电势偏移 `a(C_R)`。正向 `V_app>0` 时，势差 `u_R-u_L` 变小，内建势垒降低；反向则升高。**禁止**将左 N、右 P 的接触公式或相反的参考电极约定用于本项目；否则偏压极性与边界势会不一致。

无偏压内建电势差与其高掺杂近似为：

$$
V_{\mathrm{bi}}=V_T\big[a(C_R)-a(C_L)\big]
\ \approx\ V_T\ln\frac{N_A^{(L)}N_D^{(R)}}{n_i^2}
\quad (N_A^{(L)},N_D^{(R)}\gg n_i).
$$

右端相对左端的**内部静电势差** `+V_bi` 不等于用户指定的端电压 `V_app`：热平衡下外部端电压为零，而内部仍存在内建电势。不可混淆。

### 4.3 B 必须实现的用户提示与校验

- CLI 运行首屏：`PN orientation: P(left, x=0) -> N(right, x=L) | V_app = V_left-V_right | +V = forward bias | right = 0 V`，中文帮助同步说明。
- 在 YAML 示例最上方写该方向；Python API 文档/异常、PNG 标题/图例、收敛报告开头重复（简明即可）。
- 仅 `device.type: pn` 的严格基准要求 `C_L<0<C_R`，违反即报 `INVALID_PN_ORIENTATION` 并标出两端 `C` 的单位和实际值。泛用 `device.type: generic_1d` 允许同种 Si 的任意掺杂配置，但**不能宣称一定是左 P 右 N 的 PN 验收实例**。
- `V_app>0` 的端电压必须映射到左端 `u(0)+=V_app/V_T`，右端接触不随电压改变；不能在 `solve()` 或 `sweep()` 内重新解释符号。

***

## 5. 非均匀网格与 Poisson 有限体积离散

### 5.1 网格定义及节点控制体

采用 **0-based** 索引，`i=0,...,M−1`：

$$
0=\xi_0<\xi_1<\cdots<\xi_{M-1}=\ell,\quad
h_i=\xi_{i+1}-\xi_i>0\ (i=0,\ldots,M-2).
$$

内部节点的控制体宽度

$$
w_i=\frac{h_{i-1}+h_i}{2},\quad 1\le i\le M-2,
$$

对应物理宽度 `L_0 w_i`。左/右接触节点直接施加 Dirichlet，不组装内部 PDE。对于分段掺杂，应由**控制体内的分段积分**得到 `\bar C_i=(1/w_i)\int_{CV_i}C(\xi)\,d\xi`；`ND_cm3`、`NA_cm3` 在内部结点导出控制体平均浓度，接触结点导出接触侧的局部浓度，并用 `notes` 标明。不能把恰处阶跃点的一侧值当作整个控制体的掺杂。

### 5.2 从通量差推导 Poisson 离散

定义单元两侧势梯度：

$$
g_{i+1/2}=\frac{u_{i+1}-u_i}{h_i},\qquad
g_{i-1/2}=\frac{u_i-u_{i-1}}{h_{i-1}}.
$$

在节点控制体上积分 `u''=-(C+p-n)`，以节点值表示连续项：

$$
\boxed{
\frac{2}{h_{i-1}+h_i}
\left[\frac{u_{i+1}-u_i}{h_i}
-\frac{u_i-u_{i-1}}{h_{i-1}}\right]
=-\left(\bar C_i+p_i-n_i\right)}.
$$

因此定义 Poisson residual：

$$
\boxed{F_{u,i}=\frac{g_{i+1/2}-g_{i-1/2}}{w_i}
+\bar C_i+p_i-n_i}.
$$

**索引检查**：内部节点 `i` 的两侧间距必须为 `h_{i-1}` 和 `h_i`，因此控制体积宽度为 `(h_{i-1}+h_i)/2`。若错取下一条边的间距，非均匀网格将无法正确复现二次函数的二阶导数。只有两侧间距相同为 `h` 时，才退化为 `(u_{i+1}-2u_i+u_{i-1})/h²`。

### 5.3 网格模块的最低要求

- `uniform(M)`、`segmented(breakpoints, per_segment_nodes)` 两种确定性生成器；`M≥3`，端点精确覆盖 `[0,L]`；结区坐标必须包含在网格断点中。
- 默认工程网格 `x_nodes_cm` 严格递增，`h_edges_cm` 长度 `M-1`，`w_interior_cm` 长度 `M-2`；归一化网格为 `x_nodes_nd/h_edges_nd/w_interior_nd`。SI 网格如存在，仅在适配边界以 `x_nodes_m` 表达，不参与求解。`mesh` 与 `device` 的结构不允许被 Newton 迭代修改。
- `region` 可表达任意段 Si 掺杂；重叠、空段、负掺杂浓度、未覆盖区域、非递增网格直接在 B 的配置验证阶段报错。

### 5.4 Week 02 Poisson 教学检查必须覆盖的四层证据

独立 **SI** 验证器使用 `-d(ε dφ/dx)/dx=ρ`：节点 `φ[V]`、`ρ[C/m³]`，边通量 `εΔφ/Δx [C/m²]`，内部控制体积宽度 `[m]`。它的内部残差取 `-(flux_right-flux_left)/w-rho`，与第 5.2 节无量纲 `+(gradient_right-gradient_left)/w+C+p-n` **只差第 3.6 节明确的负尺度因子**，不是物理方向不同。

1. **直线制造解**：`φ=a+b x/L`，设置 `ρ=0`，两端精确 Dirichlet；用于检查零源项、索引、两端符号。
2. **二次制造解**：`φ=a+b x/L+c(x/L)^2`，`ρ=-2ε c/L²`；在常 `ε` 下，无论均匀或严格递增的非均匀网格，离散二阶导数都应精确到舍入误差。独立测试必须使用真实的相邻 `h_left/h_right`。
3. **正弦制造解（MMS）**：`φ=sin(πx/L)`、`ρ=επ² sin(πx/L)/L²`、两端为零；在常 `ε`、**均匀网格**解线性系统，用 `N=17,33,65,129` 比较最大电势误差和相邻误差比，二阶方法目标误差比在 `[3.8,4.2]`。不可把二次函数“恰好精确”误称为一般二阶证明。
4. **故障注入**：重复/倒序节点、非有限边界、`ε≤0`、数组形状不匹配、把 `ε_r=11.7` 错当绝对 `ε[F/m]`、只用一个固定 `dx`、Dirichlet 行误作 PDE 行；每个故障都应触发独立错误或相应回归测试。

**可选扩展**：节点 `epsilon` 调和平均是独立的验证实验；本期 Si 同质结的生产 `poisson.py` 使用 `ε_Si` 常量，是否支持跨材料/变量 `ε` 不属于 W2 或本期 Gate 的强制范围。

***

## 6. Scharfetter–Gummel（SG）载流子通量的推导

**目标**：在区间 `[\xi_i,\xi_{i+1}]` 内令无量纲电势线性、迁移率和电流近似常数，把漂移与扩散作为一个整体积分；不直接用中心差分分别离散两项。

定义 `\Delta u_i=u_{i+1}-u_i`、`h_i=\xi_{i+1}-\xi_i`，以及稳定的 Bernoulli 函数：

$$
B(z)=\begin{cases}\frac{z}{e^z-1},& z\ne0,\\1,&z=0.\end{cases}
\quad B(-z)=e^zB(z),\quad B(z)-B(-z)=-z.
$$

### 6.1 电子通量的积分

从 `j_n=m_n(n'-nu')` 出发：

$$
\frac{d(ne^{-u})}{d\xi}=\frac{j_n}{m_n}e^{-u},\quad
n_{i+1}e^{-u_{i+1}}-n_ie^{-u_i}
=\frac{j_{n,i+1/2}}{m_{n,i+1/2}}
\int_{\xi_i}^{\xi_{i+1}}e^{-u(\xi)}d\xi.
$$

用线性 `u(\xi)=u_i+(\Delta u_i/h_i)(\xi-\xi_i)` 积分：

$$
\int_{\xi_i}^{\xi_{i+1}}e^{-u(\xi)}d\xi
=h_ie^{-u_i}\frac{1-e^{-\Delta u_i}}{\Delta u_i}.
$$

整理可得电子边电流：

$$
\boxed{j_{n,i+1/2}=\frac{m_{n,i+1/2}}{h_i}
\left[n_{i+1}B(\Delta u_i)-n_iB(-\Delta u_i)\right]}.
$$

**实现约定**：电子通量统一采用上述 Bernoulli 形式，边界索引为 `i+1/2`，不得引入含义不明的中间位置下标。

### 6.2 空穴通量的积分

由 `j_p=-m_p(p'+pu')`：

$$
\frac{d(pe^u)}{d\xi}=-\frac{j_p}{m_p}e^u,
$$

在同一个区间线性插值 `u` 并积分，得到：

$$
\boxed{j_{p,i+1/2}=\frac{m_{p,i+1/2}}{h_i}
\left[p_iB(\Delta u_i)-p_{i+1}B(-\Delta u_i)\right]}.
$$

`m_n,m_p` 第一版常数；未来允许位置相关迁移率时，需正式定义边迁移率 `m_{i+1/2}` 的平均规则并增加一致性测试。

### 6.3 立即可检验的四条恒等式

1. `\Delta u=0` 时，`j_n=m_n(n_{i+1}-n_i)/h`、`j_p=-m_p(p_{i+1}-p_i)/h`（纯扩散极限）。
2. 结点热平衡 `n_{i+1}=n_ie^{\Delta u}`、`p_{i+1}=p_ie^{-\Delta u}` 时，**两条 SG 通量都为零**。
3. `n,p` 常数但 `\Delta u\to0` 时，电子与空穴通量都趋向由 `-m n u'`、`-m p u'` 给出的漂移电流。
4. 交换区间方向必须同时交换端点与电流取向，不能只交换 `B(z)` 与 `B(-z)`。

### 6.4 避免浮点消减的 Bernoulli 实现

$$
B(z)=1-\frac z2+\frac{z^2}{12}-\frac{z^4}{720}
+\frac{z^6}{30240}+O(z^8).
$$

- `|z|<10^{-3}`：统一使用上述至 `z^6` 的级数；`10^{-3}≤|z|≤50` 使用 `z/expm1(z)`。该阈值是统一的数值实现约定，不是数学上唯一有效的分段选择。
- 大正 `z`：使用 `z exp(-z)/(1-exp(-z))`，避免 `exp(z)` 溢出；大负 `z`：使用 `-z/(1-exp(z))`。
- 四支路边界固定为 `|z|<1e-3`、`z>50`、`z<−50` 与其余中间区。单测覆盖 `0`、`±10^-14`、`±10^-12`、`±10^-3` 两侧、`±1`、`±50` 两侧、`±700`、`±1000`、空数组、二维数组和非有限输入；有限输入的结果必须有限且非负（极大正数正确下溢为 `0.0` 是允许的），不可通过任意 clipping 改变物理通量。
- 若选 complex-step Jacobian，分段逻辑不得对复数直接调用非解析的 `abs` 形成错误梯度；先由真实基准点决定分段，再在该解析支路内求导，或改用 AD/解析局部导数。

### 6.5 第 2 周 Bernoulli 函数级契约

`bernoulli(z)` 允许有限实数标量、列表与 NumPy 数组；标量返回 Python `float`，向量/矩阵返回同形 `float64 ndarray`，空数组返回空数组。含 `NaN/±Inf` 必须抛 `ValueError`；输入不得被原地改写。实现时对每个布尔掩码对应的切片**单独计算**，不要在 `np.where(condition, safe, overflow_expr)` 中计算危险分支。

独立自检优先检查 `B(−z)=B(z)+z`、`B(0)=1` 和 80 位 `decimal` 参照；正常可表示量级的 Week02 精度目标为约 `3e−14` 相对误差，极大正 `z` 下溢则按浮点规则比较。T01 的 Gate 阈值为 `≤1e−11`；**第 2 周自测应先达到更严格的精度目标**，如有无法表示的真值要报告下溢原因，不用“0%”掩盖误差。

***

## 7. 连续性残差、守恒与边界 residual

在每个内部节点采用**同一条边通量**，不得由电子和空穴模块各自重新计算方向或插值。以 `r_i` 表示该节点无量纲 SRH 值：

$$
\boxed{F_{n,i}=\frac{j_{n,i+1/2}-j_{n,i-1/2}}{w_i}-r_i=0},
$$

$$
\boxed{F_{p,i}=\frac{j_{p,i+1/2}-j_{p,i-1/2}}{w_i}+r_i=0}.
$$

相加得到**离散严格恒等式**：

$$
F_{n,i}+F_{p,i}=\frac{j_{i+1/2}-j_{i-1/2}}{w_i},
\quad j=j_n+j_p.
$$

因此收敛后总边电流在整个器件上近似常数；如果关闭 SRH（`r=0`），电子和空穴各自也是常数。这同时规定了代码的 `current_conservation` 检验对象是**边电流**（长度 `M−1`），而不是误把接触结点电流作为已计算的独立结点量。

### 7.1 统一的残差布局

采用结点交错排序：

$$
X=[u_0,n_0,p_0,u_1,n_1,p_1,\ldots,u_{M-1},n_{M-1},p_{M-1}]^T\in\mathbb R^{3M}.
$$

内部每个节点的 `(F_u,F_n,F_p)` 由上述三式定义；两端 6 条 PDE 行被**Dirichlet residual** 代替，例如 `F_{u,0}=u_0-u_{L,target}`、`F_{n,0}=n_0-n_{L,target}`、`F_{p,0}=p_0-p_{L,target}`，右端同理。两种求解器必须调用**同一套** `residual(X, device, bias)`、`edge_flux(state, device)`、`contact_targets(device,bias)`。

### 7.2 状态限制与额外导出

- 任何被接受的状态都必须 `n_i>0`、`p_i>0`、各量有限；阻尼线搜索检查候选状态，失败应回退并记载，不可在迭代后无提示裁剪。
- 边电场：`E_{i+1/2}=-(V_T/L_0)(u_{i+1}-u_i)/h_i`（默认 `V/cm`）；边 `J_{n,p}=J_0j_{n,p}`（默认 `A/cm²`）。约定终端 `J_{\mathrm{terminal}}=J_{1/2}`（左端第一条网格边，取向为向右），`I[A]=A[cm²]·J_{\mathrm{terminal}}[A/cm²]`；右端边电流和沿线中位数只用于守恒交叉检查。`A` 必须输入，否则只能报告 `J(V)`、不能虚构 `I(V)`；SI 展示仅由输出适配器执行。
- 能带为可选**相对能带**：在同一 Si 中，取 `E_c(x)-E_c(0)=-[\phi(x)-\phi(0)]` eV、`E_v=E_c-E_g`（共同能量零点任意）；图例必须标“relative band energy (eV)”。

***
## 8. Gummel、Fully Coupled Newton 与 Jacobian

### 8.1 一次“求解”的唯一数学定义

对给定 `Device, Mesh, V_app`，找到满足六条接触条件及所有内部残差的 `X`：

$$
\boxed{F(X;V_{\mathrm{app}})=0},\qquad F:\mathbb R^{3M}\to\mathbb R^{3M}.
$$

`F` 是**物理真相的唯一代码来源**。禁止在 Gummel、Newton、Jacobian 与输出模块里各维护一个稍有不同的复合/接触/通量公式。子求解器可以有不同的迭代策略，但**最终成功状态必须通过完整的 `F` 复核**。

### 8.2 Gummel：分裂更新与严格退出条件

```text
以 equilibrium/上一电压解为 X^(0)
while iteration < max_iter:
    解 F_u(u | n,p 固定)=0                 → u_new
    解 F_n(n | u_new,p 固定)=0             → n_new
    解 F_p(p | u_new,n_new 固定)=0         → p_new
    对更新执行阻尼和正性检查
    构造完整 F(u_new,n_new,p_new)
    同时检查 scaled residual 与 scaled step
    若均满足：返回 converged；否则继续
超过迭代/无法接受步长：返回显式失败，附原因
```

电子/空穴方程因 SRH 分母包含其自身浓度，严格说子问题仍可能**非线性**。实现时可使用局部一维/稀疏子系统 Newton 或对 `R` 做滞后迭代；无论采用哪种，外层必须使用**未冻结的完整 residual** 判定收敛。Gummel 是基线，不代表可以使用不同物理方程。

### 8.3 Newton：稀疏耦合求解

$$
\mathcal J(X_k)=\frac{\partial F(X_k)}{\partial X},\qquad
\mathcal J(X_k)\Delta X_k=-F(X_k),
$$

$$
\boxed{X_{k+1}=X_k+\alpha_k\Delta X_k},\qquad
\alpha_k\in(0,1].
$$

每个内部结点的 residual 仅依赖其本身及相邻两结点；在交错排序下，`\mathcal J` 是**3×3 块三对角稀疏矩阵**（接触行经过 Dirichlet 替换）。组装使用 `scipy.sparse`；线性求解使用 `scipy.sparse.linalg.spsolve`，不自写 LU。求解失败、`NaN/Inf`、奇异/病态警报或线搜索达到最小步长，必须返回结构化失败信息，不得打印“converged”。

实现优先顺序：**可靠残差 → 稀疏数值 Jacobian → 独立 FD 检查 → 阻尼/收敛诊断**。数值差分的基础校验式为：

$$
\mathcal J_{:,j}^{\rm FD}
\approx\frac{F(X+\delta_je_j)-F(X-\delta_je_j)}{2\delta_j}.
$$

`\delta_j` 应按变量尺度确定，涉及 `n,p` 时必须保证试探状态为正，否则使用适合约束的单边差分/安全扰动；**不能**对所有变量一律使用同一绝对差分幅度。正式 Jacobian 采用邻接稀疏模式（可用颜色分组减少 residual 次数），B 的参考 checker 在小网格上独立计算，不复用 A 的稀疏导数代码。

### 8.4 阻尼和状态合法性

- 初始 `α=1`，不接受则依次减半，直至候选所有 `n,p>0`、数值有限且归一化残差足够下降；`α<α_min` 判为 `LINE_SEARCH_FAILED`。
- 允许设置 `max_iter`、`α_min`、`tol_residual`、`tol_update`、`linear_solver`、`jacobian_mode`；每次尝试的 `α`、残差、更新量均记录。
- Jacobian 的连续可微性由 SRH 和 Bernoulli 计算实现保障；不能对负载流子简单取绝对值修补。

### 8.5 双收敛判据的可执行定义

收敛必须同时满足残差条件与相对更新条件。由于变量量级不同，直接使用 `||F||_2/M` 或 `Δx/x` 可能误判；统一采用以下经过缩放的定义：

$$
\mathcal R_k=\max_j\frac{|F_j(X_k)|}{s_{F,j}(X_k)},\qquad
\mathcal D_k=\max_j\frac{|X_{k,j}-X_{k-1,j}|}{s_{X,j}(X_k,X_{k-1})}.
$$

行尺度（仅用于监控，**不更改物理残差**）约定：

- Poisson 内部 `s_{F,u,i}=max(1,|\bar C_i|,n_i,p_i)`；电子内行 `s_{F,n,i}=max(1,|j_{n,i±1/2}|/w_i,|r_i|)`；空穴内行同理；Dirichlet 行 `s=max(1,|target|)`。
- 更新尺度：势 `s_{X,u}=max(1,|u_k|,|u_{k-1}|)`；浓度 `s_{X,n}=max(|n_k|,|n_{k-1}|,n_floor)`、`s_{X,p}` 同理。`n_floor>0` 只是**归一化分母保护参数**，不是将真实浓度硬设为该值。
- 默认建议阈值（首周集成后可用固定 benchmark 标定并记录版本）：`tol_residual=1e-8`、`tol_update=1e-6`、`max_iter=60`、`n_floor=1e-15`。
- `converged` 的必要条件是**两者同时**满足；第 0 步已满足时 `\mathcal D_0:=0`。还须运行 11 节的正性、边界和电流守恒独立检查；判据被放松必须记录于报告，不能默默更换。
- 物理自洽验证不是 Newton 停止条件的替代品；如果 `\mathcal R` 很小但电流不守恒，应当判为验收失败并调查 residual/离散实现。

***

## 9. Equilibrium 初始化与自适应 DC 扫压

### 9.1 无偏压平衡态的更可靠初始化

采用共同平衡准费米参考 `E_F=0`，由 Boltzmann 关系：

$$
n(\xi)=e^{u(\xi)},\qquad p(\xi)=e^{-u(\xi)},\qquad np=1.
$$

将其带入 Poisson，解**单一非线性标量方程**：

$$
\boxed{u''(\xi)=-\big(C(\xi)+e^{-u(\xi)}-e^{u(\xi)}\big)},
$$

边界 `u(0)=a(C_L)`、`u(\ell)=a(C_R)`。收敛后组装 `n=e^u,p=e^{-u}` 并交给完整残差检查：`j_n≈0,j_p≈0,np≈1`。此步骤优先于直接拿阶跃掺杂分布作 Newton 初值，可减轻初始解不连续带来的不稳定。

### 9.2 自动继续法与失败回滚

```text
equilibrium success → 保存 V=0 基准
当前成功点 V_k, X_k
  ├─ 尝试 V_trial=V_k+ΔV，使用 X_k 作为初值并更新接触 BC
  ├─ 成功：记录状态/迭代/步长；可放大下一步 ΔV（上限 max_step）
  └─ 失败：保留 X_k 不变；ΔV ← shrink_factor × ΔV 并重试
       当 ΔV < min_step 或重试超限：标记失败、保存原因、停止该方向
```

`ΔV` 是**端电压**的步长，允许正向/反向分别按有符号目标处理；终点需精确命中。不得把失败迭代的临时解作为下一点的初值。`sweep` 与 `solver` 的职责分开：B 控制偏压和回滚，A 返回稳定的 `SolveResult`。

对比实验要统一网格/初值/容差，分别报告 fixed sweep 与 adaptive sweep 的成功率、总迭代数和运行时间，不能只比较最终 IV 图是否相似。

***

## 10. 软件模块、函数签名、数据契约与 Owner

### 10.1 约定的源代码结构

下列以 **`src/` 布局**表示正式 Python 包；`tests/`、`examples/`、`docs/`、`reports/` 放在**仓库根目录**，而非包内。GUI 目录名与随附源码包中的 PySide6 参考实现对齐。已有非 `src/` 仓库可采用等价布局，但不可形成两套独立 `Mesh/State/Residual`。

```text
src/minitcad/
├── core/
│   ├── material.py           # B: cm-first material + explicit units validation
│   ├── region.py             # B: Piecewise doping
│   ├── mesh.py               # B: uniform/segmented mesh + control volumes
│   ├── device.py             # B: DeviceConfig / normalized device
│   ├── contact.py            # B: left-P/right-N contact targets and hints
│   └── state.py              # B: immutable mapping, State / SolveResult
├── physics/
│   ├── normalization.py      # A: cm <-> dimensionless formulas (one core)
│   └── recombination.py      # A: midgap SRH
├── discretization/
│   ├── bernoulli.py          # A: stable B(x)
│   ├── poisson.py            # A: FV Poisson and CV doping interface
│   ├── sg.py                 # A: electron/hole edge flux
│   ├── continuity.py         # A: two continuity residuals
│   └── residual.py           # A: unique F(X, device, bias)
├── solvers/
│   ├── gummel.py             # A
│   ├── newton.py             # A
│   ├── jacobian.py           # A: sparse Jacobian generation
│   ├── linear.py             # A: call SciPy sparse solve
│   ├── damping.py            # A
│   └── convergence.py        # A: shared two-stop criteria
├── simulation/
│   ├── equilibrium.py        # A owns math; B owns calling workflow
│   ├── dc.py                 # B: single-bias orchestration
│   └── sweep.py              # B: adaptive/fixed continuation + rollback
├── verification/
│   ├── fd_jacobian_oracle.py # B: independent checker
│   ├── exact_identities.py   # B: contact, equilibrium, conservation
│   ├── reference.py          # B: refined-mesh/reference comparisons
│   └── acceptance.py         # B: Gate tests and three-file audit
├── io/
│   ├── config.py             # B: YAML schema + declared cm/SI validation
│   ├── units.py              # B: SI input adapter / optional SI output conversions
│   ├── results_csv.py        # B: one long-form CSV
│   ├── plot.py               # B: one PNG multipanel chart
│   └── report.py             # B: one convergence Markdown
├── gui/                      # A: W8 GUI Owner；只消费 B 的公共 API
│   ├── app.py                # A: PySide6 入口、--smoke-test
│   ├── main_window.py        # A: 左参数 / 右 I–V、剖面与报告
│   ├── form_model.py         # A: YAML / 表单纯数据适配；B 校验
│   ├── worker.py             # A: QThread、progress、协作式取消
│   └── plots.py              # A: 只绘制实际接受状态
├── api.py                    # B: public Python interface，W7 冻结
├── cli.py                    # B: minitcad solve / validate
└── __init__.py               # python package
apps/minitcad_gui.py            # repo root，A: Windows/desktop 启动入口
scripts/build_windows.bat       # repo root，B: W8 Windows onedir 构建
scripts/build_windows_onefile.bat # 可选单文件构建
docs/GUI_USER_GUIDE.md          # A 撰写 / B 独立验证
docs/WINDOWS_BUILD.md           # B: Windows 构建与烟雾测试
examples/pn_left_p_right_n.yaml  # repo root
examples/verify_poisson.py       # repo root, Week02 standalone SI MMS
tests/                           # repo root, A tests + B independent oracles
docs/                            # repo root, Week01 conventions / Week02 AB review
reports/                         # repo root, internal weekly evidence, NOT solve outputs
```

**协同规则**：A 对方程、符号、归一化和求解器正确性最终负责，**W8 对 GUI 代码及交互行为负责**；B 对 `Device/API/CLI/输出/可复现实验/独立验收` 最终负责，**W8 对公共接口、端到端集成和 Windows 发布验证负责**。两人必须共同评审**端电压极性、Residual 索引和 CSV 字段表**；A 的提交不能仅凭自己的测试通过，B 的集成不能在出现物理测试失败时标为成功。

### 10.2 不可临时改变的最小接口

```python
# 以下为签名契约，不是省略实现的“可运行库”。
@dataclass(frozen=True)
class Mesh1D:
    x_nodes_cm: np.ndarray       # shape (M,), increasing; physical cm grid [0,L_cm]
    h_edges_cm: np.ndarray       # shape (M-1,), cm
    w_interior_cm: np.ndarray    # shape (M-2,), cm
    x_nodes_nd: np.ndarray       # shape (M,), dimensionless
    h_edges_nd: np.ndarray       # shape (M-1,), dimensionless
    w_interior_nd: np.ndarray    # shape (M-2,), dimensionless

@dataclass(frozen=True)
class State:
    u: np.ndarray                # (M,), normalized potential
    n: np.ndarray                # (M,), normalized electron concentration > 0
    p: np.ndarray                # (M,), normalized hole concentration > 0

@dataclass(frozen=True)
class Fluxes:
    jn_edges: np.ndarray         # (M-1,), normalized conventional current
    jp_edges: np.ndarray         # (M-1,), normalized conventional current

@dataclass(frozen=True)
class SolveResult:
    bias_V: float
    state: State | None
    converged: bool
    iterations: int
    residual_norm: float
    update_norm: float
    reason: str                  # SUCCESS / MAX_ITER / ...
    iteration_history: tuple    # each record: step, alpha, residual_norm, update_norm

# A 提供给 B 的数值入口：
def residual(state: State, device: "NormalizedDevice", bias_V: float) -> np.ndarray: ...
def edge_flux(state: State, device: "NormalizedDevice") -> Fluxes: ...
def solve_gummel(device, bias_V: float, initial: State, settings) -> SolveResult: ...
def solve_newton(device, bias_V: float, initial: State, settings) -> SolveResult: ...

# B 提供给 A 的器件入口：
def parse_config(path: str) -> "DeviceConfig": ...
def load_config(path: str) -> "DeviceConfig": ...  # public alias, GUI/CLI 同源入口
def build_device(cfg: "DeviceConfig") -> "NormalizedDevice": ...
def contact_targets(device, bias_V: float) -> tuple[float, float, float, float, float, float]: ...
# 约定返回 (u_L, n_L, p_L, u_R, n_R, p_R)

# W6 明确行为、W7 冻结签名，W8 GUI/CLI 均以此为唯一运行通道：
def parse_config_dict(data: dict) -> "DeviceConfig": ...
def save_config(config: "DeviceConfig", path: str) -> None: ...
def run(
    config: "DeviceConfig",
    progress: "Callable[[float, str], None] | None" = None,
    should_stop: "Callable[[], bool]" = lambda: False,
) -> "RunResult": ...
def write_required_outputs(
    run: "RunResult", path: str, overwrite: bool = False
) -> tuple[str, str, str]: ...
# 固定返回 comparison.csv, comparison.png, convergence.md 的绝对/可解析路径
# 输出单位从 run.config.output.unit_system 获取（cm 默认 / si 可选）；不修改内核
# progress(bias_V, status) 仅表示已尝试偏压及 PASS/FAIL，不伪装成完成百分比
# should_stop() 为协作式取消；已接受状态必须完整保留；取消结果 REJECTED
```

接口严格要求：`State` 内的 `u,n,p` 无量纲；`bias_V` 的单位始终是 **V**（cm 与 SI 下电压单位相同，只在 `contact_targets` 转成 `v=V_app/V_T`）。默认物理 API 字段采用 `length_cm/area_cm2`、`intrinsic_density_cm3`、`mobility_n_cm2_per_Vs`、`mobility_p_cm2_per_Vs`、`reference_diffusivity_cm2_s`、`x_nodes_cm` 等明确后缀；SI 输入不允许混用这些字段冒充 cm，由 `units.system: si` 及其专用适配映射转换为 cm 配置。Week02 教学的 `pack_state(phi,n,p)` 若使用 `phi[V]`，**不得**不经归一化直接作为 `State.u`。`NormalizedDevice` 保存归一化网格/掺杂/迁移率/寿命及 **cm 基准尺度**（`L0_cm/ni_cm3/D0_cm2_s/J0_A_cm2/R0_cm3_s`）和温度等可复现元数据；可选 SI 字段只属于适配层。API 可选择求解器，但不得影响方程。输出 schema 由 `output.unit_system: cm|si` 确定，默认 cm；writer 只对已恢复的物理量进行单位/列名转换，不改变 `State` 或求解结果。**面向 GUI 的 `RunResult` 至少公开 `attempts`、`successful`（已接受偏压点）、`decision`（`ACCEPTED/REJECTED`）、`failure_reason`、`config`、`device`；单点尝试保留 `bias_V/status/result/runtime_s`。GUI 只能从已接受点绘制物理曲线，同时显示失败尝试和最终验收结论。**实现时 `State` 数组必须防止跨模块原地修改；`SolveResult(converged=False)` 的 `state` 不得被后续 sweep 当作已收敛起点。

### 10.3 A/B 的明确交接标准

| 交接项 | 提供者 → 接收者 | 交付证据 / 失败策略 |
|---|---|---|
| `Mesh1D`, `NormalizedDevice`, Dirichlet targets | B → A | pytest：cm 网格/掺杂、严格递增、结点/断点、左右 P/N 判别、cm↔无量纲；SI↔cm 适配往返 |
| Stable Bernoulli, edge flux, residual | A → B | 独立边流向/平衡恒等式、Poisson 非均匀二次函数回归 |
| `solve_gummel` / `solve_newton` | A → B | 同一接口、结构化异常、初值与容差显式参数、迭代历史非空 |
| 独立 FD Jacobian oracle & reference benchmark | B → A | 小网格 FD 方向导数与物理不变量；失败列明具体变量/节点 |
| Fixed/adaptive sweep `RunResult` | B ↔ A | 失败回滚和“上一成功解”单元测试；共用 bias 定义 |
| CSV / PNG / convergence writer | B → A、最终用户 | 字段/理论可用性检查；固定三文件路径；空理论值禁止填 `0` |
| **W7 GUI 公共 API 冻结** | **B → A** | `load_config/parse_config_dict/save_config/run/write_required_outputs`、`RunResult` 字段、`progress/should_stop`，样例 YAML 和异常语义通过契约测试 |
| **W8 PySide6 GUI** | **A → B** | 窗口/表单/高级 YAML、`QThread` Worker、真实 I–V/剖面/报告、协作式停止；不从 Qt 主线程求解，不向 GUI 引入新 PDE |
| **W8 Windows 发布与端到端验收** | **B ↔ A** | B 检查 GUI 与 CLI 同配置结果及三文件，执行 T16–T18；A 修复 UI/适配层缺陷；无实际 Windows 证据不能标 Windows EXE 已验收 |

### 10.4 Git 与代码评审纪律

- 每一模块一个明确 Owner；所有涉及 `F(X)`、SG、SRH、接触 BC、CSV schema 的修改由**两人共同 review**。
- 每个 PR 附：改动的方程编号、单元测试、含旧错误的回归测试、输入参数与最小复现；禁止“运行成功”的截图替代测试。
- 每周保存能够重建的 tag（例如 `w04-gate1`、`w06-gate2`、`w07-core-freeze`、`w08-release`），报告记录提交哈希；**W7 末冻结物理核心/公共接口，W8 仅允许 GUI、联调与阻塞性 bug 修复，不得加入新物理模型**。

***

## 11. 独立验证、理论基准和量化通过标准

### 11.1 必须区分的三种对照

1. **严格理论值（exact）**：接触电中性/质量作用、解析接触浓度、接触势、热平衡 `np=n_i²`、平衡零电流、Poisson 的制造解如 `u=ξ²` 的二阶导数。这些可设数值测试的硬阈值。
2. **受假设限制的理论近似（approximation）**：耗尽近似 `V_bi≈V_T ln(N_A N_D/n_i²)`、理想 Shockley 指数式（仅在给定低注入、准中性区、边界/复合模型满足其条件时启用）。图中必须加 `approx` 和适用条件，不能在任意 SRH 偏置下作为精确基准。
3. **独立数值参考（reference_numeric）**：同物理模型的更细网格、高精度容差或独立代码的解，用于一般偏压剖面和 `I(V)`；放在 CSV 的 `reference_value` 而非 `theory_value`。如未获得高可信参考，留空并在报告说明。

**可选的解析 I–V 近似（不是对任意 SRH 模型的精确解）**：在非简并、突变结、低注入、准中性区远大于少子扩散长度、耗尽区复合可忽略且端接触条件适用时，可以计算

$$
I_{\mathrm{Shockley}}(V)=I_s\big(e^{V/V_T}-1\big),\qquad
I_s=qA n_i^2\left(\frac{D_p}{L_p N_D^{(R)}}+\frac{D_n}{L_n N_A^{(L)}}\right),
$$

$$
L_n=\sqrt{D_n\tau_n},\qquad L_p=\sqrt{D_p\tau_p}.
$$

这些条件必须由 benchmark 配置和近似误差说明明确满足，才可将其写入 `theory_value` 且标 `theory_kind=approximation`；否则即使公式可计算也**不得作为当前全耦合 SRH 结果的“真实理论值”**。用于强 SRH 耗尽区复合、有限长度准中性区或高注入时，它可能严重偏离数值解。

**关键限制**：常规有偏压的完整 1D PN Poisson–DD–SRH 系统没有适用于任意参数的简单闭式解。不能为了让“理论列”看起来齐全而伪造 `φ(x),n(x),p(x),J(V)` 的解析值。

### 11.2 测试矩阵（Owner B 执行、A 必须修复物理失败）

| 测试 ID | 输入与独立预期 | 建议通过阈值 / 8 周执行窗口 |
|---|---|---|
| `T01_BERNOULLI` | `z=0, ±1e-14, ±1e-3` 两侧、`±1, ±50` 两侧、`±700, ±1000`；高精度 Decimal/解析参照；`B(z)−B(−z)=−z` | W2 独立自测：可表示值相对误差目标 `≤3e-14`；**W3–W4 Gate 1** 独立核验 `≤1e-11`，无 NaN/Inf |
| `T02_SRH_SWAP` | `n=10,p=2,τ̂_n=1,τ̂_p=4`，理论 `19/47`，拒绝错误的 `19/23` | 绝对误差 `≤1e-12`；W3 |
| `T03_CONTACT` | 多组 `C_L<0,C_R>0`；`np=1`、`p−n+C=0`、`u_L=a(C_L)+v` | 相对中性误差及 $\lvert np-1\rvert$ `≤1e-10`；W2–W4 |
| `T04_POTENTIAL_POLARITY` | 对比 `V=0` 与 `V>0`，势垒 `u_R−u_L` 减少 `v` | 接触电势偏差 `≤1e-10`（无量纲）；错误方向拒绝；W2–W4 |
| `T05_POISSON_NONUNIFORM` | `ξ=[0,0.25,0.75,1.5]`、`u=ξ²`；二阶导数恰为 `2`，旧分母错误给 `1.2` | 误差 `≤1e-11`；W2 专项、W3 集成、**W7 变网格终审** |
| `T06_SG_EQUILIBRIUM` | `n_R=n_L exp(Δu)`、`p_R=p_L exp(−Δu)`，双 SG 通量为零 | 归一化误差 `≤1e-11`；W3–W4 |
| `T07_SG_DIFFUSION` | `Δu=0` 时纯扩散 `+m_n Δn/h`、`−m_p Δp/h` | 相对/绝对误差 `≤1e-11`；W3–W4 |
| `T08_EQUILIBRIUM` | 完整 0 V PN：`np≈1`、边 `j_n,j_p≈0`、边界中性 | 质量作用偏差 `≤1e-6`；归一化电流 `≤1e-6`；**W4 Gate 1** |
| `T09_CURRENT` | 每个收敛偏压的总边电流 `j_n+j_p` 沿器件恒定 | $\max_i\lvert J_i-\mathrm{median}(J)\rvert/\max(J_0,\lvert\mathrm{median}(J)\rvert)\le 10^{-3}$；W4 初验，W6–W7 复验 |
| `T10_JACOBIAN` | B 独立中心/安全单边 FD，对比 A 的 `J·d` | $\lVert J_A d-J_{\mathrm{FD}}d\rVert_\infty/\max(1,\lVert J_{\mathrm{FD}}d\rVert_\infty)\le 10^{-4}$；W5 初测、**W6 Gate 2** |
| `T11_TWO_SOLVERS` | 相同器件/网格/偏压/容差下 Gummel vs Newton | 势峰值差 `≤1e-3 V` 且非零终端电流差（参考尺度保护）`≤1%`；**W6 Gate 2**，W7 压力复测 |
| `T12_MESH` | uniform 与结区加密网格逐级比较 `I`、势垒和峰值 `E` | 全部差值入报告；最细两级非零 `I` 相对差目标 `≤1%`；接近零改报绝对差；W6 建参考，**W7 签核** |
| `T13_SWEEP_ROLLBACK` | 人工失败注入；fixed/adaptive 均测试 | 失败态不传播；`ΔV` 正确收缩，失败原因真实记录；W5 初测、**W6 Gate 2** |
| `T14_OUTPUT` | 正常/失败、零理论/无理论、非均匀网格、双 solver；默认 cm 与可选 SI writer 的列名/quantity/unit、位置、理论/参考/误差一致性 | 固定三文件均可解析；cm↔SI 双向换算和 SI 仅输出层变换通过，缺失理论不填 `0`；W6 初版，**W7 完整验收** |
| `T15_REPRODUCIBLE` | 同配置、同版本/依赖/种子重复执行 | 关键值在声明容差内一致；运行元数据、SHA、完整配置齐全；W7–W8 |
| **`T16_GUI_CONFIG`** | W8 GUI 加载/编辑/保存 YAML；默认 cm 参数面板与高级 YAML；`units.system`、`output.unit_system`、非法温度/方向/单位输入及 SI 适配往返 | A 的 GUI 只能调用 B 的配置校验；非法输入禁止启动；保存再读取的规范化配置等价；**W8 Gate 3** |
| **`T17_GUI_RUN_CANCEL`** | GUI/CLI 相同 YAML 计算与画图；后台线程、进度、取消与失败点 | 窗口运行期间无阻塞；进度仅真实尝试；停止返回结构化 `CANCELLED/REJECTED` 并保留已接受点及三文件；图与 `RunResult` 一致；**W8 Gate 3** |
| **`T18_WINDOWS_EXE`** | Windows 10/11 x64 onedir 构建、Qt 启动烟雾测试、干净机器加载 YAML/运行输出 | EXE 可启动，依赖随目录完整分发，0→0.30 V 温和基准完成、三文件可验收；如未在 Windows 实测，**不得宣称 EXE 发布验收通过**；**W8 Gate 3** |

阈值中 `≤1%` 等是**验收目标**而非项目现在已经通过的结果。GUI/EXE 测试是本版新增的 W8 Gate 3 必测项；**数值基准通过与 Windows 打包通过分别给出判定**，不能相互代替。凡选择其他阈值，必须在项目启动前或通过双方审核的基准修订中登记；不能为了掩盖失败而事后调整。

### 11.3 特别注意“零值相对误差”

对具有有效理论值 `Q_th` 的结果：

$$
\mathrm{abs\_error}=|Q_{\rm num}-Q_{\rm th}|,
\quad
\mathrm{rel\_error}=\frac{|Q_{\rm num}-Q_{\rm th}|}{|Q_{\rm th}|}
\quad\text{仅在 }|Q_{\rm th}|>Q_{\rm floor}\text{ 时定义}.
$$

`Q_floor` 由测试规则按该物理量**实际输出单位**给出，并记录于报告；SI 输出时 `Q_floor` 和绝对误差阈值同物理量同步换算，不能原样沿用 cm 的数值。若 `Q_th=0` 或低于下限：`theory_rel_error` 在 CSV **留空**，图中标 `N/A (zero reference)`，改报告**绝对误差**（如热平衡电流密度）。不可写 `0%`，也不可除以零。独立数值参考同理，用**独立的** `reference_abs_error/reference_rel_error` 字段。

***
## 12. 固定验收 Output 契约（不能增加“可选第四份报告”）

一次 `minitcad solve CONFIG --out OUTPUT_DIR` 或 GUI 通过同一 `run()`/writer 发起的正式计算，必须在指定目录写出下列 **3 个固定文件**，并遵守本节规定的字段、内容及失败处理契约。

```text
OUTPUT_DIR/
├── comparison.csv       # 一张完整长表：全部数值 + 有根据的理论值/参考值 + 误差 + 迭代序列
├── comparison.png       # 一张图：数值/理论/误差、IV、剖面、守恒与收敛
└── convergence.md       # 一份简报：每个偏压迭代次数、残差/更新容差、验证状态、复现配置
```

**GUI 建议为每次运行创建独立 `run_时间戳` 子目录；桌面 EXE、`.vscode/`、周报、Qt 日志、构建产物均不得混入本次求解的三文件目录。**`--out` 目录不得靠覆盖已有文件制造旧数据与新图混杂：默认若三文件已存在则要求 `--overwrite`，覆盖时三个文件**原子替换**或用临时目录统一发布；失败也生成这三个文件，但在文件内标明 `FAILED/INCOMPLETE` 和原因。自动化验证根据这三个**确切文件名**进行检查，不能换成 `result_final2.csv`。

### 12.1 文件一：`comparison.csv` —— 统一完整数据 + 理论 + 误差

**格式**：UTF-8 with BOM (`utf-8-sig`，便于中文 Excel 打开)，RFC 4180 CSV，使用英文 ASCII 列名、句点小数，缺失数据为空字段而非伪数字 `0`、字符串 `NaN`。采用 **long-form tidy schema**：每行描述一个物理量在特定偏压/空间位置/迭代步的一个值；禁止把 201 个结点塞入一个单元格。

**默认 cm 模式强制列（不可随意删除或改名，后续仅可附加列）**；显式 SI 输出使用下文定义的**等价 SI schema**，不允许在一张 CSV 中混用两个单位模式：

```csv
run_id,config_sha256,solver,bias_V,record_type,quantity,location_kind,index,position_cm,iteration,unit,numeric_value,theory_value,theory_kind,theory_abs_error,theory_rel_error,reference_value,reference_kind,reference_abs_error,reference_rel_error,status,notes
```

显式 `output.unit_system: si` 时，**只在 writer** 将此处 `position_cm` 改为 `position_m`，并按下述 `quantity` 命名映射/`unit` 做物理量转换；其余表头顺序不变。两种模式都必须在 `convergence.md` 的 Physics/units 中记录 `input_unit_system`、`core_unit_system=dimensionless`、`output_unit_system`，以便自动判定 schema。

各列定义：

| 列 | 严格约定 |
|---|---|
| `run_id,config_sha256` | 唯一运行标识及完整规范化配置 SHA-256，能与报告对应 |
| `solver,bias_V` | `gummel/newton`、V（两个单位制均不换算）；每个成功或失败偏压均须出现记录 |
| `record_type` | `node / edge / sweep / benchmark / iteration / diagnostic` |
| `quantity,unit` | 量名须带相应单位后缀，默认以 cm 工程单位表示；`unit` 如 `cm^-3`、`V/cm`、`A/cm^2`、`C/cm^3`、`cm^-3/s`，无量纲填 `1`；SI 文件改写后缀及 `unit`，不得混合 |
| `location_kind,index,position_cm` | 默认以 `cm` 记录物理位置：`node` 索引 `0..M−1`；`edge` 索引 `0..M−2`、位置取区间中心；无空间位置则留空。SI 文件仅此列改名 `position_m` 且数值乘 `1e-2` |
| `iteration` | 求解迭代编号；非迭代记录留空；`iteration` 行包含 `residual_norm/update_norm/alpha` 等量 |
| `numeric_value` | 仿真数值，成功记录不能缺失或非有限；布尔量 `converged` 用 `1/0` 表示；失败未求出的量留空并标状态 |
| `theory_value,theory_kind` | 理论值只允许有明确依据时填；类别 `exact / approximation / none` |
| `theory_abs_error,theory_rel_error` | 仅相对于理论值；理论不存在则均为空；理论为零时相对误差为空 |
| `reference_value,reference_kind` | 独立数值基准值和 `refined_mesh / external_solver / none` |
| `reference_abs_error,reference_rel_error` | 仅相对于数值基准；未提供时为空 |
| `status,notes` | `PASS / FAIL / N_A`；误差阈值、近零理论、近似适用域等原因写在 `notes` |

**`numeric_value` 至少覆盖的完整数据**：

| record_type | 必须写入的 quantity | 位置信息 |
|---|---|---|
| `node` | `phi_V`, `n_cm3`, `p_cm3`, `ND_cm3`, `NA_cm3`, `net_doping_cm3`, `charge_C_cm3`, `SRH_cm3_s`, `F_poisson_nd`, `F_electron_nd`, `F_hole_nd` | 所有结点、所有成功偏压 |
| `edge` | `E_V_cm`, `Jn_A_cm2`, `Jp_A_cm2`, `Jtotal_A_cm2` | 所有边、所有成功偏压 |
| `sweep` | `I_A`, `Jterminal_A_cm2`, `applied_voltage_V`, `iteration_count`, `runtime_s`, `converged` | 每个尝试偏压，成功/失败均记录；给定面积时 `I_A=A_cm2·Jterminal_A_cm2`；面积未提供时 `I_A` 留空并在 `notes` 标 `AREA_NOT_PROVIDED`，不得虚构电流 |
| `benchmark` | `contact_n_cm3`, `contact_p_cm3`, `contact_potential_V`, `equilibrium_np_ratio_nd`, `equilibrium_zero_current_A_cm2`, `built_in_potential_V` 等适用项 | 有理论公式时写明对应 cm 工程单位；`equilibrium_np_ratio_nd=np/n_i²` 的理论值为 `1` |
| `iteration` | `residual_norm`, `update_norm`, `alpha`, `trial_status`（状态采用 `notes`/`status`） | 每轮迭代，含失败的偏压尝试 |
| `diagnostic` | `current_nonconservation_nd`, `min_n_nd`, `min_p_nd`, `boundary_max_error_nd` 等 | 每个收敛点和失败终点；此行所有带 `_nd` 的量 `unit=1` |

**两套文件 schema 的唯一单位映射**：`n_cm3/p_cm3/ND_cm3/NA_cm3/net_doping_cm3/contact_n_cm3/contact_p_cm3 → *_m3`（数值 `×1e6`）；`charge_C_cm3 → charge_C_m3`、`SRH_cm3_s → SRH_m3_s`（`×1e6`）；`E_V_cm → E_V_m`（`×1e2`）；`Jn_A_cm2/Jp_A_cm2/Jtotal_A_cm2/Jterminal_A_cm2/equilibrium_zero_current_A_cm2 → *_A_m2`（`×1e4`）；`position_cm → position_m`（`×1e-2`）。`phi_V/I_A/applied_voltage_V`、所有无量纲量、迭代与时间值不变。`numeric_value/theory_value/reference_value/theory_abs_error/reference_abs_error` 必须按该行物理单位**同步**换算；相对误差不变，`unit` 则同步切换。不得只改单位文字、不改数值，也不得只改数值、不改量名。

有有效相对能带模型时可附加 `Ec_relative_eV,Ev_relative_eV` 节点行；无模型则不得臆造绝对带边。CSV 必须包含所有实际运行的偏压（包括 adaptive 插入点）；若求解器只执行某些偏压，应如实列出而非补插值当作求解值。

**理论/误差规则**：

- 精确理论值：左/右接触 `n,p,φ`，0 V 的 `np=n_i²`、零电流、内建电势差；条件满足时才记录有条件理论近似。
- 对 0 V 的 `J_theory=0`，计算 `theory_abs_error=|J_numeric|`，`theory_rel_error` 留空，在 `notes` 写 `zero theory: relative error undefined`。
- 有偏压时各节点的 `φ,n,p` 一般**无**通用闭式解，`theory_kind=none`、`theory_value` 留空。可将细网格插值到同一物理位置的结果填 `reference_value` 并明确其来源；两者的误差**分开**计算。
- 近零阈值 `Q_floor` 要按物理单位和测试约定存入 `convergence.md`，零阈值检查前不得对绝对值小的理论值不加说明地相除。
- CSV `status=PASS` 仅代表所在记录对应的测试/偏压通过；全局是否通过由报告汇总，各文件 `run_id` 一致。

### 12.2 文件二：`comparison.png` —— 一幅真正能验收的综合图

**一个 PNG 文件，多面板（推荐 3000×2000 px 或更高，可读字号）**，默认轴单位 `x[cm]`、`E[V/cm]`、`J[A/cm²]`、`n,p[cm⁻³]`，SI 输出时仅在绘图适配层改为 `m`、`V/m`、`A/m²`、`m⁻³`，数值依单位同步换算；至少显示：

1. `I–V` 全扫压数值曲线；若理论/数值参考合法，叠加 **Theory (exact/approx)** 或 **Reference (fine mesh)**，两者用不同线型，并标条件。
2. 有代表性的已求解偏压（明确列出所选 `V`）的 `φ(x)`、`n(x),p(x)`，n/p 使用对数轴；若有对照，叠加参考/理论曲线，无则标注 `Theory: N/A`。
3. `E(x)` 和 `Jn(x), Jp(x), Jtotal(x)`；读者能够看到总电流沿空间接近常数。
4. **绝对误差/相对误差**面板：理论存在且不为零时显示 `Theory relative error`；零理论显示 `Absolute error`；有独立数值参考则单独显示 `Reference relative error`。三个概念不得合并成一个模糊“误差”。
5. `iterations vs bias` 和/或 `residual vs iteration`；最好按 Newton/Gummel 颜色区分，同时记录未收敛点。
6. 一张小型 `Benchmark summary` 表，列出主要可用量的 `Numerical / Theory / Rel. error or Abs. error / Type / PASS-FAIL`；其中理论不可用时显示 `N/A` 而不是 0。
7. 头部和图例显著标注 **`P(left) → N(right); V_app = V_left − V_right; forward > 0 V`**；标题/轴/图例**使用英文**，以便最终技术报告引用。

**显示完整性的边界**：CSV 按行保存**所有**数据；一张 PNG 无法逐格写完全部结点×偏压的数字，所以图像要求完整覆盖“数值、理论、误差”这**三类信息**、全部 `I–V` 点和代表性空间剖面，不要求将每个 CSV 单元格以不可读小字印到图片上。图中选择的偏压、参考数据来源、未展示的曲线均在图注中说明。缺失理论时不许绘制造假的理论曲线。

### 12.3 文件三：`convergence.md` —— 简短但信息足够的收敛报告

**固定结构**：

```text
1. Run summary        run_id, 代码版本/git sha, 配置 SHA256, 日期, 成败总览
2. Orientation        P(left) → N(right); V_app=V_left−V_right; +V=forward
3. Physics/units      input_unit_system(cm/si)、core_unit_system=dimensionless、output_unit_system(cm/si)；
                      生效的 cm 归一化尺度及 T/材料来源；物理结果按 output 单位显示：
                      cm→ε[F/cm],ni[cm^-3],μ[cm²/(V·s)],L0[cm],D0[cm²/s],J0[A/cm²],R0[cm^-3/s]；
                      si→ε[F/m],ni[m^-3],μ[m²/(V·s)],L0[m],D0[m²/s],J0[A/m²],R0[m^-3/s]；
                      τ[s]、VT[V] 不变，记录 SI 输入/输出适配和 cm 内核尺度的换算审计
4. Solver settings    method, Jacobian, max_iter, alpha_min, n_floor
5. Tolerances         tol_residual, tol_update, current conservation 等
6. Bias table         每个尝试 V, solver, step, success?, iterations,
                     final Rnorm, final Dnorm, runtime, fail reason
7. Validation         Bernoulli/SRH/BC/equilibrium/current/Jacobian/mesh
                     的指标数值、阈值及 PASS/FAIL/N/A
8. Error policy       exact/approx/reference 适用性，各量 Q_floor
9. Reproducibility    Python/NumPy/SciPy 版本、随机种子（若使用）,
                     完整 YAML 配置快照（以代码块嵌入）及运行命令
10. Final decision    ACCEPTED 或 REJECTED + 未通过项（如有）
```

当 `output.unit_system: si` 时，收敛报告中的物理量、基准尺度、绝对误差及 `Q_floor` 同步转换为 SI，并保留用于复现的原始 cm 归一化尺度审计记录。这仍然只是文件输出层转换，不改变数值核心或 Gate 的无量纲判据。

报告**不是第二份巨大的 CSV**：保存每偏压迭代总数、最终双判据及关键物理误差；逐次迭代完整序列由 `comparison.csv` 承担。若 `Newton`/`Gummel` 两者均运行，分别给出迭代和容差结果，不能用一个“平均收敛次数”遮盖失败。

### 12.4 三文件的自动化验收脚本必须检查

- 三个文件名完全正确、都能解码/读取、`run_id` 与输入配置哈希一致；PNG 可打开且非空，报告有所有固定章节。`units.system`、`output.unit_system` 与 CSV 表头、`quantity`/`unit`、PNG 轴、报告尺度一致；默认 cm、显式 SI 输出分别按各自确定的 schema 验收。
- 在 CSV 中，成功点的节点数量恰为 `M`，边数量 `M−1`，对应所要求的每一 quantity；无法满足就报 `DATA_INCOMPLETE`。
- 有理论值时误差能独立重算：`abs_error=abs(numeric−theory)`；合法非零理论时 `rel_error=abs_error/abs(theory)`；`theory_kind=none` 时理论值和相关误差必须留空。
- 每个偏压的 `iteration` 行数与报告迭代次数一致（约定是否计入第 0 次初始化）；默认终端 `I_A=A_cm2·Jterminal_A_cm2`，SI 文件检查 `I_A=A_m2·Jterminal_A_m2`；同一结果在双输出单位下的电流应一致。方向说明在报告和 PNG 出现。
- 失败的求解仍保留 `FAIL` 记录和最后已接受状态，但**最终决定为 `REJECTED`**，不能仅因为三文件存在就算通过。
- **单位专项回归**：cm↔SI 双向参数往返、cm↔无量纲物理量往返、同器件 cm/SI 输入适配得到同一无量纲结果；要求 SI 输出时，比较两个 writer 的所有物理数据、理论/参考/绝对误差、零值相对误差空白、空间坐标和图轴，确认转换发生在输出层，且不增加第四份文件。
- **W8 GUI/CLI 对照**：同一规范化 YAML 的 GUI 与 CLI 必须调用 B 的公共 API，检验 config SHA、偏压点、物理量、失败原因及三文件 schema 对应；GUI 切换结果页不能以显示成功替代 `convergence.md` 的 Final decision。

***

## 13. 可以直接落地的 YAML 示例（左 P / 右 N；默认 cm 工程单位）

> 文件名 `examples/pn_left_p_right_n.yaml`。以下数值保持原 SI 示例**同一物理器件**：长度 `1.0e-3 cm`、面积 `1.0e-4 cm²`、结位置 `2.5e-4 cm`、本征浓度 `1.0e+10 cm^-3`。这些是**配置启动示例**，不是已通过仿真的证据；科学计数采用显式指数符号，避免 YAML 解析器差异。

```yaml
# PN orientation: LEFT=P (x=0), RIGHT=N (x=L)
# V_app = V_left - V_right; V_app > 0 is forward bias; right terminal = 0 V.
# Input: cm engineering units; dimensionless physics core; SI only via adapters.
units:
  system: cm                  # 输入物理量必须声明，项目默认工程单位为 cm

device:
  type: pn
  orientation: p_left_n_right
  length_cm: 1.0e-3
  area_cm2: 1.0e-4
  temperature_K: 300.0       # 示例温度；不可把 300 K 当作隐含固定材料模型

material:
  name: silicon
  # 以下参数仅针对本示例 T=300 K；改 T 必须重审 n_i/迁移率
  epsilon_relative: 11.7
  intrinsic_density_cm3: 1.0e+10
  mobility_n_cm2_per_Vs: 1350.0
  mobility_p_cm2_per_Vs: 480.0
  lifetime_n_s: 1.0e-6
  lifetime_p_s: 2.0e-6
  bandgap_eV: 1.12

regions:
  - start_cm: 0.0
    end_cm: 2.5e-4
    material: silicon
    acceptor_cm3: 1.0e+16
    donor_cm3: 0.0
  - start_cm: 2.5e-4
    end_cm: 1.0e-3
    material: silicon
    acceptor_cm3: 0.0
    donor_cm3: 2.0e+15

mesh:
  type: uniform
  nodes: 201

normalization:
  reference_diffusivity_cm2_s: 10.0  # 明确指定 D0[cm²/s]；不得悄悄更换基准

contacts:
  type: ideal_ohmic
  left_terminal: p
  right_terminal: n
  right_voltage_V: 0.0
  bias_definition: left_minus_right

solver:
  method: newton            # or gummel
  jacobian: sparse_fd
  tol_residual: 1.0e-8
  tol_update: 1.0e-6
  max_iter: 60
  alpha_min: 1.0e-6
  positive_scale_floor: 1.0e-15

sweep:
  start_V: 0.0
  stop_V: 0.6
  initial_step_V: 0.05
  min_step_V: 0.005
  max_step_V: 0.1
  grow_factor: 1.25
  shrink_factor: 0.5

output:
  directory: runs/pn_example
  unit_system: cm            # 默认 cm；如须交付 SI 文件，改为 si，只有 writer 转换
  # Names are fixed: comparison.csv, comparison.png, convergence.md
```

**单位契约**：`units.system: cm` 时只能使用 `_cm/_cm2/_cm3/_cm2_per_Vs/_cm2_s` 等工程字段；`units.system: si` 时只能使用各自 `_m/_m2/_m3/_m2_per_Vs/_m2_s` 字段，并在**输入适配层**变换为等价 cm 参数后归一化。二者不得在同一 YAML 中混用。`output.unit_system` 仅控制最终文件是 cm 还是 SI，**不控制求解器、尺度公式或数值状态**。若仅有 SI 文件交付要求，保持上述 `units.system: cm` 和全部 cm 数值不动，只将 `output.unit_system: si`；writer 以第 12.1 节的 SI schema 输出。

**必须拒绝的配置**：负温度或寿命、非正 `n_i`、非正面积、负 `N_D/N_A`、region 重叠/缺口、网格未覆盖器件、PN 两端方向反转、`stop>start` 而初始步长非正、`min_step>max_step`、不可解析数值/单位、`units.system` 不合法、物理字段单位后缀与声明不一致、`output.unit_system` 不属于 `cm/si`、非法 solver 名称。新版示例要求显式 `units.system`；迁移旧 YAML 须通过明确的兼容入口指定其单位来源，不能按数值大小猜测。`generic_1d` 不自动套用 PN 的左右方向特别校验，但不得使用错误的 `p_left_n_right` 标签。

### 13.1 CLI / Python API / W8 GUI 的用户体验

```bash
minitcad solve examples/pn_left_p_right_n.yaml --out runs/pn_example
minitcad validate runs/pn_example/comparison.csv --report runs/pn_example/convergence.md
```

第二个命令只读取现有结果和报告，**不额外生成第四个验收文件**；如需重算验收图，应显式重新 `solve`。Python API：

```python
from minitcad.api import load_config, run

cfg = load_config("examples/pn_left_p_right_n.yaml")
result = run(cfg)
# result is a RunResult; its writer emits the three fixed deliverables.
```

`solve` 默认回显偏置方向、**输入 cm/SI 声明、内核无量纲、输出 cm/SI 选择**、材料范围、节点数、方法、接受或失败原因及最终三个文件的路径；绝不把反向误称正向，也不在失败时返回看似成功的 IV 曲线。

**W8 GUI 的等价入口**（启用可选 `gui` 依赖后）：

```bash
python -m minitcad.gui.app
# 可选：用于自动化启动检查，不是数值仿真结果
python -m minitcad.gui.app --smoke-test qt_smoke.txt
```

GUI 必须通过第 10.2 节的 B 公共接口加载 YAML、运行和写文件。Qt 启动烟雾测试文件是**测试目录的内部工件**，不是求解输出目录的第四个文件。Windows 打包和实际验证在第 14.9 节规定。

***

## 14. 8 周、两位负责人、完整实施任务与 GUI 交付 Gate

### 14.1 角色、最终责任与实施原则

**A：物理/离散/数值求解 Owner（W1–W7）＋ PySide6 GUI Owner（W8）。** W1–W7 必须完成 cm 尺度/无量纲公式、SRH、稳定 Bernoulli、Poisson FV、SG、两条 continuity、唯一 `F(X)`、Equilibrium、Gummel、Fully-coupled Newton、稀疏 Jacobian、阻尼/正性、强偏压、非均匀网格、网格精度与性能、派生量、双求解器比较和物理终审，并提供结构化 `SolveResult` 与可测数值入口。W8 是独立的 GUI 开发和接口联调周，A 按第 14.8 节逐日交付真实可操作的桌面程序。

**B：器件/工作流/独立验证/交付 Owner（W1–W8）。** 必须完成 Material/Region/Mesh/Device/Contact/State、`DeviceConfig`/YAML、cm 默认工程单位校验与 SI 输入/输出适配、独立物理/Jacobian oracle、fixed/adaptive continuation、细网格参考实验、三文件 writer、Python API、CLI、自动验收、复现审计、安装/用户说明、Windows 构建与发布；其中 **W7 末必须冻结 GUI 所需的配置校验、`run()`/`RunResult`、进度/取消契约、writer 与结果解释接口**，W8 只做联调、独立端到端验收和发布收口。

**共同约束**：A 负责物理真相、计算稳定性和 GUI 代码；B 负责独立反证、稳定接口、最终科学判定和发布。涉及端电压/电流方向、SRH、Poisson/SG 符号、`F(X)`、配置 schema、GUI 公开接口的变更必须双人 review。B 的 oracle 不得复制 A 的实现。**T01–T15 全部保留且阈值不降低，失败证据不得删除；W8 另执行 GUI/EXE 的 T16–T18。**

**四个两周实施阶段**：

- **基础准备 W1–W2**：统一物理/单位与 A/B 数据契约；W2 Bernoulli、Poisson 教学与正式无量纲版本必须完成。
- **Gate 1 集成 W3–W4**：SG/SRH/唯一残差 → Equilibrium/Gummel → 首个非零偏压，W4 末裁决。
- **Gate 2 集成 W5–W6**：Newton/Jacobian/稳定性 → fixed/adaptive sweep 与双求解器独立比较，W6 末裁决。
- **Gate 3 集成 W7–W8**：W7 核心终审、网格与三文件完整性、接口冻结；W8 A 开发 GUI、B 联调与 Windows 发布验证。

### 14.2 每周任务、不可推迟项与可检查证据

| 周 | A：本周必须完成 | B：本周必须完成 | 本周验收证据 / 阻塞项 |
|---|---|---|---|
| **W1** | 锁定方程符号、cm/无量纲尺度（保留 SI 适配与教学核对）、左 P 右 N 接触、SRH 寿命；列出数值压力场景与接口变量 | 建立项目骨架、`DeviceConfig`/CSV schema/报告模板、示例 cm 配置、测试与代码评审框架；写方向提示 | `math_contract.md`、初版 API/字段表、物理假设及首批回归清单；**A/B 双签** |
| **W2** | 完成 `src/minitcad/discretization/bernoulli.py`、正式无量纲 `src/minitcad/discretization/poisson.py`、`tests/test_bernoulli.py`、`tests/test_poisson.py`；用独立 SI MMS 验证器核查 Poisson 符号/尺度，并冻结 SG/continuity 下游签名 | 完成 `Mesh1D`、Region/Material/Contact/State/Device、cm 单位校验与 SI↔cm 双向适配、接触解析基准、CV 掺杂平均接口；会审并冻结 Mesh/State/残差单位契约 | **Week 02 专项门槛**：T01、T03/T04、T05、cm↔nd、cm↔SI、SI 教学↔nd、非均匀二次 MMS、均匀正弦加密、异常注入；失败不得进入 W3 集成 |
| **W3** | 完成 SG 电子/空穴、SRH、continuity、Dirichlet/唯一残差 `F(X)`；提前实现 equilibrium scalar 求解初版 | 完成独立 T01–T07 oracle、反向 PN/非法配置测试；提供 `RunResult`/CLI 最小骨架与 writer 的 long-form schema 样表 | 统一 residual、SG 平衡/纯扩散、SRH 异寿命及非均匀掺杂测试，A/B 接口回归全绿 |
| **W4** | **Equilibrium 正式求解、Gummel、0 V＋至少一个非零偏压、端电流/守恒检查**；解决首轮数值失败 | 完成 T08/T09、标准基准配置与数据快照、固定步进扫描；串通 B 的 `run()`、失败日志和 CLI 方向摘要 | **Gate 1（W4）**：真实 0 V/非零 V 解、质量作用/零电流/守恒及可重跑记录全部合格 |
| **W5** | 完成 Fully-coupled Newton、交错 `3M` residual/Jacobian 索引、稀疏矩阵装配；加入阻尼、正性与非有限状态检查、残差/更新双判据、数值异常回退；输出结构化 `SolveResult` 与完整迭代历史 | 建立独立 FD Jacobian oracle（小网格中心/安全单边差分、自动扰动与方向检查）；实现 fixed/adaptive continuation、失败注入/回滚、步长收缩/恢复；冻结 `run(config,progress,should_stop)` 与 GUI 协作式取消语义 | T10 方向导数初测、Newton 已知平衡解、T13 失败态不传播；Jacobian 形状/非零模式和接口草案经 A/B 会审 |
| **W6** | 在相同器件/网格/偏压/容差下完成 Newton/Gummel 全 bias 交叉比较；修复 Jacobian 精度、残差/更新尺度、线搜索与强偏压收敛；做稀疏性能优化和 bug fix；形成网格/压力测试数据 | 正式通过 T10/T11/T13；形成 Gummel vs Newton 一键对比、T08–T11/T13 一键回归/CI 与故障注入工具；完成细网格参考首版（T12）；**交付可调用 CLI/API 与三文件 writer 初版**；固定 progress/取消事件契约 | **Gate 2（W6）**：双 solver、稀疏 Jacobian、adaptive rollback、强偏压失败诊断、实际收敛/失败历史全部可复现 |
| **W7** | 完成非均匀 CV+SG 全链条验证、结区局部误差修正、uniform/refined 网格精度与迭代/成本对比；复核强偏压、残差/守恒/SG/BC、派生量和速度瓶颈；完成 Newton/Gummel/Jacobian/SG/网格收敛行为终审并只修科学正确性 bug | 完成 segment mesh/结区 refine、独立细网格参考与网格收敛表；正式通过 T12/T14/T15；冻结 `comparison.csv`/`comparison.png`/`convergence.md` schema、独立验收器、CLI/API/YAML、稳定输出目录、CI 全量回归、复现审计、安装/运行/用户说明、基准打包与配置 SHA；**冻结 GUI API、`RunResult`、writer** | **`w07-core-freeze`**：A 的数值核心与 B 的工作流/三文件/文档全部可独立运行并签收；只允许登记且不影响 GUI/科学正确性的非阻塞缺陷进入 W8 |
| **W8** | 完成 PySide6 GUI（参数表单＋高级 YAML＋结果分页），实现 QThread 仿真、真实进度/停止、I–V/物理剖面/报告展示和错误提示；与 B 的冻结接口联调，只修 GUI/适配与发布阻塞缺陷 | 执行 GUI/CLI 同配置一致性、T16/T17、完整 T01–T15 回归和复现审计；**在 Windows 执行 T18、onedir EXE 构建与烟雾测试**；完成演示、发布说明、release tag 与最终判定 | **Gate 3（W8）**：可用 GUI＋统一 API＋固定三文件＋最终 `ACCEPTED/REJECTED`＋可核对的 Windows/源码发布证据 |

#### 任务覆盖矩阵（按工作包核对，不允许用笼统“剩余任务”代替）

| 工作包 | A：实现/修复责任 | B：独立验证/工作流责任 | 最晚完成与证据 |
|---|---|---|---|
| 数学、单位与接触 | cm/无量纲尺度、SI 适配、`C`/电流/SRH 符号、接触/偏压公式 | `DeviceConfig`、方向/极性提示、cm 默认校验与 SI 输入适配 | W1；`math_contract.md`、符号表、双签 |
| Bernoulli + Poisson | 稳定 `B(z)`、无量纲 FV Poisson、SI MMS 对照 | Mesh/State/Device、CV 掺杂平均、单位往返与反例 | W2；四项必交文件＋Week 02 专项门槛 |
| SG + SRH + continuity + 唯一残差 | 电子/空穴 SG、异寿命 SRH、两条 continuity、Dirichlet、`F(X)` | T01–T07 独立 oracle、非法 PN/配置反例 | W3；集成测试全绿 |
| Equilibrium + Gummel + 电流 | 平衡初始化、Gummel、0 V/非零偏压、端电流 | T08/T09、固定 benchmark、`RunResult`/CLI 最小链路 | W4；Gate 1 |
| Newton + Jacobian + 阻尼 | Coupled Newton、3M 索引、稀疏装配、阻尼/正性/回退、迭代历史 | FD Jacobian、自动扰动/方向检查 | W5–W6；T10 |
| DC continuation + 失败回滚 | 强偏压收敛、更新尺度、结构化失败诊断 | fixed/adaptive sweep、步长增减、失败回滚、偏压日志、`progress/should_stop` | W5–W6；T13、Gate 2 |
| 双求解器与性能 | Gummel/Newton 全 bias 对照、稀疏优化、残差缩放、线搜索、bug fix | T11 一键对照、失败注入、可复现诊断 | W6；Gate 2 |
| 非均匀网格与细网格参考 | 非均匀 CV+SG 全链条、结区局部误差、迭代/成本随网格变化 | segment/refined mesh、独立细网格参考、网格收敛表 | W6–W7；T05/T12 |
| 物理/派生量终审 | 残差、守恒、SG、BC、强偏压、派生量、速度瓶颈与科学 bug 终审 | 独立检查器复验、失败工况与修复前后证据 | W7；`w07-core-freeze` |
| 三文件 + CLI/API/YAML | 保证求解器结果/诊断可被稳定消费 | writer、严格 schema、自动验收器、稳定输出目录、CLI/API/YAML、零理论/无理论/失败场景 | W6 初版，W7 冻结；T14 |
| 复现、文档与发布基线 | 提供稳定可调用内核与版本信息 | 配置 SHA、依赖/平台、复现审计、安装/运行/用户说明、基准打包 | W7；T15 |
| GUI + Windows 发布 | PySide6、QThread、真实结果页、取消、错误提示、GUI 文档 | GUI/CLI 一致性、T16/T17、Windows T18、onedir、最终审计与 release tag | W8；Gate 3 |

**A 的“抗压”要求是交付密度而非跳过验证**：W3–W7 同时推进数值开发、算法复现和强偏压/网格压力问题定位；W7 必须出具失效工况、残差/失败码、瓶颈与修复记录，不能用“之后 GUI 再处理”推迟。W8 GUI 工作按第 14.8 节五天执行，不设第九周缓冲；若核心阻塞未解决，直接按 Gate 的失败策略记录而非降低阈值。

**样例分级**：随附源码包中的温和掺杂 `examples/pn_left_p_right_n.yaml` 可作为可重复 GUI/CLI 冒烟基线（例如 0→0.30 V）；高掺杂 `examples/pn_high_doping_stress.yaml` 是**数值压力场景**，其物理适用域与历史失败必须如实呈现。Gate 只对**预先登记、处于模型适用范围内**的 benchmark 作数值达标裁决；已承诺范围内的失败不能靠删除偏压点或改阈值掩盖。

### 14.3 三个 Gate 的硬性退出条件与失败处理

**Gate 1（W4，窗口 W3–W4）**：第 14.4 节四项必交文件、SI MMS、非均匀二次/正弦加密、接口签核先通过；W3 的 T01–T07 集成通过；W4 真正完成 0 V 与至少一个正向非零偏压状态，T08/T09 的质量作用、接触、零电流和总电流守恒达标。A 负责修复物理错误；B 提供独立 oracle 和失败数据。失败则 Gate 1=`REJECTED`，暂停增加算法外围功能，不得声明“完成大测试点”。

**Gate 2（W6，窗口 W5–W6）**：Newton 与 Gummel 使用唯一 `F` 和同一边界；T10 Jacobian 独立差分、T11 双求解器、T13 失败回滚全部按第 11.2 节的明确阈值通过；存在强偏压失败时必须记录偏压、迭代、残差、线搜索/回滚原因，未到终点不得假报成功。B 的 `run()` 需支持 `progress`/`should_stop`；writer 至少可生成三文件初版，供 W7 冻结与 W8 GUI 消费。

**Gate 3（W8，窗口 W7–W8）**：W7 A 数值核心与 B 的配置、公共 API、T12/T14/T15 和三文件正式通过并冻结；W8 A 完成 GUI 且在真实桌面环境打开，B 对 GUI/CLI 同源输入、进度/取消、失败与成功、T16/T17、自动回归做独立验收。**无论从 GUI 还是 CLI 发起，成功和失败都必须遵守固定三个求解输出文件和同一科学判定。** Windows EXE 若列入交付，T18 必须在 Windows 上实机/干净环境测试；仅有打包脚本或 Linux 上的结果不能写为“Windows 已通过”。

**失败与回退**：每个 Gate 保留 `PASS/FAIL/N_A`、Owner、阻塞项和复现步骤；物理/三文件硬性约束未达标则最终 `REJECTED`。若 W8 时间冲突，优先保留输入校验、真实计算、取消、三文件和验收报告，允许缩减非必要 UI 美化；**不得牺牲数学正确性、A 的 GUI 必做项或 B 的独立验证**。

### 14.4 Week 02：五日执行与验收表

本节安排与正式 Python 包布局及 A/B 交接契约一致。学习环境为 VS Code + Python 3.10+ + NumPy/SciPy + pytest。**本周只承诺“稳定 Bernoulli + Poisson 离散及其独立证明”**；SG 电子/空穴、SRH、耦合 `F(X)`、Gummel、Newton 不在本周的“已完成”清单内。作为 MMS 工具调用 `scipy.sparse.linalg.spsolve` 解**独立的线性 Poisson**，不等于已实现 Gummel。

| 日期 | A：数学/编码与教师提问 | B：同步对接/独立检查 | 当天可核对的证据 |
|---|---|---|---|
| **Day 1：Bernoulli** | 推导级数、`B(-z)=B(z)+z`、正负大数渐近；问“为什么 `expm1` 仍不能直接处理 `z=1000`？”；写分段 mask | 冻结“输入为无量纲 `z`；异常必须明确”的 API 与压力样本 | `bernoulli.py`；`0`、微小数、`±1000`、数组的运行日志 |
| **Day 2：Poisson** | 手画节点/边/对偶控制体积；推导负号及两端 Dirichlet；问“`φ=(x/L)^2` 为什么对应负 `ρ`？” | 提供独立 SI 教学网格、物理浓度/电荷、接触电势验证样本（生产仍为 cm）；确认内部 `C/m³` 与边界 `V` 仅限 SI 教学 | SI MMS 验证器与正式无量纲 `poisson.py` 的签名、直线/二次测试 |
| **Day 3：非均匀与状态** | 用 `np.diff(x)` 和 `w=(h_left+h_right)/2`，拒绝重复点；演示 pack/unpack 顺序 | 交付正式 `Mesh1D` 的 `x_nodes_nd`、`h_edges_nd`、`w_interior_nd` 及控制体掺杂平均 | 非均匀二次制造解、形状测试、A/B 接口评审记录 |
| **Day 4：数值证据** | 正弦 MMS 在均匀网格求线性 Poisson；打印误差比；故意改变负号、单位、一个固定 `dx` | 不重用 A 的 Poisson 实现生成“真实解”；用解析正弦/独立线性代数交叉核对 | 收敛误差表、错误注入记录、非有限输入失败日志 |
| **Day 5：验收** | 从干净环境运行所有单测；梳理 Week03 的 SG→continuity→`F(X)` 接口 | 复核四项必交文件、cm↔无量纲、cm↔SI 双向测试与 SI 教学 oracle、版本及命令；联合签核 | `reports/week02_acceptance.md` 内部周报及双签、可复现提交 hash |

**四项必交文件（Owner A）**：

1. `src/minitcad/discretization/bernoulli.py`：稳定 `B(z)`，标量/数组与非有限输入策略完整。
2. `src/minitcad/discretization/poisson.py`：正式**无量纲** CV Poisson，接受 B 的网格与边界目标。
3. `tests/test_bernoulli.py`：独立精度、分支边界、反射恒等式、极值、空数组与异常。
4. `tests/test_poisson.py`：SI/无量纲对应、直线/二次/正弦 MMS、Dirichlet、网格形状与非法输入。

`tests/test_state_layout.py`、`examples/verify_poisson.py` 和 `.vscode/` 设置是强烈建议的**内部辅助材料**，不增加 `minitcad solve` 顶层交付文件；固定三文件契约仅约束一次**正式仿真 solve**，不禁止开发仓库含测试、周报或课程文件。

### 14.5 Week 02 可复现环境和 VS Code 操作

推荐在 VS Code 中使用 **File → Open Folder** 打开项目根目录，安装 Python 和 Pylance 扩展，在命令面板执行 `Python: Select Interpreter` 选择项目虚拟环境。Markdown 预览：Windows/Linux `Ctrl+Shift+V`，macOS `Cmd+Shift+V`；测试可通过 Testing 面板或终端运行。

```bash
# macOS / Linux；所有命令在项目根目录执行
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[dev]"
python -m pytest -q tests/test_bernoulli.py
python -m pytest -q tests/test_poisson.py
python -m pytest -q
python -m examples.verify_poisson
```

```powershell
# Windows PowerShell；若执行策略限制激活，可直接调用 .venv\Scripts\python.exe
py -3 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -e ".[dev]"
python -m pytest -q
python -m examples.verify_poisson
```

`pyproject.toml` 至少应声明 `requires-python = ">=3.10"`、NumPy、SciPy 及可安装的 `dev` 依赖 `pytest`，并按 `src/` 包布局配置 setuptools 的发现路径；`.vscode/settings.json` 至少打开 `python.testing.pytestEnabled` 并设置 `python.testing.pytestArgs=["tests"]`。需要执行项目任务时，VS Code `tasks.json` 应使用当前选中的 Python 解释器，不能在不同开发者机器上硬编码绝对路径。

**执行环境的真实输出由运行者填写。**未经实际运行的通过数量、日志和截图不得写入本项目 Gate 报告。

### 14.6 Week 02 测试矩阵与准入标准

| 测试组 | 输入 / 独立判据 | 建议验收 |
|---|---|---|
| B-01 正常与小量 | `z=0, ±1e−14, ±1e−10, ±1e−3` 两侧，`±1, ±50` 两侧 | 高精度参照：可表示范围约 `3e−14` 相对误差 |
| B-02 大量和恒等式 | `±700, ±1000, ±float64_max`，`B(−z)=B(z)+z` | 所有有限输入结果有限、非负；允许 `B(1000)=0` 正确下溢 |
| B-03 API 与异常 | 浮点标量、二维整数数组、空数组、列表、`NaN/Inf`、不可转换输入 | 标量 `float`、数组同形 `float64`、非法输入显式异常 |
| P-01 线性与边界 | `φ=a+bx/L`、`ρ=0`，刻意改变一个接触电势 | 内部残差接近舍入误差，边界返回准确 Dirichlet 差 |
| P-02 二次与非均匀 | `φ=a+bx/L+c(x/L)^2`、`ρ=−2εc/L²`；任意严格递增网格 | 均匀/非均匀均为舍入误差量级；错误分母下标的回归样例必须失败 |
| P-03 正弦与网格加密 | `φ=sin(πx/L)`，常数 `ε`、均匀网格、`N=17,33,65,129` | 最大**解误差**之比在约 `[3.8,4.2]`，保留实际值 |
| P-04 单位等价 | 先按 cm 直接归一化，并对同物理设备计算 SI 教学 `F_P^SI` 与生产 `F_u^nd` | 内点 `F_P^SI/(−q n_i[m^-3])` 与 `F_u^nd` 一致；边界另按 `V_T` 缩放；另测 cm↔SI 往返及 cm 输入/经 SI 适配的一致性 |
| P-05 非法配置 | `N<3`、重复/倒序节点、非正 `ε`、错误形状、非有限电荷和 BC | 每项显式失败，不出现可接受的伪结果 |
| I-01 A/B 交接 | Mesh/State 不可原地改写，3M 交错 pack/unpack 往返 | 值、形状、复制语义正确；接口表 A/B 签字 |

上述阈值是**课程/门槛目标**，不应宣称所有机器上绝对残差必定单调下降。采用不同物理尺度时应优先报告同量纲比较和可解释的舍入下限。**任意非常不规则的非均匀网格不自动享有均匀网格的同阶收敛声明**；本表二阶比值专用于均匀正弦解。

### 14.7 W2 → W3 的接口会审记录（双方签字）

| 字段 | 正式生产版本 | 教学 SI 验证版本 | Owner / 核对方式 |
|---|---|---|---|
| 网格 | `x_nodes_cm[M]` 及 `x_nodes_nd[M]`,`h_edges_nd[M−1]`,`w_interior_nd[M−2]`（cm 物理网格经尺度归一化） | `x_m[M]` 单位 m，`h=diff(x)`（只在 SI 教学/适配） | B 定义、A 校核严格递增与双向换算 |
| 势 | `State.u[M]` 无量纲 | `phi_V[M]`、`phi_L/R_V` | A/B 核对 `φ=V_Tu` |
| 载流子和掺杂 | 物理输入 `n_cm3,p_cm3,ND_cm3,NA_cm3`，内核 `n_nd,p_nd,Cbar_nd`；内部 CV 掺杂平均 | SI 教学 `n,p,N_D,N_A` 用 `m−3`，`rho_C_m3=q(p−n+N_D−N_A)` | A 核对符号、B 核对 cm↔SI 适配 |
| 介电常数 | `epsilon_F_cm>0` 用于 `L0_cm`（生产无量纲残差不传 ε） | `epsilon_F_m>0`，可选节点数组用于拓展实验 | B 核对 F/cm↔F/m，A 负责 SI MMS |
| 状态布局 | `X=[u0,n0,p0,...]`，无量纲 | 仅供布局练习的 `[phi0,n0,p0,...]`，SI | B 保管 `State`，A 的 pack/unpack 不另创第二份模型 |
| Poisson 行 | `F_u^nd`，所有行无量纲 | 内部 `C/m³`，边界 `V` | A 提供，两套间用第 3.6 节解析关系比较 |
| 非法状态 | 非正 `n/p`、非有限数据拒绝 | 负 ε、非法 mesh/ρ/边界拒绝 | 双方交换反例，失败有可读原因 |
| Week03 下游 | `bernoulli → SG → continuity → F(X)` | SI MMS 验证器**不参与**耦合求解 | A 主集成，B 独立测试 |

建议在 `docs/a_b_interface_week02.md` 和 `docs/numerical_conventions.md` 记录接口版本、批准人、变更日期；如果团队已有 Week01 规范，直接合并，不另建相互矛盾的第二套契约。

**Week02 周报模板**（仓库内部，不计入 solve 的三文件）：

```markdown
# Week 02 验收记录

- A / B 审核：
- Git commit、Python/NumPy/SciPy/pytest 版本、运行平台：
- 输入参数、T 与单位约定：
- 4 项主文件状态：
- 测试命令及实际结果：
- Bernoulli 独立精度/恒等式最大误差：
- 直线/二次/正弦制造解：节点数、误差、误差比：
- cm↔无量纲、cm↔SI 双向往返及 SI 教学 Poisson 对照误差:
- 非法输入、故障注入和诊断：
- 未通过项、Owner、Week03 阻塞项：
- A 签字 / B 签字：
```

验收者须能够**口头独立推导** Bernoulli 的反射恒等式、Poisson 的控制体积分及其符号；讲清 `C/m³` 与 `V` 的混合单位；指出为什么二次制造解不能单独证明二阶精度。否则只运行绿灯测试不视为知识层面完成。

***


### 14.8 W8：A 主责 GUI、B 主责接口联调的五日执行表

> **前置条件**：`w07-core-freeze` 已包含 B 经测试冻结的配置校验、API、三文件 writer、合法温和掺杂 YAML，以及 A 可复用的 solver/RunResult。随附源码包的 `minitcad/gui/`、`apps/minitcad_gui.py` 和 `docs/GUI_USER_GUIDE.md` 可直接参考和适配，不要求为了“重写”而丢掉已有正确实现。

| 日期 | A：必须交付的 GUI 代码与演示 | B：接口与独立验证任务 | 当日退出证据 |
|---|---|---|---|
| **Day 1：GUI 启动/框架** | `gui/app.py`、`main_window.py`；主窗口展示左 P 右 N、右接地、正向 `V_app>0`；读取内置 YAML，菜单/分页骨架，窗口可正常关闭 | 再次锁定 `load_config/parse_config_dict/save_config`、`run(progress,should_stop)`、`RunResult` 与 writer；提供可预测的合法/非法基准 | Qt `--smoke-test`、源码启动及方向/单位提示截图；接口签名对照表 |
| **Day 2：输入编辑/校验** | `form_model.py`、默认 cm 参数表（SI 输入模式走 B 的适配）、求解器选择、导入/保存 YAML、高级 YAML 验证并应用；多分区/非均匀配置不被简化面板悄悄覆盖 | T16：字段往返等价、负温度/错误 PN 极性/材料温度不一致、无效输出目录均由统一校验拒绝 | GUI 与 CLI 导入同配置得到相同规范化 YAML/SHA；无效配置无法进入求解 |
| **Day 3：后台运行/取消** | `worker.py`，`QThread + QObject + Signal`；运行期间窗口响应；显示真实偏压尝试/状态；`threading.Event` 协作式取消；可定位 Python/Qt 异常 | T17：成功、真实 solver 失败、失败点回滚、用户取消；检查 `CANCELLED/REJECTED` 和三个文件不丢失已接受点 | 计算不阻塞 UI；工作线程不得直接修改 QWidget；取消后报告有真实失败原因 |
| **Day 4：展示/打包** | `plots.py`：真实 I–V、代表性剖面、失败点表、验收报告页；打开输出目录；编写 GUI 手册并修复适配缺陷 | 比对 GUI/CLI 偏压、`RunResult`、CSV/PNG/report 一致性；运行所有 T01–T17；在 Windows 执行 onedir 打包及 EXE 烟雾测试首轮 | GUI 截图能回溯 CSV；单次仿真输出目录恰好三件固定文件；打包日志与缺陷列表 |
| **Day 5：完整验收/发布** | 修复 GUI/高优先接口缺陷；从干净 Python 环境启动 GUI；交付交互演示、最简使用说明与 GUI tag | 在 Windows 10/11 验证 EXE/打包目录完整性，真实 0→0.30 V 示例、负例、取消；审计复现信息、合并 release tag，给出最终判定 | `w08-release`、T16–T18 测试证据、用户手册、三文件、双人签收；若缺 Windows 实测，EXE 状态标 `NOT_VERIFIED`，不冒充通过 |

**GUI 的最小强制界面**：顶部方向/极性提示与导入/保存/运行/停止/打开结果；左侧“参数设置＋高级 YAML”；右侧“I–V & 扫压＋物理剖面＋验收报告”；状态/日志区明确显示 `READY/RUNNING/ACCEPTED/REJECTED/ERROR`。允许更美观的视觉设计，但不得减少输入校验与科学证据。

**线程与绘图约束**：计算只能在工作线程中执行，使用 Signal 把偏压点、完成和失败信息传回主线程；`should_stop()` 是安全点检查而非线程强杀。GUI 的百分比不可基于未知长度的 adaptive sweep 虚构；可用不确定进度条和真实偏压表。只对 `RunResult.successful` 画物理曲线，**不能把失败尝试绘为已接受 I–V 点**。`convergence.md` 中的最终验收判定高于窗口“画出来了”的表象。

### 14.9 A/B 接口冻结、Windows 构建与 VS Code 工作流

| 时点 | B 必须向 A 提供 | A 必须向 B 提供 | 不一致时的处理 |
|---|---|---|---|
| **W2 末** | `Mesh1D/State/NormalizedDevice`、contact、正式 cm ↔ nd 契约与 SI↔cm 双向适配 | `bernoulli`、Poisson、独立单位映射验证 | 双方签名冻结；不得在 W3 无通知改字段 |
| **W4 末** | `RunResult` 最小结果对象、固定扫描/错误码骨架 | `residual/edge_flux/equilibrium/Gummel` | T08/T09 双方交叉验证；无可复现 0 V 解不进入 W5 |
| **W6 末** | `run(config,progress,should_stop)`、三文件 writer 初版、schema/失败回滚 | Newton/Gummel、Jacobian、完整 `SolveResult` 与诊断 | T10/T11/T13 未通过直接标记 Gate 2 失败 |
| **W7 末** | 完整配置读写、已冻结 `RunResult`/`Attempt`、稳定输出路径和回调、文档/样例 | 归一化、Poisson/SG/SRH/continuity、Equilibrium/Gummel/Newton、Jacobian、阻尼、强偏压、非均匀网格、性能与物理终审均通过，并提供可供 GUI 直接调用的稳定内核 | `w07-core-freeze`；W8 不做破坏性 API 改名 |
| **W8 末** | 独立 T16–T18、Windows 打包和端到端对照结论 | 完整 PySide6 GUI、后台计算/取消、真实结果页和 GUI 文档 | 双签后发布；核心科学失败或必需桌面功能缺失则 `REJECTED` |

**在 VS Code 中统一使用本文件规定的 Markdown 与 Python 工作流**：文件采用 UTF-8、标准 `.md` 标题/表格、围栏代码块和 `$$...$$` 数学公式；用 `Ctrl+Shift+V`（macOS `Cmd+Shift+V`）预览，不转换成 HTML 或 Word。继续使用第 14.5 节的 `.venv`、pytest、`.vscode/settings.json`/`tasks.json`。W8 GUI 依赖作为**可选项**安装，核心模块不应无条件导入 Qt：

```powershell
# Windows / PowerShell：在项目根目录；建议 Python 3.11 x64
py -3.11 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -e ".[gui,dev]"
.\.venv\Scripts\python.exe -m minitcad.gui.app
.\.venv\Scripts\python.exe -m pytest -q tests/test_gui_qt_smoke.py

# Windows EXE：需在 Windows 10/11 x64 本机构建、测试
.\scripts\build_windows.bat
```

```bash
# macOS / Linux：仅验源码 GUI，不据此宣称 Windows EXE 已通过
python3 -m venv .venv
./.venv/bin/python -m pip install -e '.[gui,dev]'
./.venv/bin/python -m minitcad.gui.app
./.venv/bin/python -m pytest -q tests/test_gui_qt_smoke.py
```

GUI/打包参考入口与随附源码包一致：`src/minitcad/gui/app.py`、`apps/minitcad_gui.py`、`scripts/build_windows.bat`（推荐 onedir）；单文件 `scripts/build_windows_onefile.bat` 仅为备选。发布 onedir 时必须复制整个 `dist/MiniTCAD/` 目录，不能只复制 `MiniTCAD.exe`。若环境不满足 PySide6/Windows 前提，必须在日志中标明缺失的测试，而不能把 `pytest.importorskip` 的跳过记录算作 GUI 合格。

### 14.10 8 周排期的风险登记与最低签核文件

- **A 的关键路径 W3–W7**：SG/SRH/残差 → equilibrium/Gummel → Newton/Jacobian → 强偏压/非均匀网格终审。每次失败保留 YAML、commit、最大残差、step、线搜索原因和修复前后对照。W7 **不得承诺“后来 W8 补算”**。
- **B 的关键路径 W5–W7**：独立 Jacobian/回滚 → 网格参考 → 三文件完整验收和 API 冻结。B 若到 W7 末还没有可靠 `RunResult`/writer，A 的 GUI 会被阻塞，必须在 Gate 报告中明确责任归属。
- **高掺杂压力场景**：可用于触发算法弱点，不得超越非简并/常迁移率等本期物理假设而声称物理结论普适；对于事前定义属于接受域的工况，求解失败应如实记为失败。
- **W8 时间冲突**：优先保留 GUI 的有效输入、后台求解、正确取消、真实结果和固定三文件；缩减装饰性动画、主题、额外导出，不可砍掉科学验收。
- **签核最低材料**：`reports/week04_gate1.md`、`reports/week06_gate2.md`、`reports/week07_core_freeze.md`、`reports/week08_gate3.md` 为**仓库内部验收记录**，附命令、日志、Git SHA、配置 SHA、测试实际值/阈值/结果、A/B 签字。它们**不属于**单次 `solve` 输出目录的三文件契约。

---

## 15. 最终计算实验与评分维度

| 实验 | 保持不变 | 只改变 | 必须量化的结果 |
|---|---|---|---|
| A. Gummel vs Newton | Device、mesh、SRH、bias 列表、容差 | nonlinear solver | 每偏压迭代数、用时、最终残差、成功率、终端电流差 |
| B. Jacobian 方法 | 同 `F(X)`、采样状态、矩阵结构 | FD（独立检查）、稀疏 FD，若有可再测 AD/complex-step | 方向导数误差、生成时间、Newton 步数；注明方法可微性边界 |
| C. Uniform vs junction-refined | Device、偏压、材料、总节点数（对比组） | 节点空间分布 | `I`、势垒、峰值 `E`、与更细参考的误差、总用时 |
| D. Fixed vs adaptive sweep | 初始平衡态、目标终点、solver/容差 | 电压调步策略 | 成功点、失败重试、总迭代数和耗时、完整偏压覆盖 |
| E. Verification & reproducibility | 预先登记的测试集 | solver/mesh 的交叉设置 | T01–T15 通过率、反例和失败诊断、三文件验收复现 | 
| **F. W8 GUI ↔ CLI 等价与桌面交付** | 相同规范化 YAML、计算核心与物理验收规则 | GUI/CLI 入口及 GUI 成功/失败/取消路径 | **T16–T18**：一致性、窗口响应、真实数据展示、取消回滚、Windows 启动/打包记录 |

“像一条二极管 I–V 曲线”不是证据。最终答辩/报告应按 **Physical model → Normalization → FV/SG → Residual → Nonlinear solver → Sweep → Verification → Three deliverables → GUI/API integration（W8）** 展示；每个图中的理论/参考含义都能追溯到 CSV 和 `convergence.md`。

***

## 16. 发布前最终检查表（A/B 双签）

### A：物理与数值完整性

- [ ] 电子/空穴电流与连续性方程符号正确，`J_n+J_p` 守恒。
- [ ] SRH 为 `τ_p(n+1)+τ_n(p+1)`；不同寿命测试 `19/47` 通过。
- [ ] 非均匀 Poisson 用 `h_{i-1}+h_i`；二次函数制造解通过。
- [ ] SG 电子、空穴两式符合 `B(z)` 对称性与热平衡零通量。
- [ ] cm 工程参数直接归一化并恢复一致；SI↔cm 往返验证通过；必要的 SI 文件仅在 writer 转换；Jacobian 对唯一残差求导。
- [ ] Gummel/Newton 的 BC、SRH、当前偏压一致；无负载流子或 silent clipping。
- [ ] 最终收敛包含 scaled residual + scaled update 两项；返回迭代历史。

### B：结构与用户验收完整性（W7 基础冻结、W8 联调发布）

- [ ] YAML/API/CLI/PNG/report 都提醒“左 P 右 N，右 0 V，正向左端电压为正”。
- [ ] `type: pn` 对左右掺杂类型、单位、区域覆盖进行严格校验。
- [ ] 三个固定文件名称、默认 cm/可选 SI schema、物理单位字段与一致的 `run_id`/SHA 已自动验证；SI 输出仅作 writer 单位转换。
- [ ] `comparison.csv` 含全部节点/边/IV/迭代，以及有依据的理论值、误差和数值参考。
- [ ] 对没有严格理论解的量写 `N/A`，零值理论时不计算相对误差。
- [ ] `comparison.png` 同时包含数值、理论（如有）、误差、收敛及明确偏压方向。
- [ ] `convergence.md` 按每个尝试偏压记录迭代、容差、失败原因、网格与全部复现设置。
- [ ] 独立 Jacobian、守恒、网格验证与回滚测试覆盖，不是由 A 的实现自证。
- [ ] 用相同配置从零重复 `solve`，输出比较结果能再现。
- [ ] W7 末 `load_config/parse_config_dict/save_config/run(progress,should_stop)/write_required_outputs`、`RunResult` 的契约和版本已冻结。
- [ ] T16/T17 验证 GUI 与 CLI 的 YAML、真实计算、进度、失败、取消、三文件一致；B 不能只接受 A 的 GUI 截图。
- [ ] Windows 发布需要 T18：打包后的 onedir EXE 在目标平台确实启动和计算；未测试则明确标 `NOT_VERIFIED`。

### A：第 8 周 GUI 与接口对齐（新增必签）

- [ ] GUI 能在 PySide6 桌面窗口运行；顶部固定呈现左 P 右 N、右端接地、`V_app>0` 为正向、默认 cm 输入与显式单位；可选 SI 输入须通过 B 的适配层。
- [ ] 参数设置 / 高级 YAML / 导入导出遵循 B 的配置校验，不擅自忽略多区域或非均匀网格。
- [ ] 后台 `QThread` 调用 B 的 `run(progress,should_stop)`；GUI 保持响应，取消保留已接受状态与 `REJECTED` 三文件。
- [ ] I–V、物理剖面和验收报告只消费真实 `RunResult` 与文件；失败偏压在列表里明确显示，不伪装成功曲线。
- [ ] 已与 B 完成 T16/T17 接口交叉验收，交付源码入口、GUI 手册、Qt 冒烟测试和复现方式。

**最终验收判断**：任何硬性物理恒等式、三文件 schema、A 的 W8 GUI 必需功能或 B 的端到端接口验收失败，结论为 `REJECTED`，必须列出失败项及责任人。**数值核心、GUI 源码运行与 Windows EXE 发布分开记录实测状态**；未经 Windows T18 实测，不能声称 Windows EXE 已验收。本版允许不添加 2D 或高级迁移率模型，但**不再允许缺省 GUI**。

***

## 17. Week 02 可复制的参考实现与数值复核

以下代码为第 2 周模块实现和单元测试提供最小参考。它们**不是已经完成的 PN 耦合求解器**，也不替代第 3 周的统一 `residual()`、SRH、SG 或 B 的 `Mesh/State` 封装。生产核心必须复用第 10 节定义的网格、状态和 Dirichlet 契约。这里给出 Bernoulli、正式无量纲 Poisson 以及用于独立交叉验证的 SI Poisson MMS 函数。

### 17.1 `src/minitcad/discretization/bernoulli.py`

```python
"""Stable Bernoulli function for Scharfetter-Gummel discretization."""
from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


def bernoulli(z: ArrayLike) -> float | NDArray[np.float64]:
    """B(z)=z/(exp(z)-1), B(0)=1. Finite real inputs only.

    Return float for a scalar; a same-shape float64 array otherwise.
    Do not alter the input.
    """
    x = np.asarray(z, dtype=np.float64)
    if not np.all(np.isfinite(x)):
        raise ValueError("Bernoulli input must be finite")

    out = np.empty_like(x, dtype=np.float64)
    small = np.abs(x) < 1.0e-3
    pos = x > 50.0
    neg = x < -50.0
    middle = ~(small | pos | neg)

    t = x[small]
    t2 = t * t
    out[small] = 1.0 - t / 2.0 + t2 / 12.0 - t2 * t2 / 720.0 + t2**3 / 30240.0

    t = x[middle]
    out[middle] = t / np.expm1(t)

    t = x[pos]
    ex = np.exp(-t)
    out[pos] = t * ex / (1.0 - ex)

    t = x[neg]
    ex = np.exp(t)
    out[neg] = -t / (1.0 - ex)

    return float(out) if x.ndim == 0 else out
```

### 17.2 `src/minitcad/discretization/poisson.py`（正式无量纲 Poisson）

```python
"""Production Poisson-only residual in dimensionless silicon coordinates.

Interior F_u=(gradient_right-gradient_left)/CV+Cbar+p-n.
Boundary F_u=u-u_target. All rows are dimensionless.
"""
from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


def _vector(name: str, data: ArrayLike, size: int | None = None) -> NDArray[np.float64]:
    a = np.asarray(data, dtype=np.float64)
    if a.ndim != 1 or (size is not None and a.size != size):
        raise ValueError(f"{name} must be a 1D vector of length {size}")
    if not np.all(np.isfinite(a)):
        raise ValueError(f"{name} must contain finite values")
    return a


def poisson_residual_nd(
    xi: ArrayLike,
    u: ArrayLike,
    n: ArrayLike,
    p: ArrayLike,
    cbar: ArrayLike,
    u_left: float,
    u_right: float,
) -> NDArray[np.float64]:
    """Return M Poisson residual rows; cbar is the CV-averaged net doping."""
    xx = _vector("xi", xi)
    m = xx.size
    if m < 3:
        raise ValueError("A mesh requires at least three nodes")
    uu = _vector("u", u, m)
    nn = _vector("n", n, m)
    pp = _vector("p", p, m)
    cc = _vector("cbar", cbar, m)
    if np.any(nn <= 0.0) or np.any(pp <= 0.0):
        raise ValueError("Carrier concentrations must be positive")
    bc = np.asarray([u_left, u_right], dtype=np.float64)
    if not np.all(np.isfinite(bc)):
        raise ValueError("Boundary potentials must be finite")

    h = np.diff(xx)
    if np.any(h <= 0.0):
        raise ValueError("xi must be strictly increasing")
    cv = 0.5 * (h[:-1] + h[1:])
    gradients = np.diff(uu) / h
    f = np.empty(m, dtype=np.float64)
    f[0] = uu[0] - bc[0]
    f[-1] = uu[-1] - bc[1]
    f[1:-1] = np.diff(gradients) / cv + cc[1:-1] + pp[1:-1] - nn[1:-1]
    if not np.all(np.isfinite(f)):
        raise FloatingPointError("Poisson residual overflow: examine physical scaling")
    return f
```

**调用示例**：非均匀节点 `xi=[0, 0.25, 0.75, 1.5]`、制造势 `u=xi²`、`n=p=1`、`cbar=-2`、Dirichlet `0,2.25`，应给出近零残差。这里的 `cbar` **已是 B 在控制体积积分后提供的平均净掺杂**，不能把层间跃变点的一侧值直接塞进内部残差。

### 17.3 `examples/poisson_si_oracle.py`（**独立 SI 教学验证器**）

```python
"""Stand-alone SI Poisson oracle for constant permittivity.

This is NOT the production F(X). Interior rows have units C/m^3;
Dirichlet rows have units V. Do not take an unscaled common norm.
"""
from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


def poisson_residual_si_mms(
    x_m: ArrayLike,
    phi_V: ArrayLike,
    rho_C_m3: ArrayLike,
    epsilon_F_m: float,
    phi_left_V: float,
    phi_right_V: float,
) -> NDArray[np.float64]:
    x = np.asarray(x_m, dtype=np.float64)
    phi = np.asarray(phi_V, dtype=np.float64)
    rho = np.asarray(rho_C_m3, dtype=np.float64)
    eps = float(epsilon_F_m)
    bc = np.array([phi_left_V, phi_right_V], dtype=np.float64)
    if x.ndim != 1 or x.size < 3 or phi.shape != x.shape or rho.shape != x.shape:
        raise ValueError("x, phi, rho must be same-shape 1D arrays with >=3 nodes")
    if not np.isfinite(x).all() or not np.isfinite(phi).all() or not np.isfinite(rho).all():
        raise ValueError("x, phi, rho must be finite")
    if not np.isfinite(bc).all() or not np.isfinite(eps) or eps <= 0:
        raise ValueError("epsilon must be positive; boundary values must be finite")
    h = np.diff(x)
    if np.any(h <= 0):
        raise ValueError("x must be strictly increasing")
    flux = eps * np.diff(phi) / h
    cv = 0.5 * (h[:-1] + h[1:])
    out = np.empty(x.size, dtype=np.float64)
    out[0] = phi[0] - bc[0]
    out[-1] = phi[-1] - bc[1]
    out[1:-1] = -np.diff(flux) / cv - rho[1:-1]
    if not np.isfinite(out).all():
        raise FloatingPointError("Nonfinite SI Poisson residual")
    return out
```

**严格提醒**：SI 验证器内部结果为 `C/m³`，边界为 `V`；不可用其未缩放范数直接驱动 W3/W4 求解器。这里仅接受常数 `epsilon_F_m`，以便与生产 Si 同质模型完成严格的 `F_P^SI=-q n_i F_u^nd` 核对。变量 ε 的测试若另写，须放在独立扩展范围。

### 17.4 参考实现必须再次验证的最小执行序列

```bash
python -m pytest -q tests/test_bernoulli.py
python -m pytest -q tests/test_poisson.py
python -m examples.verify_poisson
python -m pytest -q
```

第 14.6 节的完整测试矩阵需要由项目团队落到 `tests/` 文件中；此处参考源码的单独运行不意味着上述四条命令已经在最终团队仓库通过。应尤其验证**成功路径**与**失败路径**，并报告本次真实测试命令、实测误差、所用依赖版本。

### 17.5 非耦合 Poisson 正弦制造解：可独立复算的误差关系

对长度 `L`、常数绝对介电常数 `ε`，取

$$
\phi_{\star}=\sin(\pi x/L),\quad
\rho_{\star}=\epsilon(\pi/L)^2\sin(\pi x/L),\quad \phi_\star(0)=\phi_\star(L)=0.
$$

均匀网格 `N` 节点的离散正弦特征向量满足：**数值势相对解析正弦的幅度倍数**为

$$
A_N=\frac{\pi^2}{4(N-1)^2\sin^2\!\big(\pi/[2(N-1)]\big)}.
$$

因此当解析解峰值恰为 1 V，最大电势幅值误差应与 `A_N−1` 一致；它提供了一个不依赖同一份 Poisson 代码的解析交叉检查。实验不能只展示离散 residual 几乎为零，还必须展示数值解对**连续解析解**的误差随网格加密缩小。对应的独立参考数值见附录 E.1。

---

## 附录 A. 关键手算/计算器核对样例（不依赖完整 TCAD）

1. **SRH 交换检查**：`n=10,p=2,τ̂_n=1,τ̂_p=4` ⇒ `r=(20−1)/[4(10+1)+1(2+1)]=19/47≈0.404255319148936`；旧错误顺序给 `19/23≈0.826086956521739`，必须被测试识别。
2. **非均匀 Poisson**：`ξ=(0,0.25,0.75,1.5)`、`u=ξ²`，在 `i=1` 左/右间距为 `0.25/0.50`，正确分母给 `u''=2.0`；错用右侧下一个间距 `0.75` 给 `1.2`。
3. **接触平衡**：对任意净掺杂 `C`，令 `a=asinh(C/2), n=e^a,p=e^{-a}`，理论上 `np=1`、`p−n+C=0`。左 P 要求 `a_L<0`，右 N 要求 `a_R>0`。
4. **SG 平衡**：给定 `Δu` 和 `n_R=n_L e^{Δu}`、`p_R=p_L e^{-Δu}`，用 `B(-z)=e^zB(z)` 可验证 `j_n=j_p=0`。
5. **端电流守恒**：将每个内部节点 `F_n+F_p=0`，得 `j_{i+1/2}=j_{i-1/2}`；这是**离散代数**检验，不需要依赖 I–V 曲线形状。

## 附录 B. 技术参考资料

- TU Wien（SG 电子与空穴边电流公式）：https://www.iue.tuwien.ac.at/pdf/ib_2010/BC2010_Triebl_1.pdf
- Sesame 文档（SG、空间 link/site 和离散残差）：https://sesame.readthedocs.io/en/latest/physics/discretization.html
- COMSOL Semiconductor Manual（SRH 公式及参数约定）：https://doc.comsol.com/5.3/doc/com.comsol.help.semicond/SemiconductorModuleUsersGuide.pdf
- COMSOL 6.3 `Metal Contacts`（理想欧姆接触电中性及局部平衡假设）：https://doc.comsol.com/6.3/doc/com.comsol.help.semicond/semicond_ug_semiconductor.6.63.html

> 技术参考资料用于核查物理方程和离散形式；项目接口、验收文件及通过标准由本文规定。


## 附录 C. 示例配置的 cm/无量纲数值核对（不是求解结果）

按第 13 节 **cm** YAML 的材料参数，以 `q=1.602176634×10^-19 C`、`k_B=1.380649×10^-23 J/K`、`ε_0=8.8541878128×10^-14 F/cm` 计算（仅为配置与解析尺度核对）：

| 尺度/接触基准 | 核算结果 | 物理解释 |
|---|---:|---|
| `V_T=k_B T/q` | `0.0258519998 V` | 300 K 热电势 |
| `L_0=sqrt(ε V_T/(q n_i))` | `4.08845452×10^-3 cm` | 按 **本征浓度** 定义的长度基准；相当于 `4.08845452×10^-5 m`，不是 PN 耗尽层宽度 |
| `t_0=L_0²/D_0` | `1.67154603×10^-6 s` | 参考扩散时间 |
| `J_0=qD_0n_i/L_0` | `3.91878307×10^-6 A/cm²` | 参考电流密度；相当于 `0.0391878307 A/m²` |
| `C_L`, `C_R` | `−1.0×10^6`, `+2.0×10^5` | 内核无量纲净掺杂，左 P 右 N |
| `a(C_L)`, `a(C_R)` | `−13.81551056`, `+12.20607265` | 0 V 接触势（无量纲） |
| `V_bi=V_T[a(C_R)−a(C_L)]` | `0.6727099634 V` | 0 V 时接触之间的内部静电势差 |
| `V_app=+0.6 V` 时的两端内部静电势差 | `0.0727099634 V` | 仅指指定 Dirichlet **两端势差**，不是求解得到的内部电场/电流 |

上述数值是从**cm 配置参数和解析定义**直接计算出来的，可用于 B 的 `T03_CONTACT` / cm↔SI 双向适配与 cm↔无量纲回归测试；不应误称其为完整 Poisson–DD 的数值解或可直接填到非平衡 `I(V)` 列。

---

## 附录 D. Si 参数的 cm 工程单位、SI 对照、温度适用性与参考尺度

### D.1 300 K 参考值不是通用默认值

此表列出常见 Si 工程参数量级；数值范围用于输入检查和数量级直觉，**不是**某个温度/偏压/掺杂区间以外的绝对禁用阈值。

| 物理量 | 工程参考（通常 300 K、低场或轻掺杂） | 对应 SI 表示 | 程序说明 |
|---|---|---|---|
| Si 相对介电常数 | `ε_r≈11.7` | `ε_Si≈1.03594e−10 F/m` | `ε_r` 无量纲，内核常 ε |
| 本征浓度 | `n_i≈(1.0–1.5)×10^10 cm^-3` | `≈(1.0–1.5)×10^16 m^-3` | 与当前 `T` 的明确材料参数绑定 |
| 电子迁移率 | `μ_n≈1350–1450 cm²/(V·s)` | `≈0.135–0.145 m²/(V·s)` | 每次运行常数 |
| 空穴迁移率 | `μ_p≈450–500 cm²/(V·s)` | `≈0.045–0.050 m²/(V·s)` | 每次运行常数 |
| 电子扩散系数 | `D_n≈35–37 cm²/s` | `≈(3.5–3.7)e−3 m²/s` | 以 `μ_n V_T` 为准 |
| 空穴扩散系数 | `D_p≈11–13 cm²/s` | `≈(1.1–1.3)e−3 m²/s` | 以 `μ_p V_T` 为准 |
| 典型掺杂跨度 | `10^14–10^20 cm^-3` | `10^20–10^26 m^-3` | 高掺杂不自动满足非简并假设 |
| 典型电场跨度 | `10^3–10^5 V/cm` | `10^5–10^7 V/m` | 实际电场应由自洽解导出 |
| 备选参考迁移率 | `μ_0=1000 cm²/(V·s)` | `0.1 m²/(V·s)` | 可取 `D_0(T)=μ_0V_T(T)` |

第 13 节 YAML 的 `μ_n=1350 cm²/(V·s)`、`μ_p=480 cm²/(V·s)` 分别对应 `0.135` 和 `0.048 m²/(V·s)`，符合上述数量级；**这不证明**常迁移率模型对任意高电场都正确。

### D.2 同一器件参数下，两种 D0 选法不能混淆

按 `T=300 K`、`n_i=1.0e10 cm^-3`、`ε_r=11.7`，采用一致的 cm 尺度独立计算：

| 项目 | 数值 | 含义 |
|---|---:|---|
| 热电压 `V_T` | `0.02585199979 V` | 随 `T` 更新 |
| 本征参考长度 `L_0` | `4.08845452e−3 cm` | 即 `40.88454518 μm`，与附录 C 一致 |
| 若取 `μ_0=1000 cm²/(V·s)`，`D_0=μ_0V_T` | `25.85199979 cm²/s` | 等价于 `2.585199979e−3 m²/s`，是备选参考 |
| 第 13 节 YAML 中 `D_0=10.0 cm²/s` 对应的 `μ_0` | `386.8172707 cm²/(V·s)` | 等价于 `0.0386817271 m²/(V·s)` |

两种 D0 选择不是两种物理扩散系数模型：若所有无量纲系数正确同步转换，cm 恢复（或经 writer 输出 SI）的同一物理问题应给出相同物理解；若改变参考 D0 后 I–V 变了，应排查尺度换算和程序一致性。

### D.3 以不同浓度定义的 Debye 长度

同为 300 K 和 `ε_r=11.7`，由 `sqrt(ε V_T/(q N))` 可得：

| 参考浓度 | 估算长度 | 用途 |
|---:|---:|---|
| `n_i=10^10 cm^-3` | `≈40.88 μm` | 本规范基准 `L_0` |
| `N=10^15 cm^-3` | `≈129.29 nm` | 局部屏蔽长度 |
| `N=10^16 cm^-3` | `≈40.88 nm` | 局部屏蔽长度 |
| `N=10^17 cm^-3` | `≈12.93 nm` | 局部屏蔽长度 |
| `N=10^18 cm^-3` | `≈4.09 nm` | 局部屏蔽长度 |

以上不等同于 PN 耗尽宽度；是否捕捉结区空间变化要看实际物理网格，不仅看 `L/L_0`。可在 `convergence.md` 记录结区最小 `h_cm` 和与同单位局部尺度的比值作为网格诊断，但不能单凭一项比值宣布完整方程通过。

---

## 附录 E. Week 02 制造解参考数据及错误边界

### E.1 均匀网格正弦 MMS（**仅验证 Poisson，不是 PN 完整求解**）

针对长度 `L=1 μm`、常数 `ε`、Dirichlet 两端为零的线性 Poisson 系统，计算得到：

| 节点 N | $\max\lvert\phi_{\mathrm{num}}-\phi_{\mathrm{exact}}\rvert$ (V) | 误差比（上一网格 / 当前网格） |
|---:|---:|---:|
| 17 | `3.218964440e−3` | — |
| 33 | `8.035776794e−4` | `4.005791` |
| 65 | `2.008218097e−4` | `4.001446` |
| 129 | `5.020091592e−5` | `4.000361` |

第 17.5 节的解析幅度因子 `A_N−1` 可以独立复算这些误差。Gate 报告需填写**本团队运行的实际结果**，而不能直接把此参考表当作其工作站/器件程序已通过的测试记录。

### E.2 哪些现象是假通过？

1. 与生产代码共享同一错误实现的“理论对照”：采用解析恒等式、独立数值参考或 SI oracle 排除同错相消。
2. 二次多项式残差接近零即宣称一般问题二阶收敛：三点离散对二次函数恰好精确，还需正弦等非二次解的加密实验。
3. 极大正 `z` 时 `B(z)` 下溢为零被视为 bug：它可以是正确的 `float64` 表示，不应进行物理通量剪裁。
4. 将 SI 的内部 `C/m³` 残差与边界 `V` 残差直接作一个范数：只有第 3.6 节转换后的无量纲行才能成为生产收敛判据。

### E.3 参考函数的独立校验清单

验收时应实测第 17 节的参考函数：Bernoulli 覆盖零点、`±1e−14`、`±1e−3` 与 `±50` 分段交界、`±700`、`±1000`、标量/数组、非有限异常以及高精度 Decimal 对照；Poisson 覆盖非均匀二次制造解、非法网格、SI↔无量纲残差比例、均匀正弦 MMS 误差比。将**实际运行**的输入、误差、依赖版本、通过/失败情况写入验收记录。参考数据不能替代完整 TCAD、Gummel/Newton 或 Gate 的实际测试。