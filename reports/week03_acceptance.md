# W3 验收记录（B 侧）：QF 同步与 T01–T07 独立 oracle

> 对应 `目的与规划.md` §14.2 的 W3-B1…W3-B4。
> 日期：2026-10-25；解释器 Python 3.14.6；NumPy 2.5.3；SciPy 1.18.0。
> 同步基线：A 的 `d53871e`（2026-10-21，W3 物理/离散残差）与
> `b195ce5`（2026-10-23，QF 精度重构 + 求解器/状态接口）。
> 本轮把 A 的 W2/W3 物理模块并入根 `src/minitcad` 单一包，并把 B 的 oracle
> 从旧的 `(u,n,p)` 口径改到 A 的 QF 口径。

## 1. 复现命令与结果

```bash
cd <repo-root>
python -m pytest -q
# 235 passed in 1.15s（0 skipped）

python -m pytest tests/test_oracle_bernoulli.py tests/test_oracle_srh.py \
    tests/test_oracle_poisson.py tests/test_oracle_sg.py \
    tests/test_oracle_qf_residual.py -q
# 104 passed
```

A 的 `C2_v1` 独立树（仅作对照，不被根包导入）：

```bash
cd C2_v1 && python -m pytest -q
# 4 failed, 164 passed
# 失败均为 test_qf_stability.py 调用旧 residual_qf 签名（见 §5）
```

## 2. W3-B 任务状态

| 任务 | 产物 | 状态 / 证据 |
|---|---|---|
| W3-B1 T01–T07 独立 oracle | `tests/test_oracle_bernoulli.py`、`test_oracle_srh.py`、`test_oracle_poisson.py`、`test_oracle_sg.py`、`test_oracle_qf_residual.py` | 完成；不再 skip |
| W3-B2 反向 PN/非法配置 | `tests/test_core_contact.py`、`tests/test_io_config.py` | 完成；`INVALID_PN_ORIENTATION` 带两端 `C` 值 |
| W3-B3 `RunResult` + CLI 骨架 | `src/minitcad/simulation/run_result.py`、`src/minitcad/cli.py` | 完成；`validate` 真跑，`solve` 诚实返回 3（W4-B4 接线） |
| W3-B4 CSV long-form schema | `src/minitcad/io/results_csv.py` | 完成；22 列冻结、缺失值空字段 |

## 3. 本轮同步内容

### 3.1 A 模块并入根包（单一包，§10.1）

| A 交付 | 根包位置 | 说明 |
|---|---|---|
| Bernoulli | `src/minitcad/discretization/bernoulli.py` | 四支路稳定实现 |
| Poisson | `src/minitcad/discretization/poisson.py` | 支持 `cbar` 长度 M 或 M-2 |
| SG | `src/minitcad/discretization/sg.py` | 保留物理 `electron/hole_edge_flux`，新增 `edge_flux_qf` |
| continuity | `src/minitcad/discretization/continuity.py` | 两条连续性残差 |
| 唯一残差 | `src/minitcad/discretization/residual.py` | `evaluate_qf`/`residual_qf` 为 QF 生产入口；旧 `residual` 仅为历史接口 |
| SRH | `src/minitcad/physics/recombination.py` | `srh_recombination_nd` + `srh_recombination_qf`（`expm1`） |
| 平衡初值 | `src/minitcad/simulation/equilibrium.py` | 0 V 标量求解，7 次迭代收敛 |
| `__init__` | `discretization/`、`physics/`、`simulation/` | 重新导出；`physics/__init__` 不含 A 分叉的 normalization 辅助函数 |

`C2_v1/` 保留为 A 的参考树与回归基准，根包不导入它。

### 3.2 B 的 W2 数据结构补齐（QF）

- `core/state.py`：`QFState`、`PhysicsDiagnostics`、`SolveResult.physics`。
- `core/contact.py`：`contact_targets_qf`（左端 `eta_n=eta_p=-v`，右端 0）。
- `tests/test_core_qf_state.py`：12 项（往返、只读、3M、接触、诊断、失败路径）。

### 3.3 B 的 W3 oracle 改写

- T01 Bernoulli：修正了旧测试的反射恒等式符号（应为 `B(-z)-B(z)-z=0`）。
- T02 SRH：异寿命 `19/47`、`np=1` 零速率、符号，全部通过。
- T05 Poisson：三套非均匀网格的二次制造解，内部残差 ≤ `1e-11`。
- T06/T07 SG：
  - 物理口径 `electron_edge_flux/hole_edge_flux`：平衡零电流、纯扩散闭式、迁移率线性；
  - QF 口径 `edge_flux_qf`：常数 `eta_n/eta_p` 时电流严格为 `0.0`；
  - QF 与物理路径在非平衡状态下的边通量相对误差 ≤ `1e-12`。
- 新增 QF 残差集成 oracle（`test_oracle_qf_residual.py`，6 项）：
  - QF 内部行与旧 `residual` 内部行一致（Poisson ≤ `1e-13`、continuity ≤ `1e-11`）；
  - 接触行设成 `contact_targets_qf` 后严格为 0；
  - `evaluate_qf` 的 flux/rate 与独立函数一致；
  - 平衡标量解在 QF 残差下零电流、`np=1`；
  - 精度证据：`u=40`、`eta_n=1e-16` 时物理 `np-1` 已为 0，QF 仍给出正速率（`expm1`）。

## 4. 数值门槛对照

| 验收项 | 判据 | 结果 |
|---|---|---|
| T01 Bernoulli | 相对误差 ≤ 3e-14（自测更严） | 通过 |
| T02 SRH 异寿命 | `19/47`，误差 ≤ 1e-12 | 通过 |
| T03/T04 接触/极性 | `np=1`、`p-n+C=0`、`u_L=a(C_L)+v` ≤ 1e-10 | 通过 |
| T05 非均匀 Poisson | 二次制造解内部残差 ≤ 1e-11 | 通过 |
| T06 平衡零电流 | `j_n=j_p=0`（物理与 QF 两口径） | 通过（QF 严格 0） |
| T07 纯扩散 | `j_n=m_n(n_R-n_L)/h`、`j_p=m_p(p_L-p_R)/h` | 通过 |
| I-01（QF 扩展） | 只读、3M 交错、pack/unpack 往返 | 通过 |

## 5. 遗留问题（需 A 确认 / 后续周处理）

1. **`C2_v1` 4 项失败**：`tests/test_qf_stability.py` 仍调用旧
   `residual_qf(device, bias, u, eta_n, eta_p, ref_n=..., ref_p=...)` 并期望
   `(f, flux, rate)`；当前实现为 `residual_qf(QFState, device, bias) -> ndarray`。
   根包测试不依赖它；建议 A 更新测试并在 `C2_v1` 内 `164 passed` 收口。
2. **`solvers/` 未并入根包**：A 的 `gummel.py` 属 W4-A2，且同一文件残留两个
   `solve_gummel` 定义（旧版约 419 行、QF 版约 567 行，后者覆盖前者）。
   W4 合并前必须删旧留新，并由 A 确认 `ref_n/ref_p` 锚定右接触的契约。
3. **`normalization.py` 分叉**：A 的 `C2_v1` 版把 `e0_v_per_cm` 改名 `e0_v_cm`、
   删除 `qr0_a_cm3`、新增 `normalize_cm/denormalize_cm`；根包仍是 B 的 W1 版。
   需确认正式版本后再统一，当前根包测试均通过。
4. **`io/units.py` 分叉**：A 的 `C2_v1` 简化版删除了 `UnitError`/后缀校验；
   根包保留 B 的 W2-B4 严格版，不合并 A 的简化版。
5. **`device/mesh` 校验增强**（A 版）：`area_cm2` 可空、右端 0 V、`segmented`
   必须含结区断点；建议 W4 前由 A 确认后合并。
6. **`solve` 仍未接线**：W3-B3 只交骨架，`minitcad solve` 返回 3 且不写结果文件；
   W4-B4 才接 `run()`。

## 6. 结论

- B 的 W2/W3 交付在 QF 重构后全部通过（`235 passed, 0 skipped`）。
- A 的 W2/W3 物理模块已同步进根包，B 的独立 oracle 已能在根包上直接运行。
- W3 可以进入 W4（Equilibrium 正式版 / Gummel / sweep）的前提是 A 确认
  §5 的 `ref_n/ref_p` 与 `solvers/` 合并口径；根包当前不包含 `solvers/`。
