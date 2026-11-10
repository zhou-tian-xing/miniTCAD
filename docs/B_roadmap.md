# B 侧路线图：W5 之后（W6–W8）

> 依据 `目的与规划.md` §14.2 的 W5–W8 周块与 §14.3 的三个 Gate。
> 本文只列 **B 要做的事、依赖 A 的什么、怎么验收**；A 的任务见规范原文，不在此重复。
> 当前状态：W1–W5 已完成，见 `reports/week02_acceptance.md`、`week03_acceptance.md`、
> `week04_gate1.md`、`week05_acceptance.md`。最后更新：2026-11-10。

## 0. 已有基础（W5 收口时的状态）

| 已有 | 位置 | 说明 |
|---|---|---|
| 器件/状态/单位 | `core/`、`io/units.py`、`io/config.py` | W2 冻结；公共 API 由 `tests/test_public_api.py` 钉住 |
| 独立 oracle T01–T07 | `tests/test_oracle_*.py` | Bernoulli 用 60 位 Decimal，SG 用代数恒等式 |
| Gummel 求解 | `solvers/gummel.py`、`solvers/convergence.py` | 从 `C2_v1` 并入根包；双判据收敛 |
| Newton 求解 | `solvers/newton.py`、`solvers/jacobian.py` | W5 并入根包；九色稀疏差分 + 行缩放 + Armijo 线搜索 |
| T10 oracle | `verification/fd_jacobian_oracle.py` | 独立中心差分；实测 `1.65e-11` vs 阈值 `1e-4` |
| 运行编排 | `api.py`（`run`、`write_required_outputs`）、`simulation/{dc,sweep,diagnostics}.py` | 固定步进 + 自适应回滚、协作式取消 |
| 三文件 writer | `io/{results_csv,plot,report}.py` | 初版：22 列、6 面板、十节报告 |
| 基准与快照 | `examples/pn_baseline_0_to_030.yaml`、`baseline_0_to_030_snapshot.json` | 7 偏压全收敛，`config_sha256` 校验 |

测试基线：根包 `332 passed`，`C2_v1` 参考树 `215 passed`。

## 1. W5：Newton 的独立验证与自适应继续法（已完成 2026-11-10）

**依赖 A**：`solvers/newton.py`、`jacobian.py`、`convergence.py` 的 scaled 判据实现
（W5-A1…A6），以及 `SolveResult.iteration_history` 每步含 `iteration/alpha/residual_norm/update_norm`。
A 在 `e6c8a5f` 交付，B 已并入根包并核对数值逐位一致（`reports/week05_acceptance.md` §2）。

### W5-B1 独立 FD Jacobian oracle（`src/minitcad/verification/fd_jacobian_oracle.py`）

- 对**唯一残差** `evaluate_qf` 做中心差分；边界一侧越界时退化为单边差分并把该列标记为
  `ONE_SIDED`，不静默混用两种差分。
- 步长按变量尺度取：`δ = 1e-6·max(1, |x_j|)`；`n,p` 的试探必须保持为正，越界直接判该方向失败。
- 禁止 import A 的 `jacobian.py`；只允许 import `discretization.residual`。
- 验收：T10 的 `‖J_A d − J_FD d‖∞ / max(1, ‖J_FD d‖∞) ≤ 1e-4`，
  在 0 V 与 0.3 V 两个状态、随机方向 `d` 上各测一次。

```bash
python -m pytest tests/test_t10_jacobian.py -q
```

### W5-B2 自适应继续法与回滚（`simulation/sweep.py`）

- 步长状态机：成功→`step = min(step·grow, max_step)`；失败→`step = max(step·shrink, min_step)`
  并从**最后一个已接受解**重试同一个目标偏压；`step < min_step` 才记 `REJECTED`。
  （"重试同一个目标"的实际含义见 `simulation/sweep.py` 的模块文档：把失败的一跳拆小，
  总体目标不变。）
- 每次尝试都进 `attempts`（`PASS`/`FAIL`），回滚原因写进 `Attempt.note` 与 CSV 的
  `notes`；`successful` 只含真正收敛点。
- 验收（T13）：在 0.3 V 附近注入一次强制失败，断言 `ΔV` 收缩、回滚后仍能到达 0.30 V、
  且 `attempts` 里保留了失败记录。**已完成**，见 `reports/week05_acceptance.md` §4。

### W5-B3 冻结 `run()` 与取消语义

- 现有初版已实现 `progress(bias, status)` 与 `should_stop()`；W5 只做契约测试：
  取消后 `decision=REJECTED`、`failure_reason="CANCELLED"`、已接受点与三文件保留。
- 验收：`tests/test_cancel_contract.py`（7 项，含“取消后仍能写三文件”）。**已完成**。

## 2. W6：Gate 2（双求解器与自适应扫压）

**依赖 A**：Newton 已可用（W5 并入）；A 的 W6-A1 对照实验数据。
**已知风险**：Newton 在 `ΔV > 0.02 V` 时会在 0.10–0.15 V 之间失败，靠回滚兜住
（`reports/week05_acceptance.md` §4.1、§8 第 1 项）。W6 的 T11 只在**两个求解器都收敛**
的偏压上逐点对照，不把 Newton 的继续法缺陷算成“解不一致”。

| B 任务 | 产物 | 验收 |
|---|---|---|
| W6-B1 T10/T11/T13 正式通过 | `tests/test_t10_jacobian.py`、`test_t11_dual_solver.py`、`test_t13_rollback.py` | 势峰值差 `≤1e-3 V`；非零终端电流差 `≤1%`；回滚不传播失败态 |
| W6-B2 一键对比与回归 | `scripts/compare_solvers.py`、`scripts/run_gates.py`、故障注入工具 | 一条命令跑完 T08–T11/T13 并输出 PASS/FAIL 表 |
| W6-B3 细网格参考首版（T12） | `verification/reference.py` | 均匀网格 `M=201/401/801` 的终端电流外推，误差随 `h` 单调下降 |
| W6-B4 writer 补齐 | `io/{results_csv,plot,report}.py` | `diagnostic` 行（`current_nonconservation_nd`、`min_n_nd`、`boundary_max_error_nd`）；有参考数据时画 theory/reference 误差面板 |
| W6-B5 进度/取消契约 | `tests/test_cancel_contract.py`（W5 已提前完成） | `progress` 只报已尝试偏压；取消在两次尝试之间生效 |

Gate 2 失败处理：任一 T10/T11/T13 不通过即 Gate 2 = `REJECTED`，不删偏压点、不改阈值。

## 3. W7：核心冻结（`w07-core-freeze`）

| B 任务 | 产物 | 验收 |
|---|---|---|
| W7-B1 网格生成器 | `core/mesh.py` 增加结区 refine 规则（`refine_near(c, width, ratio)`） | cm/nd 同尺度；结区断点强制落在结点；非法参数显式失败 |
| W7-B2 细网格参考与收敛表（T12/T14/T15） | `verification/reference.py`、`reports/week07_mesh_convergence.md` | 误差随加密单调下降；表内每行可重跑 |
| W7-B3 三文件 schema 冻结 | `scripts/accept_outputs.py` | 逐条检查 §12.4：列名/`quantity`/`unit`、理论空值、`run_id`/SHA 一致、PNG 可读、报告十节 |
| W7-B4 CLI/API/YAML 冻结 | `docs/api_reference.md`、`docs/user_guide.md`、CI 配置 | `tests/test_public_api.py` 扩展为含 `run`/`write_required_outputs` 的签名冻结 |
| W7-B5 GUI API 冻结 | `docs/gui_api_contract.md` | `load_config/parse_config_dict/save_config/run/write_required_outputs` 与 `RunResult` 字段固化 |

冻结后 W8 只允许 GUI、联调与阻塞性 bug 修复。

## 4. W8：GUI 联调与 Windows 发布（`w08-release`）

| B 任务 | 产物 | 验收 |
|---|---|---|
| W8-B1 GUI/CLI 一致性 + T16/T17 | `tests/test_t16_gui_cli_parity.py`、`test_t17_cancel_gui.py` | 同一 YAML 的 `config_sha256`、偏压点、终端电流、三文件 schema 与 CLI 完全一致 |
| W8-B2 全回归与复现审计 | `reports/week08_gate3.md` | T01–T15 全绿；两次干净环境运行得到同一 `config_sha256` 与同一判决 |
| W8-B3 Windows 打包与 T18 | `scripts/build_windows.bat`、`docs/WINDOWS_BUILD.md` | onedir EXE 在干净 Windows 上加载 YAML、跑完基线、写出三文件；未实测则标 `NOT_VERIFIED` |
| W8-B4 发布收口 | release notes、tag `w08-release` | 最终 `ACCEPTED/REJECTED` 与演示记录 |

## 5. 与 A 的接口依赖（B 的阻塞点）

| B 需要 | 来自 A | 阻塞的 B 任务 | 若延期的降级方案 |
|---|---|---|---|
| `solve_newton` | W5-A1…A5 | T10/T11、Gate 2 | **已交付**（`e6c8a5f`，B 并入根包） |
| `SolveResult.iteration_history` 每步字段 | W5-A6 | 迭代行 CSV、收敛面板 | **已交付**：`iteration/alpha/residual_norm/update_norm/status`（另含 `merit_L2`/`linear_backward_error`，不进 CSV） |
| `solvers/convergence.py` 的 8.5 节尺度 | W5-A5 | B 的 oracle 对照 | **已交付**：`residual_scales_qf`/`scaled_update_physical`，根包与 `C2_v1` 同名 |
| Newton 继续法边界修复（接触投影后整流） | W6-A?（待 A 登记） | T11 的偏压覆盖范围 | 用自适应回滚兜住；T11 只在都收敛的偏压上对照 |
| 结区局部误差修复 | W7-A1 | T12 网格收敛表 | 先给误差数据，结论等 A 修完再写 |

## 6. 风险与顺序

1. **Newton 的继续法边界是当前最大风险**：`ΔV > 0.02 V` 时 Newton 在 0.10–0.15 V 会失败，
   自适应回滚能兜住（实测 0.03–0.07 V 各种起点都能到 0.30 V），但**代价是尝试次数翻倍**，
   且 0.6 V 的强偏压基准可能兜不住。已把根因（接触投影造成的虚假边界层通量）与两条
   可能修法写进 `reports/week05_acceptance.md` §8，等 A 确认。
2. **强偏压失败**：`examples/pn_left_p_right_n.yaml` 的 0.6 V 终点可能失败，
   回滚必须留下失败记录；Gate 2 只对**预先登记**的基准裁决。
3. **自适应路径不做快照**：尝试序列在失败边界附近对 `bias_V` 末位敏感
   （`week05_acceptance.md` §4.2），基线只记“是否到达 `stop_V`”与数值判据。
4. **非均匀网格**：`segmented` 已强制结区落结点；refine 规则在 W7-B1，不要提前改动
   已冻结的 `Mesh1D` 字段。
5. **文档与 schema 冻结**：W7 之后改 CSV 列名会让 W8 的 GUI 与历史报告失效，
   任何改动都要走 `docs/a_b_interface_week02.md` 的变更记录。W5 只新增了 `Attempt.note`，
   CSV 列序未变。
