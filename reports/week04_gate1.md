# W4 验收记录（B 侧）：Gate 1

> 对应 `目的与规划.md` §14.2（W4-B1…W4-B4）与 §14.3 的 Gate 1 退出条件。
> 日期：2026-11-04；Python 3.14.6；NumPy 2.5.3；SciPy 1.18.0；matplotlib 3.11.1（Agg）。
> 求解器：A 的 `solve_gummel`（QF 口径），本轮把 `solvers/` 并入根包并清理。
> 基准配置：`examples/pn_baseline_0_to_030.yaml`（同一器件，0 → 0.30 V，步长 0.05 V，201 结点）。

## 1. 复现命令与结果

```bash
cd <repo-root>
python -m pytest -q                       # 291 passed（0 skipped）
python -m pytest tests/test_week04_gate1.py tests/test_api_run.py \
    tests/test_io_writer.py tests/test_baseline_snapshot.py -q   # 29 passed
PYTHONPATH=src python -m minitcad.cli solve examples/pn_baseline_0_to_030.yaml \
    --out runs/pn_baseline_0_to_030 --overwrite
# decision: ACCEPTED；7/7 偏压收敛；三文件写入 runs/pn_baseline_0_to_030/
cd C2_v1 && python -m pytest -q           # 168 passed（A 的参考树，本轮已修好旧签名测试）
```

## 2. Gate 1 判据对照

Gate 1 要求：存在真实 0 V 解与至少一个正向非零偏压解；`|np-1| ≤ 1e-6`、平衡零电流
`≤ 1e-6`、接触中性误差 `≤ 1e-10`、总电流守恒达标；T08/T09 通过；记录可重跑。

| 指标 | 实测 | 阈值 | 结论 | 证据 |
|---|---|---|---|---|
| 0 V 质量作用 `max|np-1|` | `2.220e-16` | `1e-6` | PASS | `convergence.md` §7 |
| 0 V 终端电流密度 `|J|` | `0.0`（精确） | `1e-6 A/cm²` | PASS | 同上 |
| 接触中性（两端 `u,n,p` 与 `contact_targets`） | `0.0` | `1e-10` | PASS | 同上 |
| 总电流沿空间起伏（相对 `max|J|`） | `2.970e-09` | `1e-3` | PASS | 同上 |
| 最大 scaled residual | `6.694e-10` | `1e-8` | PASS | 同上 |
| 最大 scaled update | `2.144e-07` | `1e-6` | PASS | 同上 |
| 正向偏压解 | 0.05…0.30 V 共 6 点，全部收敛 | ≥ 1 点 | PASS | 偏压表 §6 |
| T08（0 V 完整 PN） | 见 `tests/test_week04_gate1.py` | — | PASS | 3 项断言 |
| T09（总电流守恒） | 每偏压独立断言 | `1e-3` | PASS | 同上 |

正向电流随偏压的变化（`I_A = A·J`，`A = 1e-4 cm²`）：

| V_app [V] | 0.00 | 0.05 | 0.10 | 0.15 | 0.20 | 0.25 | 0.30 |
|---|---|---|---|---|---|---|---|
| `J_terminal` [A/cm²] | 0 | 8.585e-9 | 3.990e-8 | 1.906e-7 | 1.055e-6 | 6.503e-6 | 4.257e-5 |
| 迭代数 | 0 | 3 | 3 | 3 | 3 | 3 | 4 |

相邻点算出的理想因子 `n = ΔV / (V_T ln(J_k/J_{k-1}))` 落在 `1.00–1.06`
（`tests/test_week04_gate1.py::test_forward_current_is_exponential_with_ideality_near_one`
只断言 `1.0 ≤ n ≤ 1.5`，具体值随偏压小幅变化，属 SRH 复合区的正常偏离）。

## 3. 三文件（§12.1–§12.4）

一次 `minitcad solve` 写出固定三文件，`run_id` 与 `config_sha256` 在三者中一致：

| 文件 | 大小 | 关键校验 |
|---|---|---|
| `comparison.csv` | 21 185 行 | 22 列表头与 §12.1 逐字一致；`sweep` 行覆盖 7 个尝试偏压；`node` 行 `11×201×7`；`edge` 行 `4×200×7`；`benchmark` 行 9 条 |
| `comparison.png` | 6 面板 | I–V、φ(x)、n/p（对数轴）、J 与 E、收敛历史、验收摘要表；图头含方向/极性 |
| `convergence.md` | 十节 | §1–§10 齐全；第 7 节指标与上表同源 |

理论/误差列遵守 §12.1：逐点 `phi/n/p` 无闭式解 → `theory_kind=none`、理论列留空；
0 V 零电流写绝对误差、相对误差留空并在 `notes` 说明；接触浓度、接触势、内建电势、
质量作用为精确值 → `theory_kind=exact`、误差 0。

`minitcad validate <csv>` 同时接受 cm 与 SI 两种表头（仅第 9 列 `position_cm`/`position_m` 不同）。

## 4. 本轮对 A 代码的改动

| 文件 | 改动 | 验证 |
|---|---|---|
| `C2_v1/tests/test_qf_stability.py` | 4 个失败用例改用 `QFState` + `evaluate_qf`（旧 `residual_qf(device,bias,u,eta_n,eta_p,ref_n=…)` 已不存在） | `C2_v1` 168 passed |
| `C2_v1/src/minitcad/solvers/gummel.py` | 删除 148 行被三引号注释的旧版 `solve_gummel`；`convergence.py` 的重复逻辑抽到 `scaled_residual_qf`/`scaled_update`；补模块与函数文档；裸 `except:` 改为明确的四类异常；trial 点的 `exp` 溢出用 `errstate` 屏蔽（仍由 `State`/`Fluxes` 有限性检查与线搜索拒绝） | 重写前后 5 个偏压的 `iterations/residual_norm/update_norm/J/u_mid/n_sum` **逐位相同** |
| `C2_v1/src/minitcad/solvers/convergence.py` | 原文件被截断（函数体停在 `m = qf.node_count`），补全为可测试的收敛判据模块 | `C2_v1` 168 passed |
| 根包 `src/minitcad/solvers/` | 从 `C2_v1` 并入 `gummel.py`、`convergence.py`（单一包原则，§10.1） | 根包 291 passed |

## 5. B 的 W4 交付

| 任务 | 产物 | 证据 |
|---|---|---|
| W4-B1 T08/T09 | `tests/test_week04_gate1.py`（7 项） | 上表 |
| W4-B2 基准与快照 | `examples/pn_baseline_0_to_030.yaml`、`examples/baseline_0_to_030_snapshot.json`、`tests/test_baseline_snapshot.py` | 快照含 `config_sha256` 与逐偏压电流；配置被改时测试先失败 |
| W4-B3 固定步进扫描 | `src/minitcad/simulation/sweep.py`（`run_fixed_sweep`、`equilibrium_start_state`） | `tests/test_api_run.py` |
| W4-B4 `run()` + CLI + writer | `api.py`（`run`、`write_required_outputs`）、`cli.py`（`solve` 接线）、`io/results_csv.py`（行构造）、`io/plot.py`、`io/report.py`、`simulation/dc.py`、`simulation/diagnostics.py` | `tests/test_api_run.py`、`tests/test_io_writer.py`、`tests/test_cli.py` |
| 求解器并入 | `src/minitcad/solvers/{__init__,gummel,convergence}.py` | `python -m pytest -q` |

## 6. 遗留与下一周阻塞项

1. **Newton 未落地**：`solver.method: newton` 在 `run()` 中抛 `NotImplementedError`，CLI 退出码 3。
   W5-A1…A6 交付后才能做 T10/T11；B 的 FD Jacobian oracle 可以先按 `evaluate_qf` 写好主体。
2. **扫描只有固定步进**：失败即停止，没有回滚与步长收缩（W5-B2）。T13 需要它们。
3. **writer 是初版**：有理论/参考的误差面板、`diagnostic` 行（`current_nonconservation_nd`、
   `min_n_nd`、`boundary_max_error_nd`）留到 W6-B4 补齐；PNG 未画理论/参考曲线（当前没有合法参考数据）。
4. **`area_cm2` 可空**仍未合并（B 的契约要求显式给出；要改需 A 先给 writer 缺省规则）。
5. **`C2_v1` 与根包的 `normalization.py` 分叉**仍在；根包用 B 的 W1 版，A 的 `e0_v_cm` 命名未采纳。

## 7. 结论

- Gate 1 的数值判据全部达标，且由 `convergence.md` 第 7 节与
  `tests/test_week04_gate1.py` 双向记录。
- 三个固定文件可由 `minitcad solve` 一次产出，`minitcad validate` 能独立复核表头。
- 判定：**Gate 1 = PASS**（B 侧）；A 侧需确认 `solve_gummel` 的 `ref_n/ref_p` 锚定口径
  （代码用 `ref_n` 取右端、`ref_p` 取左端，残差对该选择不变，见
  `tests/test_oracle_qf_residual.py::test_internal_rows_are_invariant_to_the_qf_reference_shift`）。
