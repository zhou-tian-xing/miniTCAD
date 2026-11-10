# W5 验收记录（B 侧）：Newton 接入、T10、T13 与取消契约

> 对应 `目的与规划.md` §14.2（W5-B1…W5-B3）与 `docs/B_roadmap.md` §1。
> 日期：2026-11-10；Python 3.14.6；NumPy 2.5.3；SciPy 1.18.0；matplotlib 3.11.1（Agg）。
> 本轮同步 A 的 `e6c8a5f`（"add Newton solver and reverse bias simulation"），
> 把 `solvers/{jacobian,newton}.py` 并入根包并做**逐位一致性**核对。
> Gate 2 本身在 W6，本文件只给 T10/T13 与取消契约的证据。

## 1. 复现命令与结果

```bash
cd <repo-root>
python -m pytest -q                                   # 332 passed（W4 收口时 296）
python -m pytest tests/test_t10_jacobian.py tests/test_t13_rollback.py \
    tests/test_cancel_contract.py tests/test_newton_adapter.py -q    # 33 passed
PYTHONPATH=src python -m minitcad.cli solve examples/pn_newton_adaptive_0_to_030.yaml \
    --out runs/pn_newton_adaptive_0_to_030 --overwrite
# decision: ACCEPTED；到达 0.30 V（中途 3 次回滚，路径见 §4.2）；三文件写入 runs/pn_newton_adaptive_0_to_030/
cd C2_v1 && python -m pytest -q                       # 215 passed（A 的参考树）
```

## 2. 并入根包的 Newton 与 A 的实现逐位一致

`C2_v1/src/minitcad/solvers/newton.py` 与 `jacobian.py` 被并入 `src/minitcad/solvers/`，
补了模块/函数文档、把裸 `except` 换成显式异常类别、增加了与 `solve_gummel` 对齐的
`diagnostics` 参数。**数值行为未改**：对同一器件、同一初值链，在 0、0.05、0.10、0.15、
0.20、0.25、0.30 V 七个偏压上，两边的 `converged / iterations / residual_norm /
update_norm / reason` 完全相同，且最终 `u/n/p` 三个数组的 SHA-256 逐位相同（含 0.10、
0.15、0.20 三个失败点）。

这条结论可以用仓库里的脚本复核（两个包同名，必须分两次进程跑）：

```bash
python scripts/compare_newton_port.py run --src src        --out .tmp/root.json
python scripts/compare_newton_port.py run --src C2_v1/src  --out .tmp/c2v1.json
python scripts/compare_newton_port.py diff .tmp/root.json .tmp/c2v1.json   # → identical
```

`solvers/convergence.py` 有两处结构改动，都是为了让两个求解器共用同一套判据：

| 改动 | 原因 |
|---|---|
| 抽出 `residual_scales_qf()` | Newton 需要**逐行**尺度做行缩放，Gummel 只需要它的 `max`；两者必须同源，否则“解得很准却判成没收敛” |
| `scaled_update` → `scaled_update_physical(old, new, *, n_floor)` | 判据作用在物理量上，名字必须说明这一点；旧名保留为薄封装，W4 的调用点与测试不受影响 |

## 3. T10：独立 FD Jacobian oracle

`src/minitcad/verification/fd_jacobian_oracle.py` 对 `discretization.residual.residual_qf`
做中心差分（`h_j = 1e-6·max(1,|z_j|)`），逐列记录差分方式；**不 import** `solvers.jacobian`
（`test_oracle_only_depends_on_the_residual` 直接读源码检查 import 行）。

判据（§14.3）：`‖J_A d − J_FD d‖∞ / max(1,‖J_FD d‖∞) ≤ 1e-4`。

| 状态 | 实测比值 | 单边列数 | 阈值 | 结论 |
|---|---|---|---|---|
| 0 V（平衡解） | `1.650e-11` | 0 | `1e-4` | PASS |
| 0.30 V（正向解） | `1.650e-11` | 0 | `1e-4` | PASS |

补充：`0.30 V` 上另取 5 个独立随机方向，最大比值仍在 `1.650e-11`；抽查 17 列的
**逐列**相对误差最大为 `5.3e-07`（0 V）与 `3.6e-10`（0.30 V），也都低于 `1e-4`。
逐列误差比方向导数误差大两三个量级是正常的：方向导数被量级 `1e12` 的 Poisson 行主导，
而逐列比较的是量级只有 `1e-1` 的那些列。

误差来源已定位到**空穴连续性行**：`argmax` 的比值在所有随机方向上恒为 `1.650e-11`，
即该行的相对误差与方向无关，属差分的截断/舍入，不是实现缺陷。

## 4. T13：自适应继续法与回滚

`simulation/sweep.py` 新增 `run_adaptive_sweep()`：成功则 `step ← min(step·1.5, max_step)`，
失败则记 `FAIL` 并 `step ← max(step·0.5, min_step)`，从**最后一个已接受解**重新逼近；
缩到 `min_step` 仍失败就停止（`run()` 据此判 `REJECTED`）。

### 4.1 为什么必须有它

固定 0.05 V 步长下 Newton 走不到 0.30 V：

```text
$ PYTHONPATH=src python -m minitcad.cli solve <pn_baseline 改成 newton>
FAIL at 0.1 V: LINE_SEARCH_FAILED          → decision: REJECTED
  0.00  PASS  0 次迭代
  0.05  PASS  6 次迭代
  0.10  FAIL  3 次迭代（residual 1.28，alpha 折到 1e-3 仍不下降）
```

原因是接触投影把左端 `η_n` 一次抬 `Δv = 1.93`，内部结点仍是旧值，第一条边被算成
`1e12` 量级的虚假空穴通量；该行尺度随之变成 `1e12`，merit 失去分辨力。实测的可用阶跃
上限在 `ΔV ≈ 0.02 V`：

| 继续步长 | 0 → 0.30 V 结果 |
|---|---|
| 0.05 V | 在 0.10 V 失败 |
| 0.02 V | 全部收敛 |
| 0.01 V | 全部收敛 |

### 4.2 自适应后的实测路径

`examples/pn_newton_adaptive_0_to_030.yaml`（起点 0.05 V，`grow 1.5`、`shrink 0.5`、
`min_step 1e-3`）：

```text
0.0000 + | 0.0500 + | 0.1250 - | 0.0875 - | 0.0688 + | 0.0969 - |
0.0828 + | 0.1039 + | 0.1355 - | 0.1197 + | 0.1435 + | 0.1791 + |
0.2324 + | 0.3000 +                                   → ACCEPTED，14 次尝试，4 次回滚
```

四次失败都写在 `attempts` 里，并带 `rollback: step … -> …` 说明；CSV 的 `notes` 列同时
保留扫描层说明与求解器原因：

```text
0.125     FAIL  notes=rollback: step 0.075 -> 0.0375; LINE_SEARCH_FAILED
0.0875    FAIL  notes=rollback: step 0.0375 -> 0.01875; LINE_SEARCH_FAILED
0.096875  FAIL  notes=rollback: step 0.028125 -> 0.0140625; LINE_SEARCH_FAILED
0.135547  FAIL  notes=rollback: step 0.0316406 -> 0.0158203; LINE_SEARCH_FAILED
```

**中途偏压序列不可复现，但判决可复现**：起点步长取 0.03/0.04/0.05/0.06/0.07 时，五组
都得到 `ACCEPTED` 并到达 0.30 V（回滚 2–4 次），但尝试次数在 10–15 之间变化——Newton
在 0.10–0.15 V 这一段处在失败边界附近，`bias_V` 末位差 `1e-16` 就能改变一次尝试的成败
（本轮就因为把自适应偏压取到 12 位小数，使 0.1355 V 由通过变成失败）。
因此 W5 **不对自适应路径做快照测试**，只对“到达 `stop_V` / 回滚被记录”做断言
（`tests/test_t13_rollback.py`）。

### 4.3 T13 的四条断言

| 断言 | 证据 |
|---|---|
| 注入一次强制失败后步长收缩 | `test_rollback_shrinks_the_step_and_still_reaches_the_target` |
| 回滚后仍到达 0.30 V | 同上（`reached_V == 0.3`） |
| 失败记录保留在 `attempts` 且带 `rollback:` 说明 | 同上 + 上面的 CSV 片段 |
| 重试从最后一个已接受解出发（状态不回退） | `test_retry_starts_from_the_last_accepted_state` |
| 缩到 `min_step` 就放弃，不无限重试 | `test_exhausted_step_gives_up_instead_of_retrying_forever` |
| 固定步进语义未被改写（失败即停） | `test_fixed_sweep_still_stops_at_the_first_failure` |

## 5. W5-B3：进度与取消契约

| 契约 | 证据 |
|---|---|
| 取消 → `REJECTED` + `failure_reason="CANCELLED"` | `test_cancel_keeps_accepted_points_and_is_not_accepted` |
| 已接受的偏压点完整保留，且仍能写出三文件 | `test_cancel_still_writes_the_three_files`（读回 CSV 校验 `run_id` 一致） |
| `progress` 只报已尝试偏压，序列与 `attempts` 逐个对齐 | `test_progress_never_reports_unattempted_biases` |
| 自适应模式下取消语义一致 | `test_cancel_in_adaptive_mode` |
| `should_stop` 每次尝试只查一次（求解过程内无取消点） | `test_should_stop_is_not_consulted_inside_the_solver` |

取消是协作式的：只在两次尝试之间生效。理由写在 `tests/test_cancel_contract.py` 的模块
文档里——中途打断 Newton 会留下一个既没收敛也没记录的状态。

## 6. 双求解器一致性（T11 的前置证据，正式对照在 W6）

同一 0.25 V 解热启动、同一偏压 0.30 V：

| 指标 | 实测 | W6 判据 |
|---|---|---|
| `max|Δu|`（无量纲） | `3.351e-11`（合 `8.663e-13 V`） | 势峰值差 `≤1e-3 V` |
| `max|Δn|/n`、`max|Δp|/p` | `3.351e-11`、`3.441e-11` | — |
| 终端电流密度 | `1.086222e+01` vs `1.086222e+01`，相对差 `7.527e-12` | 非零电流相对差 `≤1%` |

## 7. 本轮对 A 代码的改动

| 文件 | 改动 | 验证 |
|---|---|---|
| 根包 `solvers/jacobian.py`、`newton.py` | 从 `C2_v1` 并入（新增）；补文档、显式异常、`diagnostics` 参数 | 七个偏压逐位一致（§2） |
| 根包 `solvers/convergence.py` | 抽出 `residual_scales_qf`；`scaled_update` → `scaled_update_physical`（旧名保留） | 根包 332 passed |
| 根包 `simulation/dc.py` | 按 `solver.method` 分派，删掉 `NotImplementedError` | `tests/test_cli.py`、`test_newton_adapter.py` |
| 根包 `simulation/sweep.py` | 新增 `run_adaptive_sweep`；`SweepOutcome` 增加 `successful`/`reached_V` | `tests/test_t13_rollback.py` |
| 根包 `simulation/run_result.py` | `Attempt` 增加 `note`（扫描层说明，与 `result.reason` 分工） | `tests/test_public_api.py` |
| 根包 `core/device.py` | `SweepConfig` 增加 `adaptive/min_step_V/max_step_V/grow_factor/shrink_factor` 并校验 | `tests/test_io_config.py` |
| 根包 `io/config.py` | `sweep` 块真正解析自适应参数（之前接受但丢弃） | 同上 |
| 根包 `io/results_csv.py` | `notes` 合并回滚说明与求解器原因；新增 `read_rows` 供验收脚本使用 | `tests/test_cancel_contract.py` |
| `scripts/compare_newton_port.py`（新增） | 分进程核对两个包的 Newton 逐位一致（`run` / `diff` 两个子命令） | 上表 §2 的 `identical` |
| `C2_v1` | **未改**：A 的 `e6c8a5f` 树自洽（已有 `residual_scales_qf`/`scaled_update_physical`） | `C2_v1` 215 passed |

## 8. 遗留与下一周

1. **Newton 在 0.10–0.15 V 的鲁棒性边界太窄**：`ΔV > 0.02 V` 就会失败，靠回滚兜住。
   这是继续法能解决的，但**根因**在 §4.1 的虚假边界层通量——建议 A 评估两点：
   （a）接触投影后先对内部结点做一次线性外插（把 `η` 从旧值整流到新边界附近）；
   （b）线搜索的 merit 改用**求解开始时冻结**的尺度，避免被瞬时大通量污染。
   B 侧已测过 (b) 单独使用无效（`0.10 V` 仍失败），所以 (a) 更可能是关键。
2. **`run()` 尚无故障注入接口**：T13 目前通过 `run_adaptive_sweep(solve=...)` 注入。
   W6-B2 的 `scripts/run_gates.py` 需要把它提升成正式参数或独立脚本。
3. **`Attempt.note` 是 W5 新增字段**：CSV schema 未变（仍写 `notes` 列），但
   `RunResult` 的 dataclass 字段表变了，已同步 `tests/test_public_api.py` 与
   `docs/a_b_interface_week02.md`。
4. **自适应路径不做快照**（§4.2）：W6 的 Gate 2 判决应基于“是否到达终点 + 数值判据”，
   不要把尝试序列写进基线文件。

## 9. 结论

- T10 达标（`1.65e-11` vs 阈值 `1e-4`，两个状态、6 个方向、0 个单边列）。
- T13 达标（注入失败 → 收缩 → 到达 0.30 V，失败记录保留）。
- W5-B3 取消契约达标，取消后仍能写出三文件。
- 两个求解器在 0.30 V 上一致到 `7.5e-12`，远优于 Gate 2 的 `1%`。
- 判定：**W5 = PASS**（B 侧）。Gate 2 仍待 W6 的 T10/T11/T13 一并裁决；
  Newton 的继续法边界（§8 第 1 项）需 A 确认是否在 W6 修。
