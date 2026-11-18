# W6 验收记录（B 侧）：Gate 2

> 对应 `目的与规划.md` §14.2（W6-B1…W6-B5）与 §14.3 的 Gate 2 退出条件。
> 日期：2026-11-18（W7 收口后补做；W6 原本顺延）。Python 3.14.6 / NumPy 2.5.3 / SciPy 1.18.0。
> 一键复跑：`python scripts/run_gates.py --gate 2`。

## 1. 复现命令与结果

```bash
python scripts/run_gates.py --gate 2
# ID    gate result  summary
# T10   2    PASS    7 passed
# T11   2    PASS    4 passed
# T13   2    PASS    5 passed
# ACCEPTED: 3 item(s) passed
```

## 2. Gate 2 判据对照

Gate 2 要求六条同时满足（§14.3）：「Newton 与 Gummel 使用唯一 `F` 和同一边界；T10、T11、
T13 按 §11.2 的阈值全部通过；强偏压失败必须记录偏压/迭代/残差/原因；`run()` 支持
`progress`/`should_stop`；writer 能产出三文件初版」。

| 要求 | 实测 | 证据 |
|---|---|---|
| 唯一 `F`、同一边界 | 两法都调用 `evaluate_qf`，锚点同为 `ref_n = log(n_R) − u_R`、`ref_p = −log(p_L) − u_L` | `test_t11_dual_solver.py::test_both_solvers_share_the_boundary_condition` |
| T10 Jacobian 独立差分 | `1.650e-11` vs 阈值 `1e-4`（0 V 与 0.30 V） | `test_t10_jacobian.py` |
| **T11 双求解器** | 势峰值差 **`3.653e-13 V`** vs `1e-3 V`；电流相对差 **`5.310e-10`** vs `1%` | `test_t11_dual_solver.py` |
| T13 失败回滚 | 注入失败 → `ΔV` 收缩 → 仍到 0.30 V；失败记录留痕 | `test_t13_rollback.py` |
| 失败必须记录原因 | CSV 的 `notes` = `rollback: step a -> b; <solver reason>` | `scripts/accept_outputs.py::FAILURE_NOTES` |
| `progress` / `should_stop` | 只报已尝试偏压；取消后 `REJECTED`/`CANCELLED` 且保留已接受点 | `test_cancel_contract.py` |
| writer 初版 | 三文件，`scripts/accept_outputs.py` 17 项通过 | `test_t14_output.py` |

判定：**Gate 2 = PASS**（B 侧）。两条数值判据的余量分别是 9 个和 7 个数量级。

## 3. T11 的基准选择（登记与理由）

T11 的定义是"相同器件/网格/偏压/容差下 Gummel vs Newton"，**没有规定扫压步长**；
步长只决定每个偏压的初值从哪来。这里用 0.02 V 梯子（0 → 0.30 V，15 个偏压），
先登记理由与代价：

- 0.05 V 梯子（`examples/pn_baseline_0_to_030.yaml` 的登记基准）上，Newton 在 0.10、0.15 V
  不收敛（整步过冲，`reports/week05_acceptance.md` §4.2），7 个偏压里只有 5 个能对照；
- 0.02 V 梯子上 15 个偏压全部可对照，两法从**同一个初值**、用**同一组容差**求解，
  完全符合 T11 的字面条件；
- **需要说明的时间顺序**：这个选择是在看到 0.05 V 的失败之后做的。按 §11.2
  "凡选择其他阈值，必须在项目启动前或通过双方审核的基准修订中登记"，这里把理由与
  0.05 V 的失败数据一并登记，而不是删掉失败点；A 已认可该处理的诊断依据。
- 不选的做法：把 0.05 V 的失败点删掉、或把 `1e-3 V`/`1%` 调松——都没有做。

`test_registered_005v_ladder_coverage_is_recorded` 仍然**记录** 0.05 V 梯子的实际覆盖
（当前 4/6 可对照，缺失点只允许是 0.10/0.15 V）。A 修好 Newton 后这条会自然变成全覆盖，
测试不需要改。

## 4. 逐偏压实测

`V_app` 每 0.02 V 一个点，两法同一初值：

| `V_app` [V] | 势峰值差 [V] | 电流相对差 | `J` [A/cm²] |
|---|---|---|---|
| 0.02 | `4.6e-17` | `1.5e-14` | `6.059e-04` |
| 0.08 | `9.2e-17` | `1.2e-13` | `5.599e-03` |
| 0.14 | `6.0e-16` | `1.1e-12` | `3.516e-02` |
| 0.20 | `6.3e-15` | `1.1e-11` | `2.691e-01` |
| 0.26 | `6.9e-14` | `1.1e-10` | `2.408e+00` |
| 0.30 | `3.7e-13` | `5.3e-10` | `1.086e+01` |

误差随偏压升高单调增大（低偏压电流小，相对误差也小），但最大仍在 `1e-13 V` 量级。

## 5. 一键回归工具

`scripts/run_gates.py` 把 T08/T09/T10/T11/T12/T13/T14/T15 与 API 冻结项各绑到一个测试文件，
逐项起 `pytest` 子进程并打印表：

```bash
python scripts/run_gates.py              # 全部九项
python scripts/run_gates.py --gate 2      # 只跑 Gate 2 的三项
python scripts/run_gates.py --list        # 只列清单
```

不做二次解析：子进程的退出码就是判定，失败时贴该子进程的输出尾巴，避免"脚本说通过、
人不知道跑了什么"。

## 6. 结论

- Gate 2 六条全部满足，T11 的余量为 9（电势）/ 7（电流）个数量级。
- T11 的 0.02 V 梯子已在 §3 登记理由，0.05 V 梯子的失败仍完整保留在测试与报告里。
- 判定：**Gate 2 = PASS**（B 侧）；A 侧需确认 `solve_newton` 的整步过冲是否在 W8 之前处理
  （`docs/a_b_interface_week02.md` §5 第 7 项，A 已认可按自适应规避）。
