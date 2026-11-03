# convergence.md 报告模板（第 12.3 节固定结构）

> **Owner**：B（W1-B4 模板，W4 起由 `run()` 自动填写）。复制本文件为 `<output.directory>/convergence.md`，
> 按十个固定章节填写；**未实际运行的通过数量、日志、截图不得写入**（§14.5）。
> 与 `comparison.csv`、`comparison.png` 共同构成一次 `solve` 的三文件契约；三文件的 `run_id` 必须一致（§12.4）。

```text
# convergence.md — <run_id>

## 1. Run summary
run_id            : <run_id>
git sha           : <git rev-parse --short HEAD>
config_sha256     : <sha256 of normalized YAML>
date / platform   : <YYYY-MM-DD> / <Windows 11 / Linux>
outcome           : ACCEPTED | REJECTED（失败点见第 6、10 节）

## 2. Orientation
P(left, x=0) -> N(right, x=L) | V_app = V_left - V_right | +V = forward bias | right = 0 V

## 3. Physics / units
input_unit_system  : cm | si
core_unit_system   : dimensionless
output_unit_system : cm | si
scales (cm)        : eps_s[F/cm]  ni[cm^-3]  mu[cm^2/(V s)]  L0[cm]  D0[cm^2/s]  J0[A/cm^2]  R0[cm^-3/s]
scales (si)        : eps_s[F/m]   ni[m^-3]   mu[m^2/(V s)]   L0[m]   D0[m^2/s]   J0[A/m^2]   R0[m^-3/s]
invariant          : tau[s], V_T[V]
material           : <name> @ <T> K, parameter_policy = explicit_at_T | named_model
T = <K>, V_T = <V>, n_i = <cm^-3>, eps_r = <->, tau_n/tau_p = <s>/<s>
SI audit           : <SI 输入经适配层换算成 cm 的实际因子与舍入记录；SI 输出仅发生在 writer>

## 4. Solver settings
method = newton | gummel ; jacobian = analytic | sparse_fd ; max_iter = <n>
alpha_min = <x> ; positive_scale_floor = <x> ; continuation = none | <策略>

## 5. Tolerances
tol_residual = <x> (scaled) ; tol_update = <x> (scaled)
current conservation tol = <x> ; contact neutrality tol = 1e-10 ; np-1 tol = <x>

## 6. Bias table
| V_app[V] | solver | step[V] | success | iterations | final Rnorm | final Dnorm | runtime[s] | fail reason |
|---|---|---|---|---|---|---|---|---|
| 0.000    | gummel | -0.05   | PASS    | <n>        | <x>         | <x>         | <t>        | -           |
| ...      |        |         | FAIL    | <n>        | <x>         | <x>         | <t>        | MAX_ITER / LINE_SEARCH / ... |

## 7. Validation
| 项 | 实测 | 阈值 | 结论 |
|---|---|---|---|
| Bernoulli（T01） | <max rel err> | 3e-14 | PASS/FAIL |
| SRH 异寿命（T02） | <err vs 19/47> | 1e-12 | PASS/FAIL |
| 接触/极性（T03/T04） | <err> | 1e-10 | PASS/FAIL |
| 非均匀 Poisson（T05） | <residual> | 1e-11 | PASS/FAIL |
| SG 平衡零电流（T06） | <max |j|> | <x> | PASS/FAIL |
| SG 纯扩散（T07） | <err> | 1e-12 | PASS/FAIL |
| 电流守恒（T09） | <spread> | <x> | PASS/FAIL/N_A |
| Jacobian 差分（T10） | <rel err> | <x> | PASS/FAIL/N_A |
| Mesh / 状态布局（I-01） | <ok> | - | PASS |

## 8. Error policy
theory_kind 适用域：exact = 接触与 0 V 平衡量；approximation = 耗尽近似等；none = 有偏压的逐点 phi/n/p。
Q_floor（cm / si）：<各量阈值>；零理论值不做相对误差除法（CSV `theory_rel_error` 留空并在 notes 说明）。

## 9. Reproducibility
python = <ver> ; numpy = <ver> ; scipy = <ver> ; pytest = <ver>
command : python -m minitcad.cli solve <yaml> --out <dir>
config snapshot:
```yaml
<完整 YAML 原文>
```

## 10. Final decision
decision : ACCEPTED | REJECTED
未通过项 : <列表，含 Owner 与阻塞的下一周任务>（无则写 none）
```

## 填写规则

1. **只写实测值**：任何 `PASS` 必须能由 `comparison.csv` 的对应行或一条测试命令复现。
2. **单位一致**：`output.unit_system: si` 时第 3、7、8 节的物理量、基准尺度、绝对误差与 `Q_floor` 一并换成 SI，
   同时保留 cm 内核尺度的审计记录；不得只改文字不改数值。
3. **失败不隐藏**：每个尝试偏压（含失败的）都要在第 6 节出现，并写 `fail reason`；
   未到终点不得假报成功（§11.3）。
4. **不是第二份 CSV**：第 6 节只放每偏压汇总，逐次迭代由 `comparison.csv` 的 `iteration` 行承担。
5. **双求解器**：Newton 与 Gummel 都跑时分别给出各自的迭代与容差结果，不用平均值掩盖失败。
