# 用户指南：安装、配置与运行

> 面向"拿到仓库就要跑出一个结果"的读者。接口细节见 `docs/api_reference.md`，
> 数值约定见 `docs/numerical_conventions.md`。

## 1. 环境与安装

Python ≥ 3.10；依赖 `numpy`、`scipy`、`PyYAML`；画图需要 `matplotlib`。

```bash
pip install -e ".[dev]"      # 核心 + 测试 + matplotlib
```

不装包也可以用，把源码目录加进路径：

```bash
export PYTHONPATH=src         # Windows PowerShell: $env:PYTHONPATH="src"
```

## 2. 两条命令

```bash
minitcad validate examples/pn_baseline_0_to_030.yaml   # 只校验配置，不求解
minitcad validate runs/out/comparison.csv              # 只校验 CSV 表头
minitcad solve   examples/pn_baseline_0_to_030.yaml --out runs/out
```

`validate` 打印启动摘要（方向、单位、结点数、`V_T`、尺度、求解器、扫压范围）；
`solve` 打印每个偏压的 `PASS`/`FAIL`、判决、`run_id`，并写出三个固定文件。

| 退出码 | 含义 |
|---|---|
| 0 | `ACCEPTED` |
| 2 | `REJECTED`（求解失败或配置非法） |
| 3 | 保留（求解器未实现） |

## 3. 三个固定输出文件

| 文件 | 内容 |
|---|---|
| `comparison.csv` | 22 列 long-form：`sweep` / `node` / `edge` / `iteration` / `benchmark` 行 |
| `comparison.png` | 6 面板：I–V、φ(x)、n/p（对数轴）、J 与 E、收敛历史、验收摘要 |
| `convergence.md` | 十节：Run summary / Orientation / Physics-units / Solver / Tolerances / Bias table / Validation / Error policy / Reproducibility / Final decision |

独立的验收器（§12.4）：

```bash
python scripts/accept_outputs.py runs/out              # 逐条打印检查项
python scripts/accept_outputs.py runs/out --quiet      # 只打印失败项
```

## 4. 配置（YAML）

方向固定为 **左 P、右 N，右端 0 V**，`V_app = V_left - V_right`，正偏压是正向。
长度类字段以 `_cm` 或 `_m` 结尾，浓度类以 `_cm3` 或 `_m3` 结尾，**后缀必须和
`units.system` 一致**；写错或漏写直接报错。

| 块 | 关键字段 |
|---|---|
| `units` | `system`: `cm`（默认）或 `si` |
| `device` | `type: pn`、`orientation: p_left_n_right`、`length_cm`、`area_cm2`、`temperature_K` |
| `material` | `name: silicon`、`epsilon_relative`、`intrinsic_density_cm3`、`mobility_*`、`lifetime_*`、`bandgap_eV` |
| `regions` | 每段 `start_cm`/`end_cm`/`acceptor_cm3`/`donor_cm3`；结区在相邻两段的交界 |
| `mesh` | 见下 |
| `normalization` | `reference_diffusivity_cm2_s` 或 `reference_mobility_cm2_per_Vs`（二选一） |
| `contacts` | `right_voltage_V` 必须是 `0.0` |
| `solver` | `method: gummel \| newton`、`tol_residual`、`tol_update`、`max_iter`、`alpha_min`、`positive_scale_floor` |
| `sweep` | `start_V`、`stop_V`、`step_V`；自适应再给 `adaptive: true`、`min_step_V`、`max_step_V`、`grow_factor`、`shrink_factor` |
| `output` | `directory`、`unit_system: cm \| si`（只影响 writer） |

### 网格三种

```yaml
mesh: {type: uniform, nodes: 201}                       # 等距
mesh: {type: segmented, breakpoints: [0, 2.5e-4, 1e-3], per_segment_nodes: [11, 191]}
mesh: {type: refined, nodes: 201, refine: 20}            # 结点数不变，结区最密
```

`refined` 的结区坐标取自 `regions`，不需要重复声明；`refine` 是"远端单元 / 结区单元"
的长度比（`1` 即均匀）。结区加密对这份器件把 201 结点的电流误差从 `1.6e-3` 压到
`1.1e-5`，见 `reports/week07_mesh_convergence.md`。

## 5. 选求解器

| | Gummel | Newton |
|---|---|---|
| 何时用 | 默认，低偏压稳 | 耦合强（高正向偏压）时更快 |
| 步长容忍 | 任意 | 需要较小的继续步长，建议配 `adaptive: true` |

固定 0.05 V 步长下 Newton 会在 0.10 V 失败；`examples/pn_newton_adaptive_0_to_030.yaml`
用自适应步长跑通到 0.30 V，失败与回滚都记录在 `attempts` 和 CSV 的 `notes` 里。

## 6. 复现与追溯

- 每个文件的每一行都带 `run_id` 与 `config_sha256`；报告另有平台/Python/求解器设置。
- 同一配置重跑得到逐位相同的物理量（`tests/test_t15_reproducible.py`）。
- **自适应扫压的偏压序列**在失败边界附近对末位敏感，可复现的是"是否到达终点"，
  不是中间的尝试列表（`reports/week05_acceptance.md` §4.2）。

## 7. 常见错误

| 现象 | 原因 |
|---|---|
| `UnitError: ... declare units explicitly` | 字段缺 `_cm`/`_m`/`_cm3`/`_m3` 后缀，或与 `units.system` 不符 |
| `InvalidPnOrientation` | 左端掺杂不是 P 型（`C_left < 0`）或右端不是 N 型 |
| `right_voltage_V` 报错 | 必须为 `0.0`，偏压只加在左端 |
| `segmented mesh must list junction ... as a breakpoint` | 分段网格必须把结区坐标列为断点 |
| `FileExistsError` | 输出目录已有三文件，加 `--overwrite` |
| `decision: REJECTED` + `LINE_SEARCH_FAILED` | 继续步长太大，减小 `step_V` 或开 `adaptive` |
