# miniTCAD

一维硅 PN 结漂移-扩散模拟器（课程作业，8 周，A/B 两人协作）。
内核为**无量纲**有限体积离散：Scharfetter–Gummel 边通量 + 中能级 SRH + 唯一残差 `F(X)`；
求解器内部使用**准费米能级**变量 `(u, eta_n, eta_p)`（A 的 `b195ce5` 精度重构），
对外仍是物理状态 `State(u, n, p)`。

规范与任务清单见 `目的与规划.md`；B 的 W2 工作清单见 `docs/W2_B任务清单与开始方案.md`。

## 目录结构

```text
src/minitcad/
  core/            B：网格 Mesh1D、分段掺杂、cm-first 材料、接触、State/QFState、NormalizedDevice、常量
  io/              B：YAML 配置校验、cm↔SI 适配、comparison.csv long-form schema
  physics/         A：热电压、归一化尺度、中能级 SRH（+ B 的 cm 尺度装配）
  discretization/  A：Bernoulli、Poisson、SG、continuity、唯一残差 evaluate_qf/residual_qf
  simulation/      equilibrium 标量初值（A）、RunResult（B）
  cli.py           B：`minitcad validate`（可用）/ `minitcad solve`（W4 接线）
tests/             B 的独立 oracle（T01–T07、I-01）与 A 模块的交接回归
examples/          cm / SI 同器件 YAML、SI MMS 验证器
docs/              数学契约、A/B 接口会审、数值约定、convergence 模板
reports/           W2 / W3 验收记录（实测命令与数值）
C2_v1/             A 的独立参考树（不被根包导入）
ai_example/        单文件教学示例（与本规范无接口关系）
```

## 运行

```bash
python -m pytest -q                                  # 全量（262 passed）
PYTHONPATH=src python -m minitcad.cli validate examples/pn_left_p_right_n.yaml
PYTHONPATH=src python -m minitcad.cli validate examples/pn_left_p_right_n_si.yaml
PYTHONPATH=src python -m examples.verify_poisson     # SI MMS + 均匀正弦加密误差比
```

若未 `pip install -e .`，用 `PYTHONPATH=src` 或依赖 `tests/conftest.py` 的路径注入。

## 当前进度（B 侧）

| 周 | B 的交付 | 状态 | 证据 |
|---|---|---|---|
| W1 | 包骨架、`pyproject.toml`、示例 YAML、`.vscode`、convergence 模板 | 完成 | `pyproject.toml`、`docs/convergence_template.md` |
| W2 | `Mesh1D`、分段掺杂 `C̄`、材料/温度策略、接触、`State`、`NormalizedDevice`、cm↔SI、接口会审、公共 API 冻结 | 完成 | `reports/week02_acceptance.md`、`tests/test_public_api.py` |
| W3 | T01–T07 独立 oracle、反向 PN、`RunResult`+CLI 骨架、22 列 CSV schema、QF 口径核对 | 完成 | `reports/week03_acceptance.md` |
| W4 | `sweep`、`run()` 与 `solve` 接线、T08/T09 | 未开始 | 需 A 的 `solvers/gummel.py` 合并 |

**单位契约**：输入 `cm` 或 `si`（适配层换算成 cm）→ 内核无量纲 → writer 按 `output.unit_system` 输出。
**方向契约**：`P(left, x=0) -> N(right, x=L)`，`V_app = V_left - V_right`，`+V` 为正向，右端 0 V。

## 与 A 的同步状态

- 远端 `origin/main` = `98ff013`，已包含 A 的 `d53871e`（W3 物理/离散）与
  `b195ce5`（准费米能级重构）；B 的 W2/W3 在 QF 口径上复核通过。
- 待 A 处理：`C2_v1/tests/test_qf_stability.py` 4 项仍用旧 `residual_qf` 签名；
  `solvers/gummel.py` 内旧 `solve_gummel` 与 QF 版重复定义，W4 合并前须删旧留新。
  详见 `reports/week03_acceptance.md` §5 与 `docs/a_b_interface_week02.md` §5。
