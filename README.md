# miniTCAD

一维硅 PN 结漂移-扩散模拟器（课程作业，8 周，A/B 两人协作）。
内核为**无量纲**有限体积离散：Scharfetter–Gummel 边通量 + 中能级 SRH + 唯一残差 `F(X)`；
求解器内部使用**准费米能级**变量 `(u, eta_n, eta_p)`（A 的 `b195ce5` 精度重构），
对外仍是物理状态 `State(u, n, p)`。

规范与任务清单见 `目的与规划.md`；B 的后续计划见 `docs/B_roadmap.md`。

## 目录结构

```text
src/minitcad/
  core/            B：网格 Mesh1D、分段掺杂、cm-first 材料、接触、State/QFState、NormalizedDevice
  io/              B：YAML 配置校验、cm↔SI 适配、comparison.csv/PNG/收敛报告 writer
  physics/         A：热电压、归一化尺度、中能级 SRH（+ B 的 cm 尺度装配）
  discretization/  A：Bernoulli、Poisson、SG、continuity、唯一残差 evaluate_qf/residual_qf
  solvers/         A：Gummel + 收敛判据（Newton 为 W5）
  simulation/      B：dc 单点编排、sweep 扫压、diagnostics 物理量展开；A：equilibrium 初值
  api.py           B：run() / write_required_outputs()（GUI 与 CLI 的唯一入口）
  cli.py           B：minitcad validate / solve
tests/             B 的独立 oracle（T01–T09）、API/CLI/writer 契约、基准快照
examples/          cm / SI 同器件 YAML、W4 基准配置与数据快照、SI MMS 验证器
docs/              数学契约、A/B 接口会审、数值约定、convergence 模板、B 路线图
reports/           W2 / W3 / W4(Gate 1) 验收记录
C2_v1/             A 的独立参考树（不被根包导入）
ai_example/        单文件教学示例（与本规范无接口关系）
```

## 运行

```bash
python -m pytest -q                              # 全量（291 passed）
python -m pytest tests/test_week04_gate1.py -q   # T08/T09

PYTHONPATH=src python -m minitcad.cli validate examples/pn_baseline_0_to_030.yaml
PYTHONPATH=src python -m minitcad.cli solve examples/pn_baseline_0_to_030.yaml \
    --out runs/pn_baseline_0_to_030
# decision: ACCEPTED；写出 comparison.csv / comparison.png / convergence.md
PYTHONPATH=src python -m minitcad.cli validate runs/pn_baseline_0_to_030/comparison.csv

PYTHONPATH=src python -m examples.verify_poisson   # SI MMS + 正弦加密误差比
```

若未 `pip install -e .`，用 `PYTHONPATH=src` 或依赖 `tests/conftest.py` 的路径注入。
PNG 需要 matplotlib（`pyproject.toml` 的 `dev` extras 已包含）。

## 当前进度（B 侧）

| 周 | B 的交付 | 状态 | 证据 |
|---|---|---|---|
| W1 | 包骨架、`pyproject.toml`、示例 YAML、`.vscode`、convergence 模板 | 完成 | `pyproject.toml`、`docs/convergence_template.md` |
| W2 | `Mesh1D`、分段掺杂 `C̄`、材料/温度策略、接触、`State`、`NormalizedDevice`、cm↔SI、接口会审、公共 API 冻结 | 完成 | `reports/week02_acceptance.md`、`tests/test_public_api.py` |
| W3 | T01–T07 独立 oracle、反向 PN、`RunResult`+CLI 骨架、22 列 CSV schema、QF 口径核对 | 完成 | `reports/week03_acceptance.md` |
| W4 | T08/T09、基准配置+快照、固定步进扫描、`run()`/CLI 接线、三文件 writer 初版 | 完成（Gate 1 = PASS） | `reports/week04_gate1.md` |
| W5 | FD Jacobian oracle、自适应扫压与回滚、取消契约冻结 | 未开始（依赖 A 的 Newton） | `docs/B_roadmap.md` §1 |
| W6 | T10/T11/T13、双求解器对比脚本、细网格参考、writer 补齐 | 未开始 | `docs/B_roadmap.md` §2 |
| W7 | refine 网格、T12/T14/T15、schema/GUI API 冻结、CI | 未开始 | `docs/B_roadmap.md` §3 |
| W8 | GUI/CLI 一致性、T16/T17、Windows 打包与 T18 | 未开始 | `docs/B_roadmap.md` §4 |

**单位契约**：输入 `cm` 或 `si`（适配层换算成 cm）→ 内核无量纲 → writer 按 `output.unit_system` 输出。
**方向契约**：`P(left, x=0) -> N(right, x=L)`，`V_app = V_left - V_right`，`+V` 为正向，右端 0 V。

## 与 A 的同步状态

- 远端 `origin/main` 已包含 A 的 `d53871e`（W3 物理/离散）与 `b195ce5`（准费米能级重构）；
  B 的 W2–W4 在 QF 口径上复核通过，`C2_v1` 参考树 `168 passed`。
- 待 A 处理：`solvers/newton.py`（W5-A1…A6）；`area_cm2` 可空与 `normalization.py`
  命名分叉需要 A 给出正式口径（见 `docs/a_b_interface_week02.md` §5）。
