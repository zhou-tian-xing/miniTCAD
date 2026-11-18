# W7 验收记录：非均匀网格终审、三文件验收与核心冻结

> 对应 `目的与规划.md` §14.2（W7-B1…W7-B5）与 §14.3 的 T12/T14/T15。
> 日期：2026-11-18；Python 3.14.6 / NumPy 2.5.3 / SciPy 1.18.0 / matplotlib 3.11.1。
> 交付物：本报告 + tag `w07-core-freeze`。
> **注意**：本轮按用户安排**先做 W7、W6 顺延**，Gate 2（T11 双求解器正式对照）仍未裁决，
> 见 §8。

## 1. 复现命令

```bash
cd <repo-root>
python -m pytest -q                                   # 385 passed
python -m pytest tests/test_t12_mesh.py tests/test_t14_output.py \
    tests/test_t15_reproducible.py tests/test_w07_freeze.py -q    # 52 passed
python -m minitcad.cli solve examples/pn_baseline_0_to_030.yaml --out runs/w07
python scripts/accept_outputs.py runs/w07             # 17/17 passed
python scripts/run_mesh_convergence.py                # T12 两张表
cd C2_v1 && python -m pytest -q                       # 215 passed（A 的参考树）
```

## 2. W7-B1：结区 refine 网格生成器

`core/mesh.py` 新增 `refined(length_cm, junction_cm, nodes, l0_cm, *, refine=20.0)`，
YAML 侧是 `mesh: {type: refined, nodes: 201, refine: 20}`：

- **结点总数固定**（= `nodes`），可与同结点数的 `uniform` 直接比较（§15 实验 C 的前提）；
- 结区**必定落在结点上**；`refine=1` 与 `uniform` 逐位相同（对照组因此是干净的）；
- `refine` 定义为"远端单元 / 结区单元"的长度比。**不用公比**：`M≈200` 时公比 `1.15`
  对应 1000 倍跨度，而 `1.03` 才对应 20 倍，按公比容易选错（本项第一版就踩了这个坑，
  生成出 `h_min = 8.8e-14 cm` 的病态网格，已改）；
- 非法参数（`refine<1`、结区不在内部、结点数不足）一律显式失败。

## 3. W7-B2：T12 / T14 / T15

### 3.1 T12 网格收敛（明细见 `reports/week07_mesh_convergence.md`）

参考量：uniform 801，`J = 4.249672901e-05 A/cm²`，`V_app = 0.30 V`，Gummel。

| mesh | M | \|dJ\|/J | \|dE\|/E |
|---|---|---|---|
| uniform | 201 | 1.646e-03 | 4.803e-02 |
| uniform | 401 | 2.116e-04 | 1.697e-02 |
| uniform | 801 | 参考 | 参考 |
| refined (refine=10) | 201 | **3.442e-05** | **1.462e-04** |
| refined (refine=20) | 201 | **1.103e-05** | 6.800e-03 |
| refined (refine=50) | 201 | 3.135e-05 | 1.242e-02 |

- 电流误差随加密**单调下降**，最细两级相对差 `2.12e-4 ≤ 1%` → **T12 判据通过**；
- 结点数不变时，结区加密把电流误差压到 `1.1e-5`（**150 倍**），耗时不变；
- `refine=50` 反而回退：结点总数固定，堆在结区越多、体内越粗，而正向电流由体内 SRH
  复合积分主导。**`refine≈20` 是这份器件在 201 结点下的平衡点**，不是越大越好；
- **峰值电场不单调**（`refine=20/50` 的 `E_peak` 超过 801 参考值）：掺杂是阶跃，电场峰
  在真值上是奇异的，加密只会让它继续长。T12 要求"全部差值入报告"，不要求它收敛；
  这里如实记录，并明确**不能拿 `E_peak` 当收敛判据**。

### 3.2 T14 三文件验收

`scripts/accept_outputs.py`（17 项检查）**不 import 任何 writer**，把冻结的 schema 直接
写死在验收器里，否则改 writer 时验收器会跟着一起改，等于没验。

覆盖：三文件存在非空、PNG 头可解且尺寸非零、报告十节齐、22 列表头（cm/SI 两版）、
`run_id`/`config_sha256`/`solver` 全表一致、每个 `record_type` 的 `quantity`/`unit` 白名单、
理论/误差列自洽（`abs=|num-theory|`、非零理论才有 `rel`、`none` 必须留空）、
每个收敛偏压的 node/edge 行数（`M` / `M-1`）、迭代行三元组完整、判决与最后一次尝试一致、
`FAIL` 行必须带原因。

| 场景 | 结果 |
|---|---|
| cm + Gummel + 固定步进 | 17/17 |
| SI 输出（`position_m` 表头） | 17/17 |
| Newton + 自适应（中途 4 次回滚） | 17/17 |
| Newton 固定步进（停在 0.1 V，`REJECTED`） | 18/18（含 `--expect REJECTED`） |

**验收器发现了两个真实的 schema 缺陷**（都已修）：

1. 结构行的 `theory_kind` 留空而不是 `none` → "没有理论"出现两种写法，验收器无法区分
   "本行无理论"与"writer 忘了填"。已改为 `_blank_row` 缺省 `theory_kind="none"`。
2. `applied_voltage_V` 在 0 V 时写了 `theory_rel_error=0.0`，而理论值也是 0 → `0/0`
   被写成了 0。已改为理论值为 0 时留空。

`tests/test_t14_output.py` 另外逐条**伪造损坏产物**（缺文件、坏 PNG、少结点行、错单位、
改数值不动误差、状态与判决矛盾），确认每种违规都会被抓到——只测"正常产物能通过"
等于没测验收器。

### 3.3 T15 可复现性

`tests/test_t15_reproducible.py`：同一配置连跑两次，`config_sha256`、判决、偏压点、
迭代数、残差、以及 `u/n/p` 三元组**逐位相同**；CSV 除 `run_id` 与 `runtime_s` 外逐行相同；
报告含平台/Python/判决元数据。

自适应扫压**不做**逐位复现断言：其尝试序列在失败边界附近对 `bias_V` 末位敏感
（`reports/week05_acceptance.md` §4.2），可复现的是"是否到达终点"。

## 4. W7-B3：三文件 schema 冻结

冻结内容写进验收器与 `docs/api_reference.md` §4：三个文件名、22 列顺序、各
`record_type` 的 `quantity`/`unit` 白名单、理论/误差列规则、报告十节标题。
`w07-core-freeze` 之后改任何一项都必须同时改验收器与接口会审记录。

## 5. W7-B4：CLI / API / YAML 冻结

| 冻结项 | 证据 |
|---|---|
| 包级 6 个导出与签名 | `tests/test_w07_freeze.py::test_package_exports_are_frozen`、`test_public_signatures_are_frozen` |
| `RunResult`/`Attempt`/`SweepConfig`/`MeshSpec` 字段表 | `test_result_fields_are_frozen` |
| CLI 退出码与首屏方向提示 | `test_cli_exit_codes_are_frozen` |
| YAML 未知字段拒绝 | `test_yaml_schema_rejects_unknown_fields` |
| 配置指纹覆盖全部求解参数 | `test_fingerprint_covers_the_mesh_and_sweep_options` |
| `save_config` 往返等价 | `test_save_config_round_trips_to_the_same_fingerprint`（cm）、`..._normalises_si_input_to_cm` |

本轮新增的公共函数是 `save_config`（W8 的 GUI 要用），并按 §7 的冻结清单把它纳入
`minitcad.__all__`。

**指纹缺陷修复（重要）**：`config_payload` 此前漏了 `mesh.refine` 与
`sweep.adaptive/min_step_V/max_step_V/grow_factor/shrink_factor`，即"只改自适应开关"的
两次运行会共用同一个 `config_sha256`。已在冻结前补全，**代价是所有旧指纹失效**
（基准快照已重算，`examples/baseline_0_to_030_snapshot.json` 里记了原因）。

文档：`docs/api_reference.md`、`docs/user_guide.md`；CI 见 `.github/workflows/ci.yml`
（Ubuntu + Windows × Python 3.10/3.12，跑根包、`C2_v1`、基准三文件验收与 T12 表）。

## 6. W7-B5：GUI API 契约

`docs/gui_api_contract.md` 冻结了 A 在 W8 需要的全部内容：只调五个函数、`run()` 的线程
规则（同步阻塞，必须放 `QThread`）、结果页每个元素的数据来源、成功/失败/取消三条路径
各自的三文件与判决、写文件的返回顺序、单位约定，以及 T16/T17 各自对应本文件的哪几条。

取消了"取消要立刻生效"的误解：`should_stop` 只在两次尝试之间检查，一次尝试的延迟是
毫秒~秒级。

## 7. 本轮改动清单

| 文件 | 改动 |
|---|---|
| `core/mesh.py` | 新增 `refined()`（结区加密、结点数固定、`refine=1` 等价 `uniform`） |
| `core/device.py` | `MeshSpec` 增 `refine` 与 `kind="refined"`；`build()` 接收结区坐标；结区落结点校验扩展到 refined |
| `io/config.py` | 解析 `mesh.type=refined`；`save_config` 新增；`config_payload` 补 `mesh.refine` 与自适应五项 |
| `io/results_csv.py` | `theory_kind` 缺省 `none`；0 V 的 `theory_rel_error` 留空 |
| `verification/reference.py` | 新增：网格样本、误差表、单调性判据 |
| `scripts/accept_outputs.py` | 新增：§12.4 独立验收器（17 项） |
| `scripts/run_mesh_convergence.py` | 新增：T12 表格生成 |
| `docs/{api_reference,user_guide,gui_api_contract}.md` | 新增：冻结文档 |
| `tests/test_{t12_mesh,t14_output,t15_reproducible,w07_freeze}.py` | 新增 52 项测试 |
| `examples/baseline_0_to_030_snapshot.json` | 重算 `config_sha256` 并记录原因 |
| `minitcad/__init__.py` | 导出 `load_config/parse_config_dict/save_config/config_fingerprint` |

## 8. 未完成与顺延项

1. **Gate 2 已补做（2026-11-18）**：T11 在 0.02 V 梯子上 15/15 通过（势峰值差 `3.7e-13 V`、
   电流相对差 `5.3e-10`），Gate 2 = PASS，见 `reports/week06_gate2.md`。
   冻结本身不依赖它，但交付时两者都可引用。
2. **Newton 在 0.10–0.15 V 失败，原因未查清**（已证伪“台阶即原因”“Δv 太大”两个假设），
   `docs/a_b_interface_week02.md` §5 第 7 项。
3. **A 的 `reverse_sweep.py` 未并入**根包，W6 若要覆盖反向偏压需先定口径。
4. **T05 的 W7 变网格终审**属 A 的 W7-A1，B 侧只提供 `refined` 网格与验收器。

## 9. 结论

- W7-B1…W7-B5 全部交付，测试 `385 passed`（W5 收口 333），`C2_v1` `215 passed`。
- T12 = **PASS**（电流单调收敛、最细两级 `2.12e-4`、差值完整入报告；`E_peak` 非单调
  已解释）；T14 = **PASS**（4 类场景各 17–18 项检查，并修复了它发现的两处 schema 缺陷）；
  T15 = **PASS**（逐位复现，自适应路径按设计不做逐位断言）。
- 核心与公共接口冻结，tag `w07-core-freeze`。W8 只允许 GUI、联调与阻塞性 bug 修复。
- **W6 / Gate 2 仍欠**，需在 W8 之前补做或明确记录为未通过。
