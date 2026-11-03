# W2 验收记录（B 侧，QF 重构后复验）

> 对应 `docs/W2_B任务清单与开始方案.md` §7 的 DoD 清单。
> 复验日期：2026-10-25；加固日期：2026-11-03；解释器 Python 3.14.6；NumPy 2.5.3；SciPy 1.18.0；PyYAML 6.0.3。
> 复验原因：A 于 `d53871e`（2026-10-21）/`b195ce5`（2026-10-23）把求解器核心改成
> 准费米能级（QF）表示，`State/Contact/SolveResult` 的数据结构受影响；本记录是
> 同步 A 的 `b195ce5` 之后对 W2 交付的重新验收。
> 加固原因：2026-11-03 复核 §5.1/§5.3/§5.4 时发现三处校验缺口（分段网格结区未强制落结点、
> 结点数可被静默截断、右端偏压未强制为 0 V），已在本轮补齐，见 §6。

## 1. 复现命令

```bash
cd <repo-root>
python -m pytest -q                 # 253 passed
python -m pytest tests/test_core_mesh.py tests/test_core_region.py \
    tests/test_core_material.py tests/test_core_device.py \
    tests/test_core_contact.py tests/test_core_state.py \
    tests/test_core_qf_state.py tests/test_io_units.py tests/test_io_config.py -q
PYTHONPATH=src python -m minitcad.cli validate examples/pn_left_p_right_n.yaml
PYTHONPATH=src python -m examples.verify_poisson
# SI MMS & ND scaling: PASS
# Sine refinement error ratios (17/33/65/129 nodes): 4.00579125, 4.00144626, 4.00036147
```

## 2. DoD 清单

- [x] `src/minitcad/` 可导入，`pyproject.toml` 为 src 布局；`minitcad` 控制台入口已登记（W3 起）。
- [x] `Mesh1D` 提供 `uniform/segmented` 与 cm、nd 两套字段，`M≥3`、严格递增、结区在断点内（`tests/test_core_mesh.py`）。
- [x] `core/material.py` 为 cm-first，含 `_cm/_cm2/_cm3` 后缀命名与温度一致性校验（`tests/test_core_material.py`）。
- [x] `region.py` 支持分段掺杂与控制体平均 `C̄_i`，非法配置全部报错（`tests/test_core_region.py`）。
- [x] `io/units.py` 四类往返通过，无单位输入被拒绝（`tests/test_io_units.py`、`test_io_config.py`）。
- [x] `contact_targets` 与 `State` 可被 A 直接调用；`NormalizedDevice` 含 cm 基准尺度（`tests/test_core_contact.py`、`test_core_device.py`）。
- [x] T03、T04、I-01、cm↔nd、cm↔SI、异常注入全部通过。
- [x] `docs/a_b_interface_week02.md` 已按 QF 重构更新（B 已签，A 待签）。
- [x] `docs/numerical_conventions.md` 已建立（§14.7 要求与接口会审记录配套的单位/符号/容差契约）。
- [x] `docs/convergence_template.md`（§12.3 十章节模板）与 `.vscode/settings.json`（§14.5）已补齐。
- [x] `reports/week02_acceptance.md` 即本文件；W3 记录见 `reports/week03_acceptance.md`。

## 3. QF 重构带来的 W2 变更

| 变更 | 位置 | 兼容性 |
|---|---|---|
| 新增 `QFState(u, eta_n, eta_p, ref_n, ref_p)`，只读、3M 交错、M≥3 | `src/minitcad/core/state.py` | 新增，不影响 `State` 调用 |
| 新增 `PhysicsDiagnostics(fluxes, recombination_nd, residual_qf_nd)` | 同上 | 新增 |
| `SolveResult` 末尾新增 `physics: PhysicsDiagnostics \| None = None` | 同上 | 位置参数兼容；仅 `converged=True` 可带 |
| 新增 `contact_targets_qf(device, bias_V, *, ref_n, ref_p)` | `src/minitcad/core/contact.py` | 新增；物理 `contact_targets` 不变 |
| 新增 `tests/test_core_qf_state.py`（12 项） | `tests/` | QF 往返、只读、3M、边界、诊断形状 |

## 4. 数值证据（摘录）

- 接触：`np=1`、`p-n+C=0`、`u_L=a(C_L)+v` 误差 ≤ `1e-10`（T03/T04）。
- 接触 QF 反解：由 `(u,eta_n,eta_p)` 重建的 `n,p` 与物理 `contact_targets` 相对误差 ≤ `1e-14`。
- `NormalizedDevice.L0_cm` 与附录 C 一致（`tests/test_core_device.py`）。
- 全量测试：**253 passed, 0 skipped**（2026-10-25 基线 235 项 + 2026-11-03 新增 18 项校验/契约/精度测试）。

## 5. 未通过项 / 遗留

无 W2 阻断项。A 的 `C2_v1` 参考树自身仍有 4 项失败（旧 `residual_qf` 签名），
不影响根包；详见 `reports/week03_acceptance.md` §5 与
`docs/a_b_interface_week02.md` §5。

## 6. 2026-11-03 加固：三处被静默接受的非法配置

复核 §5.1/§5.3/§5.4 时发现根包（B 的单一生产包）与 A 在 `C2_v1` 中的加强校验存在差异，
以下三项此前**不会**失败，现在全部显式拒绝，合法配置的数值结果不变。

| 缺口 | 依据 | 现在的行为 | 测试 |
|---|---|---|---|
| `mesh.type: segmented` 未强制结区落结点 | §5.1「结区坐标必须包含在网格断点中」、W2-B1 DoD | `build_device` 要求每个内部 region 边界既在 `breakpoints_cm` 中、又落在 `x_nodes_cm` 上（相对容差 `1e-12`），否则报 `ValueError` 并列出结区坐标与最近结点 | `test_core_device.py::test_segmented_mesh_must_include_junction_breakpoint`、`::test_segmented_mesh_puts_the_junction_on_a_node`、`test_io_config.py::test_segmented_yaml_without_the_junction_is_rejected` |
| 结点数可被静默截断（`nodes: 20.5` → 20） | §5.4 | `uniform`/`segmented` 与 YAML 解析层都要求整数结点数（拒绝 `20.5`、`True`） | `test_core_mesh.py::test_non_integer_node_count_is_rejected`、`test_core_device.py::test_mesh_spec_requires_integer_node_counts`、`test_io_config.py::test_non_integer_node_count_in_yaml_is_rejected` |
| `contacts.right_voltage_V` 非 0 未被拒绝 | §4.2「右端 0 V 参考」 | `DeviceConfig` 与配置解析层都拒绝非零右端电压，并说明偏压只由左端给出 | `test_core_device.py::test_right_terminal_must_be_grounded`、`test_io_config.py::test_nonzero_right_voltage_is_rejected` |
| `w_interior` 与相邻两条边不匹配时未校验 | §5.1 `w_i=(h_{i-1}+h_i)/2` | `Mesh1D` 构造时校验 `w = (h[:-1]+h[1:])/2`，错取下标显式失败 | `test_core_mesh.py::test_inconsistent_control_volume_width_is_rejected` |
| 接触浓度在 float64 边界上未拒绝 | §4.2 非法状态 | `contact_concentrations_nd` 对非有限 `C` 与不可表示结果报 `ValueError`；极端有限 `C`（`1.5e308`）实测仍满足 `np=1` | `test_core_contact.py::test_extreme_contact_doping_is_rejected_or_stays_representable` |

`uniform` 网格的结区允许落在控制体内部（由 `cbar_nd` 的控制体积分吸收），
这条差异已用 `test_core_device.py::test_uniform_mesh_keeps_junction_inside_a_control_volume` 明确固定。
