# W2 验收记录（B 侧，QF 重构后复验）

> 对应 `docs/W2_B任务清单与开始方案.md` §7 的 DoD 清单。
> 复验日期：2026-10-25；解释器 Python 3.14.6；NumPy 2.5.3；SciPy 1.18.0；PyYAML 6.0.3。
> 复验原因：A 于 `d53871e`（2026-10-21）/`b195ce5`（2026-10-23）把求解器核心改成
> 准费米能级（QF）表示，`State/Contact/SolveResult` 的数据结构受影响；本记录是
> 同步 A 的 `b195ce5` 之后对 W2 交付的重新验收。

## 1. 复现命令

```bash
cd <repo-root>
python -m pytest -q                 # 235 passed
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
- 全量测试：**235 passed, 0 skipped**（W2 原有 118 项 + W3 oracle/QF/CLI/CSV 等）。

## 5. 未通过项 / 遗留

无 W2 阻断项。A 的 `C2_v1` 参考树自身仍有 4 项失败（旧 `residual_qf` 签名），
不影响根包；详见 `reports/week03_acceptance.md` §5 与
`docs/a_b_interface_week02.md` §5。
