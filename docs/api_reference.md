# 公共 API 参考（`w07-core-freeze` 冻结）

> 冻结范围：`minitcad` 包的 6 个顶层函数、`RunResult`/`Attempt`/`SweepConfig`/`MeshSpec`
> 的字段表、CLI 退出码、YAML 字段名。冻结后只允许修阻塞性 bug；改名或改顺序必须走
> `docs/a_b_interface_week02.md` 的变更记录。
> 契约测试：`tests/test_w07_freeze.py`、`tests/test_public_api.py`。

## 1. 顶层导出

```python
import minitcad
```

| 名字 | 说明 |
|---|---|
| `minitcad.run(config, progress=None, should_stop=_never_stop)` | 跑一次完整扫压，返回 `RunResult` |
| `minitcad.write_required_outputs(run_result, path, overwrite=False)` | 写出三个固定文件，返回绝对路径 |
| `minitcad.load_config(path)` | 读 YAML → `DeviceConfig` |
| `minitcad.parse_config_dict(data)` | 读已解析的字典 → `DeviceConfig`（GUI 面板走这里） |
| `minitcad.save_config(config, path)` | `DeviceConfig` → 规范化 YAML（GUI 保存走这里） |
| `minitcad.config_fingerprint(config)` | 规范化配置的 SHA-256 |

`config` 是 `minitcad.core.device.DeviceConfig`；GUI 不需要自己构造它，用
`parse_config_dict`/`load_config` 即可。

## 2. `run()`

```python
run(config, progress=None, should_stop=<默认不取消>) -> RunResult
```

- 偏压表来自 `config.sweep`：`adaptive=False` 时按 `step_V` 等差；`adaptive=True` 时由
  继续法动态生成。
- `progress(bias_V, status)`：**每次尝试之后**回调，`status ∈ {"PASS", "FAIL"}`。
  它表示"这个偏压试完了"，不是进度百分比。
- `should_stop()`：每次尝试**之前**检查一次。取消是协作式的，正在进行的求解不会被打断。
- 0 V 平衡解由 A 的标量 Poisson 提供，作为第一个偏压的热启动。

判决只看**是否到达 `sweep.stop_V`**：

| 情形 | `decision` | `failure_reason` |
|---|---|---|
| 到达终点 | `"ACCEPTED"` | `None` |
| 被 `should_stop` 取消 | `"REJECTED"` | `"CANCELLED"` |
| 停在失败点 | `"REJECTED"` | `"FAIL at 0.1 V: LINE_SEARCH_FAILED"` 形式 |

自适应模式下中途回滚的 `FAIL` 尝试**仍留在 `attempts` 里**，但它们不阻止 `ACCEPTED`。

## 3. `RunResult`

| 字段 | 类型 | 说明 |
|---|---|---|
| `config` | `DeviceConfig` | 输入配置（原对象） |
| `device` | `NormalizedDevice` | 归一化后的器件（网格、尺度、掺杂） |
| `attempts` | `tuple[Attempt, ...]` | 按尝试顺序，含失败与回滚 |
| `successful` | `tuple[Attempt, ...]` | 只含 `PASS` |
| `decision` | `str` | `"ACCEPTED"` / `"REJECTED"` |
| `failure_reason` | `str \| None` | 见上表 |
| `run_id` | `str` | `时间戳-短随机`，写进三个文件 |
| `config_sha256` | `str` | 规范化配置指纹 |

`Attempt`：`bias_V`、`status`、`runtime_s`、`result`（`SolveResult`）、`note`
（扫描层的回滚说明，与 `result.reason` 分工）。

`SolveResult`：`bias_V`、`state`、`converged`、`iterations`、`residual_norm`、
`update_norm`、`reason`、`iteration_history`、`physics`。
失败时 `state is None`、`physics is None`。

## 4. `write_required_outputs()`

```python
write_required_outputs(run_result, path, overwrite=False) -> tuple[str, str, str]
```

按固定顺序写 `comparison.csv`、`comparison.png`、`convergence.md` 并返回绝对路径。
目录里已有同名文件且 `overwrite=False` 时抛 `FileExistsError`——求解输出是可审计证据，
静默覆盖会让报告里的 `run_id` 与文件内容对不上。验收用
`python scripts/accept_outputs.py <dir>`。

## 5. 异常

| 情形 | 异常 |
|---|---|
| YAML/字典字段非法、单位混用、未知字段 | `minitcad.io.units.UnitError`（`ValueError` 子类） |
| 器件/材料/网格/区域校验失败 | `ValueError`（`InvalidPnOrientation` 等） |
| 输出目录已有文件 | `FileExistsError` |
| 求解失败 | **不抛异常**，看 `RunResult.decision` |

## 6. CLI 退出码

| 码 | 含义 |
|---|---|
| 0 | `ACCEPTED` |
| 2 | `REJECTED`（求解失败或配置非法） |
| 3 | 保留：`run()` 抛 `NotImplementedError`（当前无此路径） |

CLI 首屏固定打印方向/极性：`P(left, x=0) -> N(right, x=L) | V_app = V_left - V_right |
+V = forward bias | right = 0 V`。
