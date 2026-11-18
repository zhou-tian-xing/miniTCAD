# GUI ↔ 公共 API 契约（`w07-core-freeze` 冻结）

> A 在 W8 写 PySide6 GUI 时只需依赖本文件列出的东西。任何一项不够用，先改本文件并通知
> B，**不要在 GUI 里另写一套求解或输出**。
> 契约测试：`tests/test_w07_freeze.py`、`tests/test_cancel_contract.py`、`tests/test_io_writer.py`。

## 1. GUI 只调这五个函数

```python
import minitcad

config = minitcad.load_config(path)          # 或 parse_config_dict(dict)  ← 表单改名走这里
config = minitcad.parse_config_dict(data)    # 已解析的 YAML 字典 / 面板值
minitcad.save_config(config, target_path)    # 保存；读回来指纹相同

result = minitcad.run(config, progress=cb, should_stop=stop)   # 后台线程里跑
paths = minitcad.write_required_outputs(result, out_dir, overwrite=True)
```

- `config_sha256` 可以直接算：`minitcad.config_fingerprint(config)`，用于界面上显示
  "这次运行的配置指纹"。
- **不要在 GUI 里构造 `DeviceConfig`**：字段校验（单位后缀、结区落结点、右端 0 V）
  都在 `parse_config_dict` 里，绕过它就等于绕过 §5 的全部校验。
- 不要在 GUI 里读 CSV 再自己算电流；用 `result.attempts[i].result.physics`。

## 2. `run()` 的线程规则

`run()` 是**同步阻塞**的，必须放在 `QThread`/`QThreadPool` 里跑，否则窗口会卡住。

```python
class Runner(QThread):
    progressed = Signal(float, str)          # (bias_V, "PASS"|"FAIL")
    finished_run = Signal(object)            # RunResult

    def __init__(self, config):
        super().__init__()
        self.config = config
        self._stop = False

    def request_stop(self):
        self._stop = True

    def run(self):
        self.finished_run.emit(
            minitcad.run(self.config,
                         progress=lambda bias, status: self.progressed.emit(bias, status),
                         should_stop=lambda: self._stop)
        )
```

约束：

1. `progress` 在**工作线程**里被调用，只能通过信号转到界面，不能直接碰控件。
2. `should_stop` 是协作式的：只在两次尝试**之间**检查，正在进行的求解不会被打断。
   所以取消的响应延迟 ≈ 一个偏压点的求解时间（本机毫秒~秒级），不是"立刻"。
3. 取消后 `decision == "REJECTED"`、`failure_reason == "CANCELLED"`，**已接受的偏压点与
   三个文件仍然可用**——不要因为取消就清空结果页。

## 3. 结果页要展示什么

| 界面元素 | 数据来源 |
|---|---|
| 最终判决（大字） | `result.decision`，颜色取 `result.accepted` |
| 失败原因 | `result.failure_reason`（`None` 时留空） |
| 偏压表 | `result.attempts`：`bias_V`、`status`、`result.iterations`、`result.reason`、`note` |
| I–V 曲线 | `attempt.result.physics.fluxes`，终端取最后一条边 |
| 本征/派生量 | `Attempt` → `simulation.diagnostics.fields_for_attempt(config, device, attempt)` |
| 三文件 | `write_required_outputs` 的返回值（绝对路径，可直接给"打开文件夹"） |
| 追溯信息 | `result.run_id`、`result.config_sha256` |

`decision` 与三个文件必须一致：`convergence.md` 的 `## 10. Final decision` 就是
`result.decision`，界面不能显示成与报告不同的结论。

## 4. 失败与取消的三条路径

| 路径 | `decision` | `failure_reason` | 文件 |
|---|---|---|---|
| 成功到达 `stop_V` | `ACCEPTED` | `None` | 三文件 |
| 停在失败点 | `REJECTED` | `FAIL at <V> V: <reason>` | 三文件（含 `FAIL` 行） |
| 用户取消 | `REJECTED` | `CANCELLED` | 三文件（只含已接受点） |

三者在 GUI 上都应能"另存/打开输出目录"——**取消和失败也要能交文件**，这是 §12.4 的
要求，不是可选项。

## 5. 写文件

```python
write_required_outputs(result, str(out_dir), overwrite=True)
```

- 固定三个文件名：`comparison.csv`、`comparison.png`、`convergence.md`；
- `overwrite=False` 时目标目录已有同名文件会抛 `FileExistsError`——GUI 应先问用户；
- 目录不存在会自动创建（含父目录）；
- 返回 `(csv, png, md)` 的绝对路径，顺序固定。

## 6. 单位

- 输入由 YAML 的 `units.system` 声明；GUI 面板默认 cm，要支持 SI 就把 `_cm` 字段换成
  `_m` 变体（见 `docs/user_guide.md` §4）。
- 输出单位由 `output.unit_system` 决定，**只影响 writer**，不影响任何物理量。
- `save_config` 一律按 cm 写；用 SI 打开的配置保存后 `input_unit_system` 变为 `cm`
  （物理不变，指纹随之改变，这是有意的）。

## 7. W8 的验收对着这张表

| 测试 | 用本文件的哪几条 |
|---|---|
| T16 GUI 配置 | §1（只调五个函数）、§6（单位往返）、§1 的 `save_config` 往返等价 |
| T17 GUI 运行/取消 | §2（线程与取消语义）、§3（结果页数据来源）、§4（三条路径）、§5（写文件） |
| GUI/CLI 一致性 | 同一 YAML 下 `config_sha256`、`attempts`、终端电流、三文件 schema 必须与 CLI 完全相同 |
