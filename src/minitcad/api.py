"""B 的公共运行入口（§10.2）：`run()` 与 `write_required_outputs()`。

`run()` 的输入是 `DeviceConfig`（YAML 或 GUI 表单都落到它），输出 `RunResult`；
CLI 与 GUI 只通过这两个函数发起正式计算，保证同一配置得到同一结果与同一组文件。

取消语义：`should_stop()` 在每次尝试偏压之前检查；取消后已接受的偏压点完整保留，
`decision=REJECTED`、`failure_reason="CANCELLED"`（§10.2 规定取消结果不是 ACCEPTED）。
"""

from __future__ import annotations

import pathlib
from datetime import datetime
from math import isclose
from typing import Callable
from uuid import uuid4

from minitcad.core.device import DeviceConfig
from minitcad.io.config import config_fingerprint
from minitcad.simulation.run_result import (
    DECISION_ACCEPTED,
    DECISION_REJECTED,
    RunResult,
    STATUS_PASS,
)
from minitcad.simulation.sweep import equilibrium_start_state, run_fixed_sweep

Progress = Callable[[float, str], None]
ShouldStop = Callable[[], bool]

REQUIRED_FILES = ("comparison.csv", "comparison.png", "convergence.md")


def _never_stop() -> bool:
    """`should_stop` 的默认实现：不取消。用函数而不是 lambda，便于签名自省与序列化。"""
    return False


def _new_run_id(config: DeviceConfig) -> str:
    """`时间戳-短随机`；与 `config_sha256` 一起保证三文件可以互相追溯（§12.4）。"""
    stamp = datetime.now().strftime("%Y%m%dT%H%M%S")
    return f"{stamp}-{uuid4().hex[:8]}"


def run(
    config: DeviceConfig,
    progress: Progress | None = None,
    should_stop: ShouldStop = _never_stop,
) -> RunResult:
    """按 `config.sweep` 的偏压表求解，返回 `RunResult`。

    0 V 平衡解由 A 的标量 Poisson 提供，作为第一个偏压的热启动；若配置的起点不是 0 V，
    它仍作为初值，由第一个外层迭代带到目标偏压。

    判决只看**最终是否到达 `stop_V`**：固定步进下这等价于“每个偏压都通过”，自适应下
    允许中途回滚，但回滚的失败尝试仍留在 `attempts` 里（`§11.3` 不允许把失败掩盖成成功）。
    未到达则 `decision=REJECTED`，`failure_reason` 记录最后一个失败偏压与求解器的原因。
    """
    from minitcad.core.device import build_device
    from minitcad.simulation.sweep import run_adaptive_sweep

    device = build_device(config)
    initial = equilibrium_start_state(device)
    if config.sweep.adaptive:
        outcome = run_adaptive_sweep(
            device,
            config.sweep,
            initial,
            config.solver,
            progress=progress,
            should_stop=should_stop,
        )
    else:
        outcome = run_fixed_sweep(
            device,
            config.sweep.biases(),
            initial,
            config.solver,
            progress=progress,
            should_stop=should_stop,
        )

    attempts = outcome.attempts
    successful = tuple(attempt for attempt in attempts if attempt.status == STATUS_PASS)
    reached = outcome.reached_V
    if outcome.cancelled:
        decision, reason = DECISION_REJECTED, "CANCELLED"
    elif reached is not None and isclose(reached, config.sweep.stop_V, rel_tol=0.0, abs_tol=1e-12):
        decision, reason = DECISION_ACCEPTED, None
    else:
        decision = DECISION_REJECTED
        if not attempts:
            reason = "NO_BIAS_ATTEMPTED"
        else:
            last = attempts[-1]
            detail = last.result.reason if last.result is not None else "NO_RESULT"
            rollback = f"; {last.note}" if last.note else ""
            reason = f"{last.status} at {last.bias_V:.6g} V: {detail}{rollback}"

    return RunResult(
        config=config,
        device=device,
        attempts=attempts,
        successful=successful,
        decision=decision,
        failure_reason=reason,
        run_id=_new_run_id(config),
        config_sha256=config_fingerprint(config),
    )


def write_required_outputs(
    run_result: RunResult, path: str, overwrite: bool = False
) -> tuple[str, str, str]:
    """写固定的三个文件，返回 `(comparison.csv, comparison.png, convergence.md)` 的绝对路径。

    `overwrite=False` 时目录里已有同名文件就抛 `FileExistsError`：求解输出是可审计证据，
    静默覆盖会让报告里的 `run_id` 与文件内容对不上（§12.4）。
    """
    from minitcad.io.plot import write_comparison_png
    from minitcad.io.report import write_convergence_md
    from minitcad.io.results_csv import write_comparison_csv

    directory = pathlib.Path(path)
    existing = [name for name in REQUIRED_FILES if (directory / name).exists()]
    if existing and not overwrite:
        raise FileExistsError(
            f"{directory} already contains {existing}; pass overwrite=True to replace them"
        )
    directory.mkdir(parents=True, exist_ok=True)

    csv_path = write_comparison_csv(run_result, directory / REQUIRED_FILES[0])
    png_path = write_comparison_png(run_result, directory / REQUIRED_FILES[1])
    md_path = write_convergence_md(run_result, directory / REQUIRED_FILES[2])
    return tuple(str(p.resolve()) for p in (csv_path, png_path, md_path))  # type: ignore[return-value]
