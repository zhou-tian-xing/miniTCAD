"""T12：独立细网格参考与网格收敛表（§15 实验 C）。

对同一器件、同一偏压，只改变结点的**空间分布**（总结点数可不同也可以相同），比较
终端电流 `I`、势垒与峰值电场，并用最细网格作为参考量出各级误差。

"独立"的含义：这里不调用 writer、不读 CSV，直接消费 `api.run()` 的 `RunResult` 与
`solve_*` 返回的 `PhysicsDiagnostics`，与三文件链路互不依赖。
"""

from __future__ import annotations

import dataclasses
import time
from dataclasses import dataclass

import numpy as np

from minitcad.api import run
from minitcad.core.device import DeviceConfig, MeshSpec

# 相对差的判据（T12：最细两级非零 `I` 的相对差目标 ≤1%）。
DEFAULT_TOLERANCE = 1e-2


@dataclass(frozen=True, slots=True)
class MeshSample:
    """一次网格-偏压实验的结果。"""

    label: str
    kind: str
    nodes: int
    refine: float | None
    bias_V: float
    converged: bool
    current_density_A_cm2: float
    current_A: float
    barrier_V: float
    peak_field_V_cm: float
    iterations: int
    runtime_s: float
    min_spacing_cm: float = 0.0


def _config_with_mesh(config: DeviceConfig, mesh: MeshSpec) -> DeviceConfig:
    return dataclasses.replace(config, mesh=mesh)


def sample(
    config: DeviceConfig,
    mesh: MeshSpec,
    *,
    label: str,
    bias_V: float | None = None,
) -> MeshSample:
    """在同一配置上换一份网格跑一次扫压，取 `bias_V`（默认 `sweep.stop_V`）处的量。"""
    target = float(config.sweep.stop_V if bias_V is None else bias_V)
    started = time.perf_counter()
    result = run(_config_with_mesh(config, mesh))
    elapsed = time.perf_counter() - started

    attempt = None
    for candidate in result.successful:
        if abs(candidate.bias_V - target) <= 1e-12:
            attempt = candidate
    if attempt is None or attempt.result is None or attempt.result.state is None:
        return MeshSample(
            label=label,
            kind=mesh.kind,
            nodes=result.device.mesh.node_count,
            refine=mesh.refine if mesh.kind == "refined" else None,
            bias_V=target,
            converged=False,
            current_density_A_cm2=float("nan"),
            current_A=float("nan"),
            barrier_V=float("nan"),
            peak_field_V_cm=float("nan"),
            iterations=0,
            runtime_s=elapsed,
            min_spacing_cm=float("nan"),
        )

    device = result.device
    state = attempt.result.state
    fluxes = attempt.result.physics.fluxes
    currents = fluxes.jn_edges + fluxes.jp_edges
    scales = device.scales

    barrier_V = float(np.max(state.u) - np.min(state.u)) * scales.v0_v
    field = np.abs(np.diff(state.u)) / device.mesh.h_edges_cm * scales.v0_v
    return MeshSample(
        label=label,
        kind=mesh.kind,
        nodes=device.mesh.node_count,
        refine=mesh.refine if mesh.kind == "refined" else None,
        bias_V=target,
        converged=True,
        current_density_A_cm2=float(currents[-1]) * scales.j0_a_cm2,
        current_A=float(currents[-1]) * scales.j0_a_cm2 * float(config.area_cm2 or 0.0),
        barrier_V=barrier_V,
        peak_field_V_cm=float(np.max(field)),
        iterations=int(attempt.result.iterations),
        runtime_s=elapsed,
        min_spacing_cm=float(device.mesh.h_edges_cm.min()),
    )


def relative_difference(value: float, reference: float) -> float:
    """`|v - ref| / max(|ref|, floor)`；参考量接近零时退化为绝对差（T12 的“接近零改报绝对差”）。"""
    if not np.isfinite(value) or not np.isfinite(reference):
        return float("inf")
    scale = abs(reference)
    if scale < 1e-12:
        return abs(value - reference)
    return abs(value - reference) / scale


def convergence_table(samples: list[MeshSample], reference: MeshSample) -> list[dict]:
    """把每个样本与参考量的差整理成表行。"""
    rows = []
    for item in samples:
        rows.append(
            {
                "label": item.label,
                "kind": item.kind,
                "nodes": item.nodes,
                "refine": item.refine,
                "bias_V": item.bias_V,
                "converged": item.converged,
                "current_density_A_cm2": item.current_density_A_cm2,
                "barrier_V": item.barrier_V,
                "peak_field_V_cm": item.peak_field_V_cm,
                "iterations": item.iterations,
                "runtime_s": item.runtime_s,
                "current_rel_error": relative_difference(
                    item.current_density_A_cm2, reference.current_density_A_cm2
                ),
                "barrier_rel_error": relative_difference(item.barrier_V, reference.barrier_V),
                "peak_field_rel_error": relative_difference(
                    item.peak_field_V_cm, reference.peak_field_V_cm
                ),
            }
        )
    return rows


def is_monotone(errors: list[float]) -> bool:
    """误差随网格加密单调下降（容许相等，不容许回升）。"""
    return all(later <= earlier for earlier, later in zip(errors, errors[1:]))
