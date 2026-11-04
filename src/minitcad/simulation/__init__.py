"""仿真编排：平衡初值、扫压、继续法与 ``RunResult``。

W3 只有 A 的平衡标量求解初版和 B 的 ``RunResult`` 骨架；``solve`` 在 W4 才接线。
"""

from minitcad.simulation.dc import solve_single_bias
from minitcad.simulation.diagnostics import BiasFields, benchmark_fields, fields_for_attempt
from minitcad.simulation.equilibrium import solve_equilibrium_scalar
from minitcad.simulation.run_result import Attempt, RunResult
from minitcad.simulation.sweep import SweepOutcome, equilibrium_start_state, run_fixed_sweep

__all__ = [
    "Attempt",
    "RunResult",
    "solve_equilibrium_scalar",
    "solve_single_bias",
    "run_fixed_sweep",
    "SweepOutcome",
    "equilibrium_start_state",
    "BiasFields",
    "fields_for_attempt",
    "benchmark_fields",
]
