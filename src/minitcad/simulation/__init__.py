"""仿真编排：平衡初值、扫压、继续法与 ``RunResult``。

W3 只有 A 的平衡标量求解初版和 B 的 ``RunResult`` 骨架；``solve`` 在 W4 才接线。
"""

from minitcad.simulation.equilibrium import solve_equilibrium_scalar
from minitcad.simulation.run_result import Attempt, RunResult

__all__ = ["Attempt", "RunResult", "solve_equilibrium_scalar"]
