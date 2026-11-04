"""求解器：Gummel（W4）与收敛判据；Newton 在 W5 由 A 补齐。

`solve_gummel` 是 B 的 `run()` 在当前阶段唯一可用的数值入口，签名与
`docs/a_b_interface_week02.md` §3 冻结的口径一致：输入输出都是物理 `State`。
"""

from minitcad.solvers.convergence import (
    is_converged,
    scaled_residual_qf,
    scaled_update,
)
from minitcad.solvers.gummel import solve_gummel

__all__ = ["solve_gummel", "scaled_residual_qf", "scaled_update", "is_converged"]
