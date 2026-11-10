"""求解器：Gummel 与 Newton（两者共用 `convergence` 的判据与残差尺度）。

`solve_gummel` 与 `solve_newton` 的签名逐参数对齐（`docs/a_b_interface_week02.md` §3），
输入输出都是物理 `State`，因此 `simulation/dc.py` 只按 `solver.method` 选入口，
扫描、取消、写文件、GUI 的代码对两种方法完全一致——Gate 2 的双求解器对照才有意义。
"""

from minitcad.solvers.convergence import (
    is_converged,
    residual_scales_qf,
    scaled_residual_qf,
    scaled_update,
    scaled_update_physical,
)
from minitcad.solvers.gummel import solve_gummel
from minitcad.solvers.jacobian import jacobian_qf
from minitcad.solvers.newton import solve_newton

__all__ = [
    "solve_gummel",
    "solve_newton",
    "jacobian_qf",
    "scaled_residual_qf",
    "scaled_update",
    "scaled_update_physical",
    "residual_scales_qf",
    "is_converged",
]
