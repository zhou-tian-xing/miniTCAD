"""离散模块（A 的 W2/W3 交付）：Bernoulli、Poisson、SG、continuity 与唯一残差。

生产入口是 ``evaluate_qf`` / ``residual_qf``（QF 口径）；``residual`` 只保留为历史/回归接口。
"""

from minitcad.discretization.bernoulli import bernoulli
from minitcad.discretization.poisson import poisson_residual_nd, poisson_discretization
from minitcad.discretization.sg import (
    edge_flux,
    edge_flux_qf,
    electron_edge_flux,
    hole_edge_flux,
)
from minitcad.discretization.continuity import continuity_residual_nd
from minitcad.discretization.residual import (
    ResidualEvaluation,
    evaluate_qf,
    residual,
    residual_qf,
)

__all__ = [
    "bernoulli",
    "poisson_residual_nd",
    "poisson_discretization",
    "edge_flux",
    "edge_flux_qf",
    "electron_edge_flux",
    "hole_edge_flux",
    "continuity_residual_nd",
    "residual",
    "residual_qf",
    "evaluate_qf",
    "ResidualEvaluation",
]
