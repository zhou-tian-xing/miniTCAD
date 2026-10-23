from minitcad.core.contact import contact_targets_qf
from minitcad.core.device import NormalizedDevice
from minitcad.core.state import QFState,State
from minitcad.discretization.residual import ResidualEvaluation
def scaled_residual_qf(
        qf : QFState,
        device : NormalizedDevice,
        bias_V : float,
        evaluation : ResidualEvaluation,
        physical : State | None = None
) -> float:
    if physical is None:
        physical = qf.to_physical()
    m = qf.node_count