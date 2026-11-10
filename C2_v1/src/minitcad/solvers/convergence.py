#统一的残差和迭代接口
import numpy as np
from numpy.typing import NDArray
from minitcad.core.contact import contact_targets_qf
from minitcad.core.device import NormalizedDevice
from minitcad.core.state import QFState,State
from minitcad.discretization.residual import ResidualEvaluation
def residual_scales_qf(
        qf : QFState,
        device : NormalizedDevice,
        bias_V : float,
        evaluation : ResidualEvaluation,
        physical : State | None = None
) -> NDArray:
    if physical is None:
        physical = qf.to_physical()
    m = qf.node_count
    if m != device.mesh.node_count or m != physical.node_count:
        raise ValueError("state/device node_count mismatch")
    if (evaluation.f.shape != (3 * m,) or evaluation.recombination_nd.shape != (m,) or
        evaluation.fluxes.jn_edges.shape != (m - 1,) or evaluation.fluxes.jp_edges.shape != (m - 1,)):
        raise ValueError("Invalid QF evaluation shapes")
    f = evaluation.f.reshape(m,3)
    flux = evaluation.fluxes
    r = evaluation.recombination_nd
    w = device.mesh.w_interior_nd
    scale = np.empty((m,3),dtype = np.float64)
    scale[1:-1,0] = np.maximum.reduce((
        np.ones(m - 2),
        np.abs(device.cbar_nd),
        physical.n[1:-1],
        physical.p[1:-1]
    ))
    scale[1:-1,1] = np.maximum.reduce((
        np.ones(m - 2),
        np.abs(flux.jn_edges[:-1]) / w,
        np.abs(flux.jn_edges[1:]) / w,
        np.abs(r[1:-1])
    ))
    scale[1:-1,2] = np.maximum.reduce((
        np.ones(m - 2),
        np.abs(flux.jp_edges[:-1]) / w,
        np.abs(flux.jp_edges[1:]) / w,
        np.abs(r[1:-1])
    ))
    targets = contact_targets_qf(device = device,bias_V = bias_V,ref_n = qf.ref_n,ref_p = qf.ref_p)
    scale[0,:] = np.maximum(
        1.0,
        np.abs(targets[:3])
    )
    scale[-1,:] = np.maximum(
        1.0,
        np.abs(targets[3:])
    )
    if not np.isfinite(scale).all() or np.any(scale <= 0):
        raise FloatingPointError("NONFINITE_RESIDUAL_SCALE")
    return scale.reshape(3 * m)
def scaled_residual_qf(
    qf : QFState,
    device : NormalizedDevice,
    bias_V : float,
    evaluation : ResidualEvaluation,
    physical = None
) -> float:
    scale = residual_scales_qf(qf,device,bias_V,evaluation,physical)
    norm = float(np.max(np.abs(evaluation.f) / scale))
    if not np.isfinite(norm):
        raise FloatingPointError("NONFINITE_SCALED_RESIDUAL")
    return norm
def scaled_update_physical(
    old : State,
    new : State,
    n_floor : float = 1e-15
) -> float:
    if old.node_count != new.node_count:
        raise ValueError("old/new node_count mismatch")
    if not np.isfinite(n_floor) or n_floor <= 0:
        raise ValueError("n_floor must be finite and positive")
    m = old.node_count
    su = np.maximum.reduce((
        np.ones(m),
        np.abs(old.u),
        np.abs(new.u)
    ))
    sn = np.maximum.reduce((
        np.full(m,n_floor),
        old.n,
        new.n
    ))
    sp = np.maximum.reduce((
        np.full(m, n_floor),
        old.p,
        new.p
    ))
    norm = float(max(
        np.max(np.abs(new.u - old.u) / su),
        np.max(np.abs(new.n - old.n) / sn),
        np.max(np.abs(new.p - old.p) / sp)
    ))
    if not np.isfinite(norm):
        raise FloatingPointError("NONFINITE_PHYSICAL_UPDATE")
    return norm