#组装discretization的全部代码
#利用准fermi能级提高了精度,2026/10/21，以下所有文件，包括求解器核心部分的接口，以qf为准。原residual只是为了保留历史测试接口
#以下所有文件，包括求解器核心部分的接口，以qf为准
import numpy as np
from numpy.typing import NDArray
from minitcad.core.contact import contact_targets,contact_targets_qf
from minitcad.core.device import NormalizedDevice
from minitcad.core.state import State,Fluxes,QFState
from minitcad.discretization.continuity import continuity_residual_nd
from minitcad.discretization.poisson import poisson_residual_nd
from minitcad.discretization.sg import edge_flux,edge_flux_qf
from minitcad.physics.recombination import srh_recombination_nd,srh_recombination_qf
from dataclasses import dataclass
@dataclass(frozen = True,slots = True)
class ResidualEvaluation:
    f : NDArray[np.float64]
    fluxes : Fluxes
    recombination_nd : NDArray[np.float64]
def residual(state:State,device:NormalizedDevice,bias_V:float) -> NDArray[np.float64]:
    #仅作为历史遗留接口，不承担实际作用
    mesh = device.mesh
    m = mesh.node_count
    if (m != state.node_count):
        raise ValueError("State and NormalizedDevice must have the same node_count")
    targets = contact_targets(device = device,bias_V = bias_V)
    fu = poisson_residual_nd(mesh.x_nodes_nd,state.u,state.n,state.p,
                             device.cbar_nd,targets[0],targets[3])
    fluxes = edge_flux(state = state,device = device)
    rate = srh_recombination_nd(state.n,state.p,device.tau_n_nd,device.tau_p_nd)
    (fn,fp) = continuity_residual_nd(fluxes = fluxes,recombination_nd = rate,
                                     w_interior_nd = mesh.w_interior_nd)
    f = np.empty((m,3),dtype = np.float64)
    f[:,0] = fu
    f[1:-1,1] = fn
    f[1:-1,2] = fp
    f[0,1] = state.n[0] - targets[1]
    f[-1,1] = state.n[-1] - targets[4]
    f[0,2] = state.p[0] - targets[2]
    f[-1,2] = state.p[-1] - targets[5]
    if not np.isfinite(f).all():
        raise FloatingPointError("couples residual is not finite")
    return f.reshape((3*m,))
def evaluate_qf(
    state : QFState,
    device : NormalizedDevice,
    bias_V : float
) -> ResidualEvaluation:
    m = device.mesh.node_count
    if state.node_count != m:
        raise ValueError("QFState/device node_count mismatch")
    state_physical = state.to_physical()
    targets = contact_targets_qf(device,bias_V,ref_n = state.ref_n,ref_p = state.ref_p)
    fu = poisson_residual_nd(
        device.mesh.x_nodes_nd,
        state.u,
        state_physical.n,
        state_physical.p,
        device.cbar_nd,
        targets[0],
        targets[3]
    )
    flux = edge_flux_qf(
        state.u,
        state.eta_n,
        state.eta_p,
        device.mesh.h_edges_nd,
        device.mu_n_nd,
        device.mu_p_nd,
        ref_n = state.ref_n,
        ref_p = state.ref_p
    )
    rate = srh_recombination_qf(
        state.u,
        state.eta_n,
        state.eta_p,
        device.tau_n_nd,
        device.tau_p_nd,
        ref_n = state.ref_n,
        ref_p = state.ref_p
    )
    fe,fh = continuity_residual_nd(
        flux,
        rate,
        device.mesh.w_interior_nd
    )
    f = np.empty((m,3),dtype = np.float64)
    f[:,0] = fu
    f[1:-1,1] = fe
    f[1:-1,2] = fh
    f[0,1] = state.eta_n[0] - targets[1]
    f[-1,1] = state.eta_n[-1] - targets[4]
    f[0,2] = state.eta_p[0] - targets[2]
    f[-1,2] = state.eta_p[-1] - targets[5]
    if not np.isfinite(f).all():
        raise FloatingPointError(
            "nonfinite quasi-Fermi residual"
        )
    return ResidualEvaluation(
        f = f.reshape((3 * m,)),
        fluxes = flux,
        recombination_nd = rate
    )
def residual_qf(
        state : QFState,
        device :NormalizedDevice,
        bias_V : float
) -> NDArray[np.float64]:
    return evaluate_qf(state,device,bias_V).f