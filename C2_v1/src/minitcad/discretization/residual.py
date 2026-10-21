#组装discretization的全部代码
#利用准fermi能级提高了精度,2026/10/21
import numpy as np
from numpy.typing import NDArray
from minitcad.core.contact import contact_targets
from minitcad.core.device import NormalizedDevice
from minitcad.core.state import State
from minitcad.discretization.continuity import continuity_residual_nd
from minitcad.discretization.poisson import poisson_residual_nd
from minitcad.discretization.sg import edge_flux,edge_flux_qf
from minitcad.physics.recombination import srh_recombination_nd,srh_recombination_qf
def residual(state:State,device:NormalizedDevice,bias_V:float) -> NDArray[np.float64]:
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
def residual_qf(
    device,
    bias_V,
    u,
    eta_n,
    eta_p,
    *,
    ref_n=0.0,
    ref_p=0.0
):
    u = np.asarray(u, dtype=np.float64)
    eta_n = np.asarray(eta_n,dtype=np.float64)
    eta_p = np.asarray(eta_p,dtype=np.float64)
    m = device.mesh.node_count
    if (
        u.shape != (m,)
        or eta_n.shape != (m,)
        or eta_p.shape != (m,)
    ):
        raise ValueError(
            "state/device quasi-Fermi shape mismatch"
        )
    n = np.exp((u + ref_n) + eta_n)
    p = np.exp((-u - ref_p) - eta_p)
    targets = device.contact_targets(bias_V)
    fu = poisson_residual_nd(
        device.mesh.x_nodes_nd,
        u,
        n,
        p,
        device.cbar_nd,
        targets[0],
        targets[3]
    )
    flux = edge_flux_qf(
        u,
        eta_n,
        eta_p,
        device.mesh.h_edges_nd,
        device.mu_n_nd,
        device.mu_p_nd,
        ref_n = ref_n,
        ref_p = ref_p
    )
    rate = srh_recombination_qf(
        u,
        eta_n,
        eta_p,
        device.tau_n_nd,
        device.tau_p_nd,
        ref_n = ref_n,
        ref_p = ref_p
    )
    fe, fh = continuity_residual_nd(
        flux,
        rate,
        device.mesh.w_interior_nd
    )
    f = np.empty((m,3),dtype = np.float64)
    f[:,0] = fu
    f[1:-1,1] = fe
    f[1:-1,2] = fh
    f[0,1] = n[0] - targets[1]
    f[-1,1] = n[-1] - targets[4]
    f[0,2] = p[0] - targets[2]
    f[-1,2] = p[-1] - targets[5]
    if not np.isfinite(f).all():
        raise FloatingPointError(
            "nonfinite quasi-Fermi residual"
        )
    return f,flux,rate