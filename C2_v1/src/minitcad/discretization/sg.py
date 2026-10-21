#利用Scharfetter–Gummel离散求电流
#修订于2026/10/19：采用了更为稳健的数值算法
#利用准fermi能级提高了精度,2026/10/21
from minitcad.discretization.bernoulli import bernoulli
from minitcad.core.state import State,Fluxes
from minitcad.core.device import NormalizedDevice
import numpy as np
from numpy.typing import ArrayLike,NDArray
def _checked_input(name:str,value:ArrayLike) -> NDArray[np.float64]:
    if np.iscomplexobj(value):
        raise ValueError(f"{name} must be real-valued")
    x = np.asarray(value)
    if (x.ndim > 1):
        raise ValueError(f"{name} must be a scalar or a 1D array")
    if not np.isfinite(x).all():
        raise ValueError(f"{name} must contain finite values")
    return x
def _edge_inputs(
        left:ArrayLike,
        right:ArrayLike,
        delta_u:ArrayLike,
        h:ArrayLike,
        mobility:ArrayLike
) -> tuple[NDArray[np.float64],...]:
    names = ("left density", "right density", "delta_u", "h", "mobility")
    values = (left,right,delta_u,h,mobility)
    arrays = tuple(
        _checked_input(name,value)
        for name,value in zip(names,values)
    )
    try:
        l,r,du,spacing,m = np.broadcast_arrays(*arrays)
    except ValueError as e:
        raise ValueError("SG edge inputs have incompatible shapes") from e
    for name,a in (("left density",l),("right density",r),
                    ("h",spacing), ("mobility",m)):
        if np.any(a <= 0.0):
            raise ValueError(f"{name} must be strictly positive")
    return (l, r, du, spacing, m)
def _stable_sg(
        left:ArrayLike,
        right:ArrayLike,
        du:ArrayLike,
        h:ArrayLike,
        mobility:ArrayLike,
        electron:bool
):
    d_log = np.log(right) - np.log(left)
    if electron:
        qf_step = d_log - du
    else:
        qf_step = d_log + du
    result = np.empty_like(qf_step,dtype = np.float64)
    positive = qf_step > 0
    nonpositive = ~ positive
    if electron:
        if np.any(positive):
            result[positive] = - (mobility[positive] / h[positive] * right[positive]
                                  * bernoulli(du[positive]) * np.expm1(-qf_step[positive]))
        if np.any(nonpositive):
            result[nonpositive] = (mobility[nonpositive] / h[nonpositive] * left[nonpositive]
                                   * bernoulli(-du[nonpositive]) * np.expm1(qf_step[nonpositive]))
    else:
        if np.any(nonpositive):
            result[nonpositive] = - (mobility[nonpositive] / h[nonpositive] * left[nonpositive]
                                      * bernoulli(du[nonpositive]) * np.expm1(qf_step[nonpositive]))
        if np.any(positive):
            result[positive] = (mobility[positive] / h[positive] * right[positive]
                                 * bernoulli(-du[positive]) * np.expm1(-qf_step[positive]))
    if not np.isfinite(result).all():
        raise ValueError("SG edge current is not finite")
    if result.ndim == 0:
        return float(result)
    else:
        return result
def electron_edge_flux(
        n_left:ArrayLike,
        n_right:ArrayLike,
        delta_u:ArrayLike,
        h:ArrayLike,
        m_n:ArrayLike
) -> float | NDArray[np.float64]:
    nl,nr,du,spacing,m = _edge_inputs(n_left,n_right,delta_u,h,m_n)
    jn = np.empty_like(du,dtype = np.float64)
    zero = (du == 0)
    if np.any(zero):
        jn[zero] = (m[zero] / spacing[zero]) * (nr[zero] - nl[zero])
    if np.any(~zero):
        select = ~ zero
        jn[select] = _stable_sg(nl[select],nr[select],du[select],spacing[select],m[select],True)
    #jn = (m / spacing) * (nr * bernoulli(du) - nl * bernoulli(-du))这个式子在高精度下存在关键问题
    if not np.isfinite(jn).all():
        raise ValueError("electron edge flux is not finite")
    if jn.ndim == 0:
        return float(jn)
    else:
        return np.asarray(jn,dtype = np.float64)
def hole_edge_flux(
        p_left:ArrayLike,
        p_right:ArrayLike,
        delta_u:ArrayLike,
        h:ArrayLike,
        m_p:ArrayLike
) -> float | NDArray[np.float64]:
    pl,pr,du,spacing,m = _edge_inputs(p_left,p_right,delta_u,h,m_p)
    jp = np.empty_like(du,dtype = np.float64)
    zero = (du == 0)
    if np.any(zero):
        jp[zero] = (m[zero] / spacing[zero]) * (pl[zero] - pr[zero])
    if np.any(~zero):
        select = ~ zero
        jp[select] = _stable_sg(pl[select],pr[select],du[select],spacing[select],m[select],False)
    #jp = (m / spacing) * (pl * bernoulli(du) - pr * bernoulli(-du))
    if not np.isfinite(jp).all():
        raise ValueError("hole edge flux is not finite")
    if jp.ndim == 0:
        return float(jp)
    else:
        return np.asarray(jp,dtype = np.float64)
def edge_flux(state:State,device:NormalizedDevice):
    u = _checked_input("state.u",state.u)
    n = _checked_input("state.n",state.n)
    p = _checked_input("state.p",state.p)
    if (u.ndim != 1 or u.size < 3 or
         n.shape != u.shape or p.shape != u.shape):
        raise ValueError("state.u, state.n, state.p must have the same 1D shape (M,), M >= 3")
    h = _checked_input("device.mesh.h_edges_nd",device.mesh.h_edges_nd)
    if h.shape != (u.size - 1,):
        raise ValueError("device.mesh.h_edges_nd must have shape (M-1,)")
    du = np.diff(u)
    jn = electron_edge_flux(n[:-1],n[1:],du,h,device.mu_n_nd)
    jp = hole_edge_flux(p[:-1],p[1:],du,h,device.mu_p_nd)
    return Fluxes(
        jn_edges = np.asarray(jn),
        jp_edges = np.asarray(jp)
    )
def edge_flux_qf(
        u,
        eta_n,
        eta_p,
        h,
        mu_n,
        mu_p,
        *,
        ref_n = 0.0,
        ref_p = 0.0
):
    #fn = ref_n + eta_n
    #fp = ref_p + eta_p
    u = np.asarray(u,dtype = np.float64)
    eta_n = np.asarray(eta_n,dtype = np.float64)
    eta_p = np.asarray(eta_p,dtype = np.float64)
    h = np.asarray(h,dtype = np.float64)
    if (u.ndim != 1 or u.size < 3 or 
        eta_n.shape != u.shape or eta_p.shape != u.shape or h.shape != (u.size - 1,)):
        raise ValueError("invalid quasi-Fermi flux shapes")
    if (not all(np.isfinite(a).all() for a in (u,eta_n,eta_p,h))
        or np.any(h <= 0)
        or not np.isfinite([mu_n,mu_p,ref_n,ref_p]).all()
        or mu_n <= 0 or mu_p <= 0):
        raise ValueError("invalid quasi-Fermi flux inputs")
    n = np.exp((u + ref_n) + eta_n)
    p = np.exp((-u - ref_p) - eta_p)
    du = np.diff(u)
    dfn = np.diff(eta_n)
    dfp = np.diff(eta_p)
    bp = bernoulli(du)
    bn = bernoulli(-du)
    jn = np.empty_like(du)
    jp = np.empty_like(du)
    pos_n = (dfn > 0.0)
    jn[pos_n] = (
        -(mu_n / h[pos_n]) * n[1:][pos_n] * bp[pos_n] * np.expm1(-dfn[pos_n])
    )
    jn[~pos_n] = (
        (mu_n / h[~pos_n]) * n[:-1][~pos_n] * bn[~pos_n] * np.expm1(dfn[~pos_n])
    )
    pos_p =(dfp > 0.0)
    jp[pos_p] = (
        -(mu_p / h[pos_p]) * p[:-1][pos_p] * bp[pos_p] * np.expm1(-dfp[pos_p])
    )
    jp[~pos_p] = (
        (mu_p / h[~pos_p]) * p[1:][~pos_p] * bn[~pos_p] * np.expm1(dfp[~pos_p])
    )
    return Fluxes(
        jn_edges = jn,
        jp_edges = jp
    )