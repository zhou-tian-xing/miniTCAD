#gummel迭代求解
#利用准fermi能级提高了精度,2026/10/21
from __future__ import annotations
import numpy as np
from numpy.typing import NDArray,ArrayLike
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve
from minitcad.core.state import State,SolveResult
from minitcad.core.device import NormalizedDevice
from minitcad import discretization,physics
from minitcad.discretization.sg import edge_flux_qf
from minitcad.discretization.residual import residual_qf
def _state(u:NDArray,fn:NDArray,fp:NDArray) -> State:
    n = np.exp(u + fn)
    p = np.exp(-u - fp)
    return State(u,n,p)
def _quasi_fermi(s:State) -> tuple[NDArray,NDArray]:
    return (np.log(s.n) - s.u,-np.log(s.p) - s.u)
def _solve_tridiagonal(diagonal,left,right,rhs):
    matrix = diags(
        (left[1:],diagonal,right[:-1]),
        offsets = (-1,0,1),
        format = "csc"
    )
    x = np.asarray(spsolve(matrix,rhs),dtype = np.float64)
    if (x.shape != rhs.shape or not np.isfinite(x).all()):
        raise FloatingPointError("LINEAR_SOLVE_FAILED: invalid linear solution")
    return x
def _poisson_step(device:NormalizedDevice,u0,fn,fp,targets,*,max_iter = 35):
    h = device.mesh.h_edges_nd
    w = device.mesh.w_interior_nd
    m = device.mesh.node_count
    u = u0.copy()
    u[0],u[-1] = targets[0],targets[3]
    def evaluate(potential):
        s = _state(potential,fn,fp)
        f = discretization.poisson_residual_nd(device.mesh.x_nodes_nd,
                                               s.u,s.n,s.p,device.cbar_nd,targets[0],targets[3])
        scale = np.maximum.reduce((
            np.ones(m - 2),
            abs(device.cbar_nd),
            s.n[1:-1],
            s.p[1:-1]
        ))
        r = float(np.max(abs(f[1:-1]) / scale))
        return s,f,r
    for i in range(max_iter):
        s,f,rn = evaluate(u)
        if rn <= 1e-10:
            return u
        diag = (
            -(1.0 / h[:-1] + 1.0 / h[1:]) / w
            - s.n[1:-1]
            - s.p[1:-1]
        )#M - 2,内部节点数
        left,right = np.empty(m - 2,dtype = np.float64),np.empty(m - 2,dtype = np.float64)
        left[1:] = 1.0 / (w[1:] * h[1:-1])#M - 3
        right[:-1] = 1.0 / (w[:-1] * h[1:-1])
        step = _solve_tridiagonal(diag,left,right,-f[1:-1])
        alpha = 1.0
        accepted = False
        while alpha >= 1e-6:
            trial = u.copy()
            trial[1:-1] += alpha * step
            try:
                _,_,trial_rn = evaluate(trial)
                if trial_rn < rn or trial_rn <= 1e-10:
                    u = trial
                    accepted = True
                    break
            except (ValueError, FloatingPointError, OverflowError):
                pass
            alpha *= 0.5
        if not accepted:
            raise ArithmeticError("POISSON_LINE_SEARCH_FAILED")
    raise ArithmeticError("POISSON_MAX_INNER")
def _electron_step(device:NormalizedDevice,u,n_lag,p_fixed,targets):
    h = device.mesh.h_edges_nd
    w = device.mesh.w_interior_nd
    du = np.diff(u)
    denominator = device.tau_n_nd * (p_fixed[1:-1] + 1) + device.tau_p_nd * (n_lag[1:-1] + 1)
    bpos = discretization.bernoulli(du)
    bneg = discretization.bernoulli(-du)
    right = - (device.mu_n_nd * bpos[1:]) / (w * h[1:])
    left = - (device.mu_n_nd * bneg[:-1]) / (w * h[:-1])
    diagonal = (device.mu_n_nd * (bneg[1:] / h[1:] + bpos[:-1] / h[:-1]) / w
                + p_fixed[1:-1] / denominator)
    rhs = 1.0 / denominator
    rhs[0] -= left[0] * targets[1]
    rhs[-1] -= right[-1] * targets[4]
    ni = _solve_tridiagonal(diagonal = diagonal,left = left,
                            right = right,rhs = rhs)
    if np.any(ni <= 0):
        raise ArithmeticError("NONPOSITIVE_ELECTRON_SUBSTEP")
    n = np.empty_like(u)
    n[0] = targets[1]
    n[-1] = targets[4]
    n[1:-1] = ni
    return n
def _hole_step(device:NormalizedDevice,u,p_lag,n_fixed,targets):
    h = device.mesh.h_edges_nd
    w = device.mesh.w_interior_nd
    du = np.diff(u)
    denominator = device.tau_n_nd * (p_lag[1:-1] + 1) + device.tau_p_nd * (n_fixed[1:-1] + 1)
    bpos = discretization.bernoulli(du)
    bneg = discretization.bernoulli(-du)
    right = - (device.mu_p_nd * bneg[1:]) / (w * h[1:])
    left = - (device.mu_p_nd * bpos[:-1]) / (w * h[:-1])
    diagonal = (device.mu_p_nd * (bpos[1:] / h[1:] + bneg[:-1] / h[:-1]) / w
                + n_fixed[1:-1] / denominator)
    rhs = 1.0 / denominator
    rhs[0] -= left[0] * targets[2]
    rhs[-1] -= right[-1] * targets[5]
    pi = _solve_tridiagonal(diagonal = diagonal,left = left,
                            right = right,rhs = rhs)
    if np.any(pi <= 0):
        raise ArithmeticError("NONPOSITIVE_HOLE_SUBSTEP")
    p = np.empty_like(u)
    p[0] = targets[2]
    p[-1] = targets[5]
    p[1:-1] = pi
    return p
def _scaled_residual(
    device:NormalizedDevice,
    bias_V:float,
    s:State
) -> float:
    m = s.node_count
    F = discretization.residual(s,device,bias_V).reshape((m,3))
    flux = discretization.edge_flux(s,device)
    r = physics.srh_recombination_nd(
        s.n,
        s.p,
        device.tau_n_nd,
        device.tau_p_nd
    )
    w = device.mesh.w_interior_nd
    t = device.contact_targets(bias_V)
    scale = np.ones((m,3),dtype = float)
    scale[1:-1,0] = np.maximum.reduce((
        np.ones(m - 2),
        abs(device.cbar_nd),
        s.n[1:-1],
        s.p[1:-1]
    ))
    scale[1:-1,1] = np.maximum.reduce((
        np.ones(m - 2),
        abs(flux.jn_edges[:-1]) / w,
        abs(flux.jn_edges[1:]) / w,
        abs(r[1:-1])
    ))
    scale[1:-1, 2] = np.maximum.reduce((
        np.ones(m - 2),
        abs(flux.jp_edges[:-1]) / w,
        abs(flux.jp_edges[1:]) / w,
        abs(r[1:-1])
    ))
    scale[0] = [max(1, abs(a)) for a in t[:3]]
    scale[-1] = [max(1, abs(a)) for a in t[3:]]
    return float(np.max(abs(F) / scale))
def _scaled_update(
    old:State,
    new:State,
    floor:float,
) -> float:
    su = np.maximum.reduce((
        np.ones(old.node_count),
        abs(old.u),
        abs(new.u)
    ))
    sn = np.maximum.reduce((
        np.full(old.node_count, floor),
        old.n,
        new.n
    ))
    sp = np.maximum.reduce((
        np.full(old.node_count, floor),
        old.p,
        new.p
    ))
    return float(max(
        np.max(abs(new.u - old.u) / su),
        np.max(abs(new.n - old.n) / sn),
        np.max(abs(new.p - old.p) / sp)
    ))
#费米能级部分
def _scaled_residual_qf(
        device:NormalizedDevice,
        bias_V,
        u,
        fn,
        fp,*,
        ref_n = 0.0,
        ref_p = 0.0
):
    m = u.size
    f,flux,rate = residual_qf(device,bias_V,u,fn,fp,ref_n = ref_n,ref_p = ref_p)
    n = np.exp(u + ref_n + fn)
    p = np.exp(-u - ref_p - fp)
    w = device.mesh.w_interior_nd
    t = device.contact_targets(bias_V = bias_V)
    scale = np.empty((m,3))
    scale[1:-1,0] = np.maximum.reduce((
        np.ones(m - 2),
        abs(device.cbar_nd),
        n[1:-1],
        p[1:-1]
    ))
    scale[1:-1,1] = np.maximum.reduce((
        np.ones(m - 2),
        abs(flux.jn_edges[:-1]) / w,
        abs(flux.jn_edges[1:]) / w,
        abs(rate[1:-1])
    ))
    scale[1:-1,2] = np.maximum.reduce((
        np.ones(m - 2),
        abs(flux.jp_edges[:-1]) / w,
        abs(flux.jp_edges[1:]) / w,
        abs(rate[1:-1])
    ))
    scale[0, : ] = np.maximum(
        1.0,
        abs(np.asarray(t[:3]))
    )
    scale[-1, : ] = np.maximum(
        1.0,
        abs(np.asarray(t[3:]))
    )
    return float(np.max(abs(f) / scale))
def _qf_carrier_step(
    device:NormalizedDevice,
    u,
    eta_n,
    eta_p,
    *,
    electron,
    ref_n = 0.0,
    ref_p = 0.0,
    max_iter = 30,
    tol = 1e-11
):
    h = device.mesh.h_edges_nd
    w = device.mesh.w_interior_nd
    du = np.diff(u)
    bp = discretization.bernoulli(du)
    bn = discretization.bernoulli(-du)
    fn = eta_n.copy()
    fp = eta_p.copy()
    nlag = np.exp((u + ref_n) + fn)
    plag = np.exp((-u - ref_p) - fp)
    denominator = (
        device.tau_n_nd * (plag[1:-1] + 1.0) + device.tau_p_nd * (nlag[1:-1] + 1.0)
    )
    def eval_at(eta):
        if electron:
            fn2 = eta
            fp2 = fp
        else:
            fn2 = fn
            fp2 = eta
        n2 = np.exp((u + ref_n) + fn2)
        p2 = np.exp((-u - ref_p) - fp2)
        flux = edge_flux_qf(
            u,
            fn2,
            fp2,
            h,
            device.mu_n_nd,
            device.mu_p_nd,
            ref_n = ref_n,
            ref_p = ref_p,
        )
        r = np.expm1(
            (ref_n - ref_p)
            + (fn2[1:-1] - fp2[1:-1])
        ) / denominator
        if electron:
            j = flux.jn_edges
            residual = np.diff(j) / w - r
            m_n = device.mu_n_nd
            left = - (
                m_n * bn[:-1] * n2[:-2]
                / (w * h[:-1])
            )
            right = - (
                m_n * bp[1:] * n2[2:]
                / (w * h[1:])
            )
            diag = (
                m_n * n2[1:-1]
                * (
                    bn[1:] / h[1:]
                    + bp[:-1] / h[:-1]
                ) / w
                + n2[1:-1] * p2[1:-1] / denominator
            )
        else:
            j = flux.jp_edges
            residual = np.diff(j) / w + r
            m_p = device.mu_p_nd
            left = - (
                m_p * bp[:-1] * p2[:-2]
                / (w * h[:-1])
            )
            right = - (
                m_p * bn[1:] * p2[2:]
                / (w * h[1:])
            )
            diag = (
                m_p * p2[1:-1]
                * (
                    bp[1:] / h[1:]
                    + bn[:-1] / h[:-1]
                ) / w
                + n2[1:-1] * p2[1:-1] / denominator
            )
        scale = np.maximum.reduce((
            np.ones_like(residual),
            abs(j[1:]) / w,
            abs(j[:-1]) / w,
            abs(r)
        ))
        rn = float(
            np.max(abs(residual) / scale)
        )
        return residual,diag,left,right,rn
    current = fn if electron else fp
    for inner in range(max_iter):
        f, diag, left, right, rn = eval_at(current)
        if rn < tol:
            return current
        delta = _solve_tridiagonal(
            diag,
            left,
            right,
            f,
        )
        alpha = 1.0
        good = False
        for _ in range(30):
            trial = current.copy()
            trial[1:-1] += alpha * delta
            try:
                *_,newrn = eval_at(trial)
                if newrn < rn or newrn < tol:
                    current = trial
                    good = True
                    break
            except:
                pass
            alpha *= 0.5
        if not good:
            name = "ELECTRON" if electron else "HOLE"
            raise ArithmeticError(
                f"{name}_INNER_LINE_SEARCH_FAILED "
                f"at {rn:.3g}"
            )
    name = "ELECTRON" if electron else "HOLE"
    raise ArithmeticError(
        f"{name}_INNER_MAX_ITER rn={rn:.3g}"
    )
"""
def solve_gummel(
    device:NormalizedDevice,
    bias_V:float,
    initial_state:State,
    *,
    tol_residual = 1e-8,
    tol_update = 1e-6,
    max_iter = 60,
    alpha_min = 1e-6,
    n_floor = 1e-15
) -> SolveResult:
    if initial_state.node_count != device.mesh.node_count:
        raise ValueError("device/state node_count mismatch")
    for name,x in [
        ("bias_V",bias_V),
        ("tol_residual",tol_residual),
        ("tol_update",tol_update),
        ("alpha_min",alpha_min),
        ("n_floor",n_floor)
    ]:
        if not np.isfinite(x).all() or (name != "bias_V" and x <= 0):
            raise ValueError(f"invaild {name}")
    if (alpha_min > 1):
        raise ValueError("alpha_min must be <= 1")
    if (
        isinstance(max_iter,bool)
        or not isinstance(max_iter,(int,np.integer))
        or max_iter < 1
    ):
        raise ValueError("max_iter must be positive integer")
    t = device.contact_targets(bias_V = bias_V)
    u = initial_state.u.copy()
    n = initial_state.n.copy()
    p = initial_state.p.copy()
    u[0],u[-1] = t[0],t[3]
    n[0],n[-1] = t[1],t[4]
    p[0],p[-1] = t[2],t[5]
    old = State(u,n,p)
    history : list[dict] = []
    try:
        rn = _scaled_residual(device,bias_V,old)
    except:
        return SolveResult(
            bias_V,None,False,0,np.inf,np.inf,"INVALID_INITIAL_STATE",tuple(history)
        )
    if rn <= tol_residual:
        return SolveResult(
            bias_V,old,True,0,rn,0.0,"SUCCESS",tuple(history)
        )
    dn = np.inf
    for it in range(1,max_iter + 1):
        try:
            (fn,fp) = _quasi_fermi(old)
            ustar = _poisson_step(
                device,old.u,fn,fp,t
            )
            intermediate = _state(ustar,fn,fp)
            nstar = _electron_step(
                device,
                ustar,
                intermediate.n,
                intermediate.p,
                t
            )
            pstar = _hole_step(
                device,
                ustar,
                intermediate.p,
                nstar,
                t
            )
            fnstar = np.log(nstar) - ustar
            fpstar = -np.log(pstar) - ustar
        except Exception as e:
            history.append({
                "iteration": it,
                "status": "SUBSTEP_FAILED",
                "reason": str(e)
            })
            return SolveResult(
                bias_V,None,False,it,
                rn,dn,
                "SUBSTEP_FAILED: " + str(e),
                tuple(history)
            )
        alpha = 1.0
        accepted = False
        while alpha >= alpha_min:
            try:
                trial = _state(
                    old.u + alpha * (ustar - old.u),
                    fn + alpha * (fnstar - fn),
                    fp + alpha * (fpstar - fp)
                )
                trial_rn = _scaled_residual(
                    device,bias_V,trial
                )
                trial_dn = _scaled_update(
                    old,trial,n_floor
                )
                good = (
                    np.isfinite(trial_rn)
                    and (
                        trial_rn < rn
                        or trial_rn <= tol_residual
                    )
                )
            except:
                good = False
                trial_rn = np.inf
                trial_dn = np.inf
            history.append({
                "iteration": it,
                "alpha": float(alpha),
                "residual_norm": float(trial_rn),
                "update_norm": float(trial_dn),
                "status": (
                    "ACCEPTED" if good else "REJECTED"
                )
            })
            if good:
                old = trial
                rn = trial_rn
                dn = trial_dn
                accepted = True
                break
            alpha *= 0.5
        if not accepted:
            return SolveResult(
                bias_V,None,False,it,
                rn,dn,
                "LINE_SEARCH_FAILED",
                tuple(history)
            )
        if rn <= tol_residual and dn <= tol_update:
            return SolveResult(
                bias_V,old,True,it,
                rn,dn,
                "SUCCESS",
                tuple(history)
            )
    return SolveResult(
        bias_V,None,False,max_iter,
        rn,dn,
        "MAX_ITER",
        tuple(history)
    )
"""
#以下是vibecoding
def solve_gummel(
    device: NormalizedDevice,
    bias_V: float,
    initial_state: State,
    *,
    tol_residual=1e-8,
    tol_update=1e-6,
    max_iter=60,
    alpha_min=1e-6,
    n_floor=1e-15,
    diagnostics=None,
) -> SolveResult:
    if initial_state.node_count != device.mesh.node_count:
        raise ValueError(
            "device/state node_count mismatch"
        )
    for name, x in [
        ("bias_V", bias_V),
        ("tol_residual", tol_residual),
        ("tol_update", tol_update),
        ("alpha_min", alpha_min),
        ("n_floor", n_floor),
    ]:
        if (
            not np.isfinite(x).all()
            or (name != "bias_V" and x <= 0)
        ):
            raise ValueError(f"invalid {name}")
    if alpha_min > 1:
        raise ValueError("alpha_min must be <= 1")
    if (
        isinstance(max_iter, bool)
        or not isinstance(max_iter, (int, np.integer))
        or max_iter < 1
    ):
        raise ValueError(
            "max_iter must be positive integer"
        )
    if diagnostics is not None:
        diagnostics.clear()
    t = device.contact_targets(bias_V)
    u = initial_state.u.copy()
    fn, fp = _quasi_fermi(initial_state)
    ref_n = np.log(t[4]) - t[3]
    ref_p = -np.log(t[2]) - t[0]
    fn -= ref_n
    fp -= ref_p
    u[0], u[-1] = t[0], t[3]
    fn[0] = np.log(t[1]) - t[0] - ref_n
    fn[-1] = 0.0
    fp[0] = 0.0
    fp[-1] = -np.log(t[5]) - t[3] - ref_p
    history = []
    def save_diagnostics():
        if diagnostics is not None:
            ff, currents, rate = residual_qf(
                device,
                bias_V,
                u,
                fn,
                fp,
                ref_n=ref_n,
                ref_p=ref_p,
            )
            diagnostics.update({
                "u": u.copy(),
                "eta_n": fn.copy(),
                "eta_p": fp.copy(),
                "ref_n": float(ref_n),
                "ref_p": float(ref_p),
                "jn_edges": currents.jn_edges.copy(),
                "jp_edges": currents.jp_edges.copy(),
                "rate": rate.copy(),
                "residual": ff.copy(),
            })
    try:
        old = _state(
            u,
            fn + ref_n,
            fp + ref_p,
        )
        rn = _scaled_residual_qf(
            device,
            bias_V,
            u,
            fn,
            fp,
            ref_n=ref_n,
            ref_p=ref_p,
        )
    except Exception:
        return SolveResult(
            bias_V,
            None,
            False,
            0,
            np.inf,
            np.inf,
            "INVALID_INITIAL_STATE",
            tuple(history),
        )
    if rn <= tol_residual:
        save_diagnostics()
        return SolveResult(
            bias_V,
            old,
            True,
            0,
            rn,
            0.0,
            "SUCCESS",
            tuple(history),
        )
    dn = np.inf
    for it in range(1, max_iter + 1):
        try:
            ustar = _poisson_step(
                device,
                u,
                fn + ref_n,
                fp + ref_p,
                t,
            )
            fnstar = _qf_carrier_step(
                device,
                ustar,
                fn,
                fp,
                electron=True,
                ref_n=ref_n,
                ref_p=ref_p,
            )
            fpstar = _qf_carrier_step(
                device,
                ustar,
                fnstar,
                fp,
                electron=False,
                ref_n=ref_n,
                ref_p=ref_p,
            )
        except Exception as e:
            history.append({
                "iteration": it,
                "status": "SUBSTEP_FAILED",
                "reason": str(e),
            })
            return SolveResult(
                bias_V,
                None,
                False,
                it,
                rn,
                dn,
                "SUBSTEP_FAILED: " + str(e),
                tuple(history),
            )
        alpha = 1.0
        accepted = False
        while alpha >= alpha_min:
            try:
                tu = u + alpha * (ustar - u)
                tf = fn + alpha * (fnstar - fn)
                tp = fp + alpha * (fpstar - fp)
                trial = _state(
                    tu,
                    tf + ref_n,
                    tp + ref_p,
                )
                trial_rn = _scaled_residual_qf(
                    device,
                    bias_V,
                    tu,
                    tf,
                    tp,
                    ref_n=ref_n,
                    ref_p=ref_p,
                )
                trial_dn = _scaled_update(
                    old,
                    trial,
                    n_floor,
                )
                good = (
                    np.isfinite(trial_rn)
                    and (
                        trial_rn < rn
                        or trial_rn <= tol_residual
                    )
                )
            except Exception:
                good = False
                trial_rn = np.inf
                trial_dn = np.inf
            history.append({
                "iteration": it,
                "alpha": float(alpha),
                "residual_norm": float(trial_rn),
                "update_norm": float(trial_dn),
                "status": (
                    "ACCEPTED"
                    if good
                    else "REJECTED"
                ),
            })
            if good:
                old = trial
                u = tu
                fn = tf
                fp = tp
                rn = trial_rn
                dn = trial_dn
                accepted = True
                break
            alpha *= 0.5
        if not accepted:
            return SolveResult(
                bias_V,
                None,
                False,
                it,
                rn,
                dn,
                "LINE_SEARCH_FAILED",
                tuple(history),
            )
        if (
            rn <= tol_residual
            and dn <= tol_update
        ):
            save_diagnostics()
            return SolveResult(
                bias_V,
                old,
                True,
                it,
                rn,
                dn,
                "SUCCESS",
                tuple(history),
            )
    return SolveResult(
        bias_V,
        None,
        False,
        max_iter,
        rn,
        dn,
        "MAX_ITER",
        tuple(history),
    )