#gummel迭代求解
#利用准fermi能级提高了精度,2026/10/21
from __future__ import annotations
import numpy as np
from numpy.typing import NDArray,ArrayLike
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve
from minitcad.core.state import State,SolveResult,QFState,PhysicsDiagnostics
from minitcad.core.contact import contact_targets_qf
from minitcad.core.device import NormalizedDevice
from minitcad import discretization,physics
from minitcad.discretization.sg import edge_flux_qf
from minitcad.discretization.residual import residual_qf,ResidualEvaluation,evaluate_qf
'''
def _state(u:NDArray,fn:NDArray,fp:NDArray) -> State:
    n = np.exp(u + fn)
    p = np.exp(-u - fp)
    return State(u,n,p)
def _quasi_fermi(s:State) -> tuple[NDArray,NDArray]:
    return (np.log(s.n) - s.u,-np.log(s.p) - s.u)
'''
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
        n = np.exp(potential + fn)
        p = np.exp(-potential - fp)
        s = State(u = potential,n = n,p = p)
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
        qf : QFState,
        device : NormalizedDevice,
        bias_V : float,
        evaluation : ResidualEvaluation,
        physical : State | None = None
):
    if physical is None:
        physical = qf.to_physical()
    m = qf.node_count
    f = evaluation.f.reshape(m,3)
    flux = evaluation.fluxes
    rate = evaluation.recombination_nd
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
        np.abs(flux.jn_edges[1:]) / w,
        np.abs(flux.jn_edges[:-1]) / w,
        np.abs(rate[1:-1])
    ))
    scale[1:-1,2] = np.maximum.reduce((
        np.ones(m - 2),
        np.abs(flux.jp_edges[1:]) / w,
        np.abs(flux.jp_edges[:-1]) / w,
        np.abs(rate[1:-1])
    ))
    targets = contact_targets_qf(device = device,bias_V = bias_V,ref_n = qf.ref_n,ref_p = qf.ref_p)
    scale[0] = np.maximum(
        1.0,
        np.abs(targets[:3])
    )
    scale[-1] = np.maximum(
        1.0,
        np.abs(targets[3:])
    )
    rn = float(np.max(np.abs(f) / scale))
    if not np.isfinite(rn):
        raise FloatingPointError("NONFINITE_SCALED_RESIDUAL")
    return rn
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
            ref_p = ref_p
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
def _project_contacts(
        qf : QFState,
        targets : tuple
) -> QFState:
    u = qf.u.copy()
    eta_n = qf.eta_n.copy()
    eta_p = qf.eta_p.copy()
    u[0],u[-1] = targets[0],targets[3]
    eta_n[0],eta_n[-1] = targets[1],targets[4]
    eta_p[0],eta_p[-1] = targets[2],targets[5]
    return QFState(
        u = u,
        eta_n = eta_n,
        eta_p = eta_p,
        ref_n = qf.ref_n,
        ref_p = qf.ref_p
    )
def _interpolate_qf(
    old : QFState,
    residual : QFState,
    alpha : float
) -> QFState:
    if (
        old.ref_n != residual.ref_n
        or old.ref_p != residual.ref_p
    ):
        raise ValueError(
            "QF reference changed during iteration"
        )
    return QFState(
        u = old.u + alpha * (
            residual.u - old.u
        ),
        eta_n = old.eta_n + alpha * (
            residual.eta_n - old.eta_n
        ),
        eta_p = old.eta_p + alpha * (
            residual.eta_p - old.eta_p
        ),
        ref_n = old.ref_n,
        ref_p = old.ref_p
    )
def _physics_from_evaluation(ev : ResidualEvaluation) -> PhysicsDiagnostics:
    return PhysicsDiagnostics(
        fluxes = ev.fluxes,
        recombination_nd = ev.recombination_nd,
        residual_qf_nd = ev.f
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
def solve_gummel(
    device : NormalizedDevice,
    bias_V : float,
    initial_state : State,
    *,
    tol_residual : float = 1e-8,
    tol_update : float = 1e-6,
    max_iter : int = 60,
    alpha_min : float = 1e-6,
    n_floor : float = 1e-15,
    diagnostics : dict | None = None,
) -> SolveResult:
    """
    QF-based Gummel solver.
    Public input:
        State(u, n, p)
    Internal unknowns:
        QFState(u, eta_n, eta_p, ref_n, ref_p)
    Public output:
        SolveResult containing physical State and
        stable PhysicsDiagnostics.
    """
    if initial_state.node_count != device.mesh.node_count:
        raise ValueError(
            "device/state node_count mismatch"
        )
    for name,value in (
        ("bias_V",bias_V),
        ("tol_residual",tol_residual),
        ("tol_update",tol_update),
        ("alpha_min",alpha_min),
        ("n_floor",n_floor)
    ):
        if isinstance(value,(bool,np.bool_)) or not isinstance(value,(int,float,np.integer,np.floating)) or not np.isfinite(value) or (name != "bias_V" and value <= 0):
            raise ValueError(f"invalid {name}")
    if alpha_min > 1.0:
        raise ValueError("alpha_min must be <= 1")
    if (isinstance(max_iter,bool) or not isinstance(max_iter,(int,np.integer)) or
        max_iter < 1
    ):
        raise ValueError("max_iter must be a positive integer")
    if diagnostics is not None:
        diagnostics.clear()
    targets = device.contact_targets(bias_V)#(uL,nL,pL,uR,nR,pR)
    ref_n = float(np.log(targets[4]) - targets[3])
    ref_p = float(-np.log(targets[2]) - targets[0])
    targets_qf = contact_targets_qf(
        device = device,
        bias_V = bias_V,
        ref_n = ref_n,
        ref_p = ref_p,
    )
    history : list[dict] = []
    def make_result(
        ok : bool,
        iterations : int,
        rn : float,
        dn : float,
        reason : str,
        qf : QFState | None = None,
        ev : ResidualEvaluation = None,
        physical : State | None = None,
    ) -> SolveResult:
        physics = None
        if ok:
            if qf is None or ev is None or physical is None:
                raise RuntimeError("Successful solve missing final state")
            physics = _physics_from_evaluation(ev)
            # Transitional debug API for older tests.
            if diagnostics is not None:
                diagnostics.update({
                    "u" : qf.u.copy(),
                    "eta_n" : qf.eta_n.copy(),
                    "eta_p" : qf.eta_p.copy(),
                    "ref_n" : float(qf.ref_n),
                    "ref_p" : float(qf.ref_p),
                    "jn_edges" : ev.fluxes.jn_edges.copy(),
                    "jp_edges" : ev.fluxes.jp_edges.copy(),
                    "rate" : ev.recombination_nd.copy(),
                    "residual" : ev.f.reshape(qf.node_count,3).copy()
                })
        return SolveResult(
            bias_V = float(bias_V),
            state = physical if ok else None,
            converged = ok,
            iterations = iterations,
            residual_norm = float(rn),
            update_norm = float(dn),
            reason = reason,
            iteration_history = tuple(history),
            physics = physics
        )
    try:
        qf = QFState.from_physical(state = initial_state,ref_n = ref_n,ref_p = ref_p)
        qf = _project_contacts(qf,targets_qf)
        physical = qf.to_physical()
        ev = evaluate_qf(state = qf,device = device,bias_V = bias_V)
        rn = _scaled_residual_qf(qf = qf,device = device,bias_V = bias_V,evaluation = ev,physical = physical)
    except (ValueError,FloatingPointError,OverflowError,ArithmeticError) as exc:
        history.append({
            "iteration" : 0,
            "status" : "INVALID_INITIAL_STATE",
            "reason" : str(exc)
        })
        return make_result(
            ok = False,
            iterations = 0,
            rn = np.inf,
            dn = np.inf,
            reason = "INVALID_INITIAL_STATE"
        )
    if rn <= tol_residual:
        return make_result(
            ok = True,
            iterations=  0,
            rn = rn,
            dn = 0.0,
            reason = "SUCCESS",
            qf = qf,
            ev = ev,
            physical = physical,
        )
    dn = np.inf
    for it in range(1, max_iter + 1):
        try:
            ustar = _poisson_step(
                device,
                qf.u,
                qf.eta_n + qf.ref_n,
                qf.eta_p + qf.ref_p,
                targets
            )
            enstar = _qf_carrier_step(
                device,
                ustar,
                qf.eta_n,
                qf.eta_p,#注意这里已经有了ref_n
                electron = True,
                ref_n = qf.ref_n,
                ref_p = qf.ref_p
            )
            epstar = _qf_carrier_step(
                device,
                ustar,
                enstar,
                qf.eta_p,
                electron = False,
                ref_n = qf.ref_n,
                ref_p = qf.ref_p
            )
            star = QFState(
                u = ustar,
                eta_n = enstar,
                eta_p = epstar,
                ref_n = qf.ref_n,
                ref_p = qf.ref_p
            )
            star = _project_contacts(star,targets_qf)
        except (ValueError,FloatingPointError,OverflowError,ArithmeticError) as exc:
            history.append({
                "iteration" : it,
                "status" : "SUBSTEP_FAILED",
                "reason" : str(exc)
            })
            return make_result(
                ok = False,
                iterations = it,
                rn = rn,
                dn = dn,
                reason = "SUBSTEP_FAILED: " + str(exc)
            )
        alpha = 1.0
        accepted = False
        while alpha >= alpha_min:
            trial_rn = np.inf
            trial_dn = np.inf
            trial_error = None
            try:
                trial_qf = _interpolate_qf(old = qf,residual = star,alpha = alpha)
                trial_physical = trial_qf.to_physical()
                trial_ev = evaluate_qf(trial_qf,device,bias_V)
                trial_rn = _scaled_residual_qf(qf = trial_qf,device = device,bias_V = bias_V,evaluation = trial_ev,physical = trial_physical)
                trial_dn = _scaled_update(old = physical,new = trial_physical,floor = n_floor)
                good = trial_rn < rn or trial_rn <= tol_residual
            except (ValueError,FloatingPointError,OverflowError,ArithmeticError) as exc:
                good = False
                trial_error = str(exc)
            history.append({
                "iteration" : it,
                "alpha" : float(alpha),
                "residual_norm" : float(trial_rn),
                "update_norm" : float(trial_dn),
                "status": (
                    "ACCEPTED"
                    if good
                    else "REJECTED"
                ),
                "reason" : trial_error
            })
            if good:
                qf = trial_qf
                physical = trial_physical
                ev = trial_ev
                rn = trial_rn
                dn = trial_dn
                accepted = True
                break
            alpha *= 0.5
        if not accepted:
            return make_result(
                ok = False,
                iterations = it,
                rn = rn,
                dn = dn,
                reason = "LINE_SEARCH_FAILED"
            )
        if rn <= tol_residual and dn <= tol_update:
            return make_result(
                ok = True,
                iterations=  it,
                rn = rn,
                dn = dn,
                reason = "SUCCESS",
                qf = qf,
                ev = ev,
                physical = physical
            )
    return make_result(
        ok = False,
        iterations = max_iter,
        rn = rn,
        dn = dn,
        reason = "MAX_ITER"
    )