#0V稳态下的电势求解
import numpy as np
from numpy.typing import ArrayLike,NDArray
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve
from minitcad.core.state import State,SolveResult
from minitcad import discretization
def solve_equilibrium_scalar(
        x_nodes_nd : ArrayLike,
        c_nodes_nd : ArrayLike,
        u_left_nd : np.float64 | float,
        u_right_nd : np.float64 | float,
        *,
        tol_residual = 1e-8,
        tol_update = 1e-6,
        max_iter = 60,
        alpha_min = 1e-6,
        n_floor = 1e-15
) -> SolveResult:
    x = np.asarray(x_nodes_nd)
    c = np.asarray(c_nodes_nd)
    if (x.ndim != 1 or x.size < 3):
        raise ValueError("invaild shape")
    m = x.size
    if (c.shape != (m,) and c.shape != (m - 2,)):
        raise ValueError("cbar must have length M or M-2")
    if not (np.isfinite(x).all() and np.isfinite(c).all()):
        raise ValueError("nonfinite input")
    h = np.diff(x)
    w = 0.5 * (h[:-1] + h[1:])
    c_in = c[1:-1] if c.size == m else c
    if (np.any(h <= 0)):
        raise ValueError("xi must be strictly increasing")
    if not np.isfinite([u_left_nd,u_right_nd]).all():
        raise ValueError("contact potentials must be finite")
    if not (np.isfinite(tol_residual) and tol_residual > 0):
        raise ValueError("invalid tol_residual")
    if not (np.isfinite(tol_update) and tol_update > 0):
        raise ValueError("invalid tol_update")
    if (
        not isinstance(max_iter,(int, np.integer))
        or isinstance(max_iter,(bool, np.bool_))
        or max_iter < 1
    ):
        raise ValueError("max_iter must be positive integer")
    if not (np.isfinite(alpha_min) and 0 < alpha_min <= 1):
        raise ValueError("invalid alpha_min")
    if not (np.isfinite(n_floor) and n_floor > 0):
        raise ValueError("invalid n_floor")
    u = np.empty(m,dtype = np.float64)
    u[1:-1] = np.arcsinh(c_in / 2.0)
    u[0] = u_left_nd
    u[-1] = u_right_nd
    def evaluate(potential:NDArray):
        n = np.exp(potential)
        p = np.exp(-potential)
        if not (np.isfinite(n).all() and np.isfinite(p).all() 
        and np.all(n > 0) and np.all(p > 0)):
            raise FloatingPointError("non positive or non finite carrier")
        f = discretization.poisson_residual_nd(xi = x,u = potential,n = n,p = p,cbar = c_in,
                                               u_left = u_left_nd,u_right = u_right_nd)
        scale = np.empty(m,dtype = np.float64)
        scale[0] = max(1.0,abs(u_left_nd))
        scale[-1] = max(u_right_nd,1.0)
        scale[1:-1] = np.maximum.reduce((np.ones(m - 2),np.abs(c_in),n[1:-1],p[1:-1]))
        rnorm = float(np.max(np.abs(f) / scale))
        if not np.isfinite(rnorm):
            raise FloatingPointError("non finite residual norm")
        state = State(
            u = potential,
            n = n,
            p = p
        )
        return state,f,rnorm
    def relative_update(old:State,new:State) -> float:
        su = np.maximum.reduce((
            (np.ones(m),
            np.abs(old.u),
            np.abs(new.u))
        ))
        sn = np.maximum.reduce(
            (np.full(m,n_floor),
            old.n,
            new.n)
        )
        sp = np.maximum.reduce(
            (np.full(m,n_floor),
            old.p,
            new.p)
        )
        return float(max(
        np.max(np.abs(new.u - old.u) / su),
        np.max(np.abs(new.n - old.n) / sn),
        np.max(np.abs(new.p - old.p) / sp),
        ))
    def make_result(
        ok,state,iterations,
        rnorm,dnorm,reason,history) -> SolveResult:
        return SolveResult(
            bias_V = 0.0,
            state = state if ok else None,
            converged = ok,
            iterations = iterations,
            residual_norm = float(rnorm),
            update_norm = float(dnorm),
            reason = reason,
            iteration_history = tuple(history),
    )
    history : list[dict] = []
    try:
        state,f,rnorm = evaluate(u)
    except(ValueError,FloatingPointError,OverflowError):
        return make_result(False,None,0,np.inf,np.inf,"invalid initial state",history)
    if rnorm <= tol_residual:
        return make_result(True,state,0,rnorm,0.0,"SUCCESS",history)
    dnorm = np.inf
    for step in range(1,max_iter + 1):
        diagonal = - ((1.0 / h[:-1] + 1.0 / h[1:]) / w + state.n[1:-1] + state.p[1:-1])
        lower = 1.0 / (w[1:] * h[1:-1])
        upper = 1.0 / (w[:-1] * h[1:-1])
        jacs = diags(
            (diagonal,lower,upper),
            offsets = (0,-1,1),
            shape = (m - 2,m - 2),
            format = "csc"
        )
        try:
            delta = spsolve(jacs,-f[1:-1])
            if (delta.shape != (m - 2,) or not np.isfinite(delta).all()):
                raise FloatingPointError("BAD_LINEAR_STEP")
        except:
            return make_result(False,state,step,rnorm,dnorm,"LINEAR_SOLVE_FAILED",history)
        alpha = 1.0
        accepted = False
        while (alpha >= alpha_min):
            candidate_u = state.u.copy()
            candidate_u[1:-1] += delta * alpha
            new_rnorm = np.inf
            new_dnorm = np.inf
            try:
                new_state,new_f,new_rnorm = evaluate(candidate_u)
                new_dnorm = relative_update(state,new_state)
                acceptable = (
                    new_rnorm < rnorm
                    or new_rnorm <= tol_residual
                )
            except:
                new_state = None
                new_rnorm = np.inf
                new_dnorm = np.inf
                acceptable = False
            history.append({
                "step":step,
                "alpha":float(alpha),
                "residual_norm":(
                    float(new_rnorm)
                    if np.isfinite(new_rnorm)
                    else None
                ),
                "update_norm":(
                    float(new_dnorm)
                    if np.isfinite(new_dnorm)
                    else None
                ),
                "trial_status":(
                    "ACCEPTED"
                    if acceptable
                    else "REJECTED"
                )
            })
            if acceptable:
                accepted = True
                break
            alpha *= 0.5
        if not accepted:
            return make_result(
                False,None,step,
                rnorm,dnorm,
                "line search fail",history
            )
        state = new_state
        f = new_f
        rnorm = new_rnorm
        dnorm = new_dnorm
        if (
            rnorm <= tol_residual
            and dnorm <= tol_update
        ):
            return make_result(
                True,state,step,
                rnorm,dnorm,
                "SUCCESS",history
            )
    return make_result(
        False,None,max_iter,
        rnorm,dnorm,
        "MAX_ITER",history)
solve_equilibrium_nd = solve_equilibrium_scalar