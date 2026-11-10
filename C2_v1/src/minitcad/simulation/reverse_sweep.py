#自适应反向偏压求解
from dataclasses import dataclass
import numpy as np
from minitcad.core.device import NormalizedDevice
from minitcad.core.state import State,SolveResult
from minitcad.solvers import solve_gummel,solve_newton
@dataclass(frozen=True,slots = True)
class ReversePoint:
    bias_V : float
    solver_name : str
    result : SolveResult
@dataclass(frozen=True,slots = True)
class ReverseSweep:
    points : tuple[ReversePoint, ...]
    reached_target : bool
    failures : tuple[dict, ...]
def _valid(
    device : NormalizedDevice,
    bias : float,
    result : SolveResult,
    rtol : float,
    utol : float
) -> bool:
    if not (result.converged and result.state is not None and result.physics is not None and
             result.residual_norm <= rtol and result.update_norm <= utol):
        return False
    state = result.state
    if not (np.isfinite(state.u).all() and np.isfinite(state.n).all() and np.isfinite(state.p).all() and
             np.all(state.n > 0) and np.all(state.p > 0)):
        return False
    targets = device.contact_targets(bias)
    actual = [
        state.u[0],
        state.n[0],
        state.p[0],
        state.u[-1],
        state.n[-1],
        state.p[-1]
    ]
    if not np.allclose(actual,targets,rtol=1e-10,atol=1e-10):
        return False
    flux = result.physics.fluxes
    total = flux.jn_edges + flux.jp_edges
    jscale = max(
        1e-10,
        float(np.max(np.abs(total)))
    )
    spread = np.max(np.abs(total - np.median(total)))
    if spread > (1e-6 * jscale + 1e-10):
        return False
    return True
def solve_reverse_sweep(
    device : NormalizedDevice,
    equilibrium_state : State,
    target_V : float,
    *,
    initial_step_V : float = 0.05,
    min_step_V : float = 0.001,
    max_step_V : float = 0.10,
    tol_residual : float = 1e-8,
    tol_update : float = 1e-6,
    max_iter : int = 60
) -> ReverseSweep:
    if not (np.isfinite(target_V) and target_V < 0):
        raise ValueError("target_V must be finite and negative")
    settings = [
        initial_step_V,
        min_step_V,
        max_step_V,
        tol_residual,
        tol_update
    ]
    if not (np.isfinite(settings).all() and 0 < min_step_V <= initial_step_V <= max_step_V 
            and tol_residual > 0 and tol_update > 0):
        raise ValueError("invalid sweep settings")
    if not isinstance(max_iter, int) or isinstance(max_iter, bool) or max_iter <= 0:
        raise ValueError("max_iter must be positive int")
    initial = solve_newton(device = device,bias_V = 0.0,initial_state = equilibrium_state,
                           tol_residual = tol_residual,tol_update = tol_update,max_iter = max_iter)
    if not _valid(device = device,bias = 0.0,result = initial,rtol = tol_residual,utol = tol_update):
        raise ArithmeticError("zero-bias state must first be validated")
    accepted = [
        ReversePoint(0.0,"newton",initial)
    ]
    failures = []
    step = initial_step_V
    at = 0.0
    while at > target_V + 1e-12:
        nxt = max(target_V,at - step)
        solution = None
        solver_name = ""
        for name,solver in (
            ("newton",solve_newton),
            ("gummel",solve_gummel)
        ):
            try:
                out = solver(
                    device = device,
                    bias_V = nxt,
                    initial_state = accepted[-1].result.state,
                    tol_residual = tol_residual,
                    tol_update = tol_update,
                    max_iter = max_iter
                )
                if _valid(
                    device = device,bias = nxt,result = out,
                    rtol = tol_residual,utol = tol_update
                ):
                    solution = out
                    solver_name = name
                    break
                failures.append({
                    "from_V" : at,
                    "tried_V" : nxt,
                    "solver" : name,
                    "reason" : out.reason,
                    "residual" : out.residual_norm
                })
            except (ValueError,ArithmeticError,FloatingPointError,OverflowError) as exc:
                failures.append({
                    "from_V" : at,
                    "tried_V" : nxt,
                    "solver" : name,
                    "reason" : repr(exc)
                })
        if solution is not None:
            at = nxt
            accepted.append(
                ReversePoint(
                    at,
                    solver_name,
                    solution
                )
            )
            if (solver_name == "newton" and solution.iterations <= 8):
                step = min(max_step_V,step * 1.25,)
        else:
            step *= 0.5
            if step < min_step_V:
                return ReverseSweep(
                    tuple(accepted),
                    False,
                    tuple(failures)
                )
    return ReverseSweep(
        tuple(accepted),
        True,
        tuple(failures)
    )