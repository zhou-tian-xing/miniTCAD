#Newton迭代
#对Newton系统进行残差行缩放，以平衡不同物理方程的数值尺度。2026/10/24
import warnings
import numpy as np
from scipy.sparse import diags
from scipy.linalg import norm as stable_l2_norm
from scipy.sparse.linalg import spsolve,MatrixRankWarning
from minitcad.core.contact import contact_targets_qf
from minitcad.core.device import NormalizedDevice
from minitcad.core.state import PhysicsDiagnostics,QFState,SolveResult,State
from minitcad.discretization.residual import ResidualEvaluation,evaluate_qf
from minitcad.solvers.convergence import scaled_residual_qf,scaled_update_physical,residual_scales_qf
from minitcad.solvers.jacobian import jacobian_qf
def _check_settings(device,bias_V,initial_state,tol_residual,tol_update,max_iter,alpha_min,n_floor):
    if initial_state.node_count != device.mesh.node_count:
        raise ValueError("device/state node_count mismatch")
    for name,value in (
        ("bias_V",bias_V),
        ("tol_residual",tol_residual),
        ("tol_update",tol_update),
        ("alpha_min",alpha_min),
        ("n_floor",n_floor)
    ):
        if (isinstance(value,(bool,np.bool_)) or not isinstance(value,(int,float,np.integer,np.floating)) or
            not np.isfinite(value) or (name != "bias_V" and value <= 0.0)):
            raise ValueError(f"invalid {name}")
    if alpha_min > 1.0:
        raise ValueError("alpha_min must be <= 1")
    if isinstance(max_iter,(bool,np.bool_)) or not isinstance(max_iter,(int,np.integer)) or max_iter < 1:
        raise ValueError("max_iter must be a positive integer")
def _project_contacts(qf:QFState,targets:tuple[np.float64,...]) -> QFState:
    z = qf.pack_qf()
    z[0:3] = targets[0:3]
    z[-3:] = targets[-3:]
    return QFState.unpack_qf(z = z,node_count = qf.node_count,ref_n = qf.ref_n,ref_p = qf.ref_p)
def _make_physics(ev : ResidualEvaluation) -> PhysicsDiagnostics:
    return PhysicsDiagnostics(fluxes = ev.fluxes,recombination_nd = ev.recombination_nd,residual_qf_nd = ev.f)
def solve_newton(
        device : NormalizedDevice,
        bias_V : float,
        initial_state : State,
        *,
        tol_residual : float = 1e-8,
        tol_update : float = 1e-6,
        max_iter : int = 60,
        alpha_min : float = 1e-6,
        n_floor : float = 1e-15
) -> SolveResult:
    _check_settings(device,bias_V,initial_state,tol_residual,tol_update,max_iter,alpha_min,n_floor)
    t = device.contact_targets(bias_V)
    ref_n = float(np.log(t[4]) - t[3])
    ref_p = float(-np.log(t[2]) - t[0])
    qf_targets = contact_targets_qf(device,bias_V,ref_n = ref_n,ref_p = ref_p)
    history : list[dict] = []
    def result(ok:bool,it:int,rn:float,dn:float,reason:str,physical : State | None = None,ev : ResidualEvaluation | None = None) -> SolveResult:
        physics = _make_physics(ev) if ok and ev is not None else None
        return SolveResult(bias_V = float(bias_V),state = physical if ok else None,converged = ok,iterations = it,
                        residual_norm = float(rn),update_norm = float(dn),reason = reason,iteration_history = tuple(history),physics = physics)
    try:
        qf = QFState.from_physical(state = initial_state,ref_n = ref_n,ref_p = ref_p)
        qf = _project_contacts(qf,qf_targets)
        physical = qf.to_physical()
        ev = evaluate_qf(state = qf,device = device,bias_V = bias_V)
        rn = scaled_residual_qf(qf = qf,device = device,bias_V = bias_V,evaluation = ev,physical = physical)
    except (ValueError,FloatingPointError,OverflowError,ArithmeticError) as exc:
        history.append({"iteration" : 0,
                        "status" : "INVALID_INITIAL_STATE",
                        "reason" : str(exc)})
        return result(False,0,np.inf,np.inf,"INVALID_INITIAL_STATE")
    if rn <= tol_residual:
        return result(True,0,rn,0.0,"SUCCESS",physical,ev)
    dn = np.inf
    for iteration in range(1, max_iter + 1):
        try:
            jac = jacobian_qf(state = qf,device = device,bias_V = bias_V)
            scales = residual_scales_qf(qf = qf,device = device,bias_V = bias_V,evaluation = ev,physical = physical)
            inv_scales = 1.0 / scales
            jac_scaled = diags(inv_scales,format = "csc") @ jac
            rhs_scaled = -ev.f * inv_scales
            with warnings.catch_warnings():
                warnings.simplefilter("error",MatrixRankWarning)
                dz = np.asarray(spsolve(jac_scaled,rhs_scaled),dtype = np.float64)
            if dz.shape != (3 * qf.node_count,) or not np.isfinite(dz).all():
                raise FloatingPointError("nonfinite Newton direction")
            linear_defect = jac @ dz + ev.f
            backward_error = np.linalg.norm(linear_defect,ord = np.inf) / max(1.0,np.linalg.norm(ev.f,ord = np.inf))
            if not np.isfinite(backward_error) or backward_error > 1e-6:
                raise FloatingPointError("linear solve backward residual too large: "
                                        f"{backward_error:.3e}")
        except (ValueError,FloatingPointError,OverflowError,ArithmeticError,RuntimeError,MatrixRankWarning) as exc:
            history.append({"iteration" : iteration,
                            "status" : "LINEAR_SOLVE_FAILED",
                            "reason" : str(exc)})
            return result(False,iteration,rn,dn,"LINEAR_SOLVE_FAILED")
        z = qf.pack_qf()
        base_scales = residual_scales_qf(qf,device,bias_V,ev,physical)
        base_merit = float(stable_l2_norm(ev.f / base_scales))
        alpha = 1.0
        accepted = False
        while alpha >= alpha_min:
            trial_rn = np.inf
            trial_dn = np.inf
            trial_merit = np.inf
            error = None
            try:
                trial = QFState.unpack_qf(z = z + alpha * dz,node_count = qf.node_count,ref_n = qf.ref_n,ref_p = qf.ref_p)
                trial_physical = trial.to_physical()
                trial_ev = evaluate_qf(state = trial,device = device,bias_V = bias_V)
                trial_rn = scaled_residual_qf(qf = trial,device = device,bias_V = bias_V,evaluation = trial_ev,physical = trial_physical)
                trial_dn = scaled_update_physical(old = physical,new = trial_physical,n_floor = n_floor)
                trial_merit = float(stable_l2_norm(trial_ev.f / base_scales))
                good = (trial_rn <= tol_residual or
                        trial_merit <= (1.0 - 1e-4 * alpha) * base_merit)
            except (ValueError,FloatingPointError,OverflowError,ArithmeticError) as exc:
                good = False
                error = str(exc)
            history.append({
                "iteration" : iteration,
                "alpha" : float(alpha),
                "residual_norm" : float(trial_rn),
                "update_norm" : float(trial_dn),
                "merit_L2" : float(trial_merit),
                "base_merit_L2" : float(base_merit),
                "linear_backward_error" : float(backward_error),
                "status" : "ACCEPTED" if good else "REJECTED",
                "reason" : error
            })
            if good:
                qf,physical,ev = trial,trial_physical,trial_ev
                rn,dn = trial_rn,trial_dn
                accepted = True
                break
            alpha *= 0.5
        if not accepted:
            return result(False,iteration,rn,dn,"LINE_SEARCH_FAILED")
        if rn <= tol_residual and dn <= tol_update:
            return result(True,iteration,rn,dn,"SUCCESS",physical,ev)
    return result(False,max_iter,rn,dn,"MAX_ITER")