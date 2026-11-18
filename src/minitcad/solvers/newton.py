"""Newton 联立迭代（QF 口径）：对 `F(z)` 直接解耦合系统。

与 `gummel.py` 的区别：Gummel 拆成三个子问题轮流解，收敛是线性的；Newton 一次解
`3M × 3M` 的联立系统，在解的邻域内二次收敛。

解之前先把每行除以 `residual_scales_qf`：Poisson 行是 `1e6` 量级、连续性行是 `1` 量级，
不缩放的话 `dz` 的精度被大行主导，连续性方向的修正会落到 `1e-16` 以下。
判据用的也是同一组尺度。

失败（线性求解、线搜索耗尽、迭代超限）一律返回 `converged=False` 且 `state=None`，
不抛异常；参数非法才抛 `ValueError`。
"""

from __future__ import annotations

import warnings

import numpy as np
from scipy.linalg import norm as stable_l2_norm
from scipy.sparse import diags
from scipy.sparse.linalg import MatrixRankWarning, spsolve

from minitcad.core.contact import contact_targets_qf
from minitcad.core.device import NormalizedDevice
from minitcad.core.state import PhysicsDiagnostics, QFState, SolveResult, State
from minitcad.discretization.residual import ResidualEvaluation, evaluate_qf
from minitcad.solvers.convergence import (
    residual_scales_qf,
    scaled_residual_qf,
    scaled_update_physical,
)
from minitcad.solvers.jacobian import jacobian_qf

# 试探点溢出是预期路径，由线搜索拒绝，不算程序错误。
_EXPECTED_FLOAT_ERRORS = (ValueError, FloatingPointError, OverflowError, ArithmeticError)

# 后向误差 `‖J dz + F‖∞ / max(1, ‖F‖∞)` 的上限：超限说明解不可信，继续线搜索只是白费。
_BACKWARD_ERROR_MAX = 1e-6


def _check_settings(
    device: NormalizedDevice,
    bias_V: float,
    initial_state: State,
    tol_residual: float,
    tol_update: float,
    max_iter: int,
    alpha_min: float,
    n_floor: float,
) -> None:
    """参数校验；非法输入抛 `ValueError`，与 `solve_gummel` 保持同一组规则。"""
    if initial_state.node_count != device.mesh.node_count:
        raise ValueError("device/state node_count mismatch")
    for name, value in (
        ("bias_V", bias_V),
        ("tol_residual", tol_residual),
        ("tol_update", tol_update),
        ("alpha_min", alpha_min),
        ("n_floor", n_floor),
    ):
        if (
            isinstance(value, (bool, np.bool_))
            or not isinstance(value, (int, float, np.integer, np.floating))
            or not np.isfinite(value)
            or (name != "bias_V" and value <= 0.0)
        ):
            raise ValueError(f"invalid {name}")
    if alpha_min > 1.0:
        raise ValueError("alpha_min must be <= 1")
    if isinstance(max_iter, (bool, np.bool_)) or not isinstance(max_iter, (int, np.integer)) or max_iter < 1:
        raise ValueError("max_iter must be a positive integer")


def _project_contacts(qf: QFState, targets: tuple[np.float64, ...]) -> QFState:
    """把两端的 QF 未知量钉到 Dirichlet 目标上：接触单元里它们是给定值，不是未知量。

    在进入迭代前置好，可以让接触行的初始残差为 0。
    """
    z = qf.pack_qf()
    z[0:3] = targets[0:3]
    z[-3:] = targets[-3:]
    return QFState.unpack_qf(z=z, node_count=qf.node_count, ref_n=qf.ref_n, ref_p=qf.ref_p)


def _physics_from_evaluation(evaluation: ResidualEvaluation) -> PhysicsDiagnostics:
    """裁出对外的诊断包；QF 迭代变量不外传。

    与 `gummel.py` 的同名函数重复是有意的：两个求解器要能各自独立替换。
    """
    return PhysicsDiagnostics(
        fluxes=evaluation.fluxes,
        recombination_nd=evaluation.recombination_nd,
        residual_qf_nd=evaluation.f,
    )


def solve_newton(
    device: NormalizedDevice,
    bias_V: float,
    initial_state: State,
    *,
    tol_residual: float = 1e-8,
    tol_update: float = 1e-6,
    max_iter: int = 60,
    alpha_min: float = 1e-6,
    n_floor: float = 1e-15,
    diagnostics: dict | None = None,
) -> SolveResult:
    """解一个偏压点 `bias_V`，输入输出都是物理状态 `State(u, n, p)`。

    签名与 `solve_gummel` 逐参数对齐，这样 `simulation/dc.py` 只需按 `solver.method`
    选入口。`diagnostics` 也同义：成功时写入内部 QF 状态与最终通量/残差。
    """
    _check_settings(
        device, bias_V, initial_state, tol_residual, tol_update, max_iter, alpha_min, n_floor
    )
    if diagnostics is not None:
        diagnostics.clear()

    targets = device.contact_targets(bias_V)
    # 锚点与 solve_gummel 相同（右端电子、左端空穴）：两边不同的话，同一物理状态会
    # 对应两组 eta，双求解器对照还得先平移。
    ref_n = float(np.log(targets[4]) - targets[3])
    ref_p = float(-np.log(targets[2]) - targets[0])
    qf_targets = contact_targets_qf(device, bias_V, ref_n=ref_n, ref_p=ref_p)

    history: list[dict] = []

    def make_result(
        ok: bool,
        iterations: int,
        rn: float,
        dn: float,
        reason: str,
        physical: State | None = None,
        evaluation: ResidualEvaluation | None = None,
        qf: QFState | None = None,
    ) -> SolveResult:
        physics = None
        if ok:
            if physical is None or evaluation is None or qf is None:
                raise RuntimeError("successful Newton solve is missing the final state")
            physics = _physics_from_evaluation(evaluation)
            if diagnostics is not None:
                diagnostics.update(
                    {
                        "u": qf.u.copy(),
                        "eta_n": qf.eta_n.copy(),
                        "eta_p": qf.eta_p.copy(),
                        "ref_n": float(qf.ref_n),
                        "ref_p": float(qf.ref_p),
                        "jn_edges": evaluation.fluxes.jn_edges.copy(),
                        "jp_edges": evaluation.fluxes.jp_edges.copy(),
                        "rate": evaluation.recombination_nd.copy(),
                        "residual": evaluation.f.reshape(qf.node_count, 3).copy(),
                    }
                )
        return SolveResult(
            bias_V=float(bias_V),
            state=physical if ok else None,
            converged=ok,
            iterations=iterations,
            residual_norm=float(rn),
            update_norm=float(dn),
            reason=reason,
            iteration_history=tuple(history),
            physics=physics,
        )

    try:
        qf = QFState.from_physical(state=initial_state, ref_n=ref_n, ref_p=ref_p)
        qf = _project_contacts(qf, qf_targets)
        physical = qf.to_physical()
        evaluation = evaluate_qf(state=qf, device=device, bias_V=bias_V)
        rn = scaled_residual_qf(qf, device, bias_V, evaluation, physical)
    except _EXPECTED_FLOAT_ERRORS as exc:
        history.append({"iteration": 0, "status": "INVALID_INITIAL_STATE", "reason": str(exc)})
        return make_result(False, 0, np.inf, np.inf, "INVALID_INITIAL_STATE")

    if rn <= tol_residual:
        return make_result(True, 0, rn, 0.0, "SUCCESS", physical, evaluation, qf)

    dn = np.inf
    for iteration in range(1, max_iter + 1):
        try:
            jac = jacobian_qf(state=qf, device=device, bias_V=bias_V)
            scales = residual_scales_qf(qf, device, bias_V, evaluation, physical)
            inverse_scales = 1.0 / scales
            # 只缩放方程、不缩放未知量：要平衡的是行之间的量级差。
            jacobian_scaled = diags(inverse_scales, format="csc") @ jac
            rhs_scaled = -evaluation.f * inverse_scales
            with warnings.catch_warnings():
                # spsolve 遇到奇异矩阵只发警告并返回垃圾值，这里当成异常处理。
                warnings.simplefilter("error", MatrixRankWarning)
                dz = np.asarray(spsolve(jacobian_scaled, rhs_scaled), dtype=np.float64)
            if dz.shape != (3 * qf.node_count,) or not np.isfinite(dz).all():
                raise FloatingPointError("nonfinite Newton direction")
            # 后向误差用未缩放的 J、F：它问的是“dz 解没解对原方程”。
            linear_defect = jac @ dz + evaluation.f
            backward_error = np.linalg.norm(linear_defect, ord=np.inf) / max(
                1.0, np.linalg.norm(evaluation.f, ord=np.inf)
            )
            if not np.isfinite(backward_error) or backward_error > _BACKWARD_ERROR_MAX:
                raise FloatingPointError(
                    f"linear solve backward residual too large: {backward_error:.3e}"
                )
        except (*_EXPECTED_FLOAT_ERRORS, RuntimeError, MatrixRankWarning) as exc:
            history.append(
                {"iteration": iteration, "status": "LINEAR_SOLVE_FAILED", "reason": str(exc)}
            )
            return make_result(False, iteration, rn, dn, "LINEAR_SOLVE_FAILED")

        z = qf.pack_qf()
        base_merit = float(stable_l2_norm(evaluation.f / scales))
        alpha = 1.0
        accepted = False
        while alpha >= alpha_min:
            trial_rn = np.inf
            trial_dn = np.inf
            trial_merit = np.inf
            error = None
            try:
                trial = QFState.unpack_qf(
                    z=z + alpha * dz,
                    node_count=qf.node_count,
                    ref_n=qf.ref_n,
                    ref_p=qf.ref_p,
                )
                trial_physical = trial.to_physical()
                trial_evaluation = evaluate_qf(state=trial, device=device, bias_V=bias_V)
                trial_rn = scaled_residual_qf(trial, device, bias_V, trial_evaluation, trial_physical)
                trial_dn = scaled_update_physical(physical, trial_physical, n_floor=n_floor)
                trial_merit = float(stable_l2_norm(trial_evaluation.f / scales))
                good = trial_rn <= tol_residual or trial_merit <= (1.0 - 1e-4 * alpha) * base_merit
            except _EXPECTED_FLOAT_ERRORS as exc:
                good = False
                error = str(exc)
            history.append(
                {
                    "iteration": iteration,
                    "alpha": float(alpha),
                    "residual_norm": float(trial_rn),
                    "update_norm": float(trial_dn),
                    # 这两项不进 CSV，但线搜索失败时要能看出是 merit 不降还是步长到底。
                    "merit_L2": float(trial_merit),
                    "base_merit_L2": float(base_merit),
                    "linear_backward_error": float(backward_error),
                    "status": "ACCEPTED" if good else "REJECTED",
                    "reason": error,
                }
            )
            if good:
                qf, physical, evaluation = trial, trial_physical, trial_evaluation
                rn, dn = trial_rn, trial_dn
                accepted = True
                break
            alpha *= 0.5
        if not accepted:
            return make_result(False, iteration, rn, dn, "LINE_SEARCH_FAILED")
        if rn <= tol_residual and dn <= tol_update:
            # 两个判据都要满足：残差够小但步长还大，说明停在平坦区域，还会漂。
            return make_result(True, iteration, rn, dn, "SUCCESS", physical, evaluation, qf)

    return make_result(False, max_iter, rn, dn, "MAX_ITER")
