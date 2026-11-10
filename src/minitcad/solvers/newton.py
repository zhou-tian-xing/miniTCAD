"""Newton 联立迭代（QF 口径）：对 `F(z)` 直接解耦合系统（W5-A1…A5）。

与 `gummel.py` 的分工：Gummel 把耦合方程拆成三个标量/三对角子问题轮流解，每步的内存与
时间都很省，但收敛是线性的——电流越大（正向偏压越高）耦合越强，迭代次数越多，最终在
某个偏压处失败。Newton 一次解 `3M × 3M` 的联立线性系统，在解的邻域内二次收敛，
代价是 Jacobian 的装配与稀疏分解。

**行缩放为什么不只是判据**：Poisson 行的量级是 `1e6`（掺杂浓度），连续性行是 `1` 左右。
若直接解 `J dz = -F`，`dz` 的精度被大行主导，连续性方向的修正会落到 `1e-16` 以下，
表现为“线性求解成功但外层不下降”。所以这里把每行除以 `residual_scales_qf` 再解，
使所有行同量级；判据用的也是同一组尺度，两者不会分叉。

**线搜索**：Newton 方向只有在解的邻域内才保证下降。这里用 scaled-residual 的 L2 范数做
merit 函数，要求 `merit` 至少按 `1 - 1e-4*alpha` 下降（Armijo 条件的宽松版）；不合格就把
`alpha` 减半，直到 `alpha_min`。`n` 在试探点上溢出是正常的，被拒绝的试探不写入状态。

**失败语义**：线性求解失败（矩阵奇异/条件数过高/后向误差过大）、线搜索耗尽、迭代超限
都返回 `converged=False` 且 `state=None`；`iteration_history` 保留每次试探的记录，
失败原因写在 `reason`。本层不抛异常，参数非法才抛 `ValueError`。
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

# 试探点溢出、解出非有限值都是预期路径：线搜索会拒绝它，改由状态构造与
# scaled_residual_qf 的有限性检查判定失败。
_EXPECTED_FLOAT_ERRORS = (ValueError, FloatingPointError, OverflowError, ArithmeticError)

# 线性解的**后向误差**上限：`‖J dz + F‖∞ / max(1, ‖F‖∞)`。
# 超限说明线性系统病态到解不可信，此时继续线搜索只会白费评估，直接判失败更快也更诚实。
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
    """把两端的 QF 未知量钉到 Dirichlet 目标上。

    `u` 与 `eta_n/eta_p` 在接触处是给定值而不是未知量；`F` 的前 3 行与后 3 行就是
    “当前值 - 目标值”。把它们在**进入迭代前**就置为目标值，可以让初始残差直接为 0，
    也避免 Newton 在第一次迭代里把接触行的修正摊到内部结点上。
    """
    z = qf.pack_qf()
    z[0:3] = targets[0:3]
    z[-3:] = targets[-3:]
    return QFState.unpack_qf(z=z, node_count=qf.node_count, ref_n=qf.ref_n, ref_p=qf.ref_p)


def _physics_from_evaluation(evaluation: ResidualEvaluation) -> PhysicsDiagnostics:
    """与 `gummel.py` 的同名函数等价：裁出对外诊断包，不外传 QF 迭代变量。

    故意各自保留一份而不是共享：`solvers` 的两个文件应能独立替换（A 侧随时可能改动
    其中一个），共享一个 6 行辅助函数会把两者重新绑在一起。
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

    签名与 `solve_gummel` 逐参数对齐（`docs/a_b_interface_week02.md` §3），
    这样 `simulation/dc.py` 可以只按 `solver.method` 选择入口，扫描、取消、写文件的
    代码完全不变。`diagnostics` 的语义也与 Gummel 相同：成功时写入内部 QF 状态与最终
    通量/残差，供双求解器对照核对两边收敛到的是同一个解。
    """
    _check_settings(
        device, bias_V, initial_state, tol_residual, tol_update, max_iter, alpha_min, n_floor
    )
    if diagnostics is not None:
        diagnostics.clear()

    targets = device.contact_targets(bias_V)
    # 锚点与 `solve_gummel` 完全一致（右端电子、左端空穴）。两边必须用同一锚点，
    # 否则同一物理状态会对应两组不同的 eta，双求解器对照要先做一次平移才能比较。
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
            # 行缩放写成 D^-1 J 而不是 J（只对未知量缩放）：要平衡的是“方程之间的量级差”，
            # 变量本身的尺度已经由行缩放的逆作用回 dz 上。
            jacobian_scaled = diags(inverse_scales, format="csc") @ jac
            rhs_scaled = -evaluation.f * inverse_scales
            with warnings.catch_warnings():
                # 奇异矩阵会让 spsolve 只发 MatrixRankWarning 并返回垃圾值，
                # 这里把它当成异常，与其它数值失败走同一条退出路径。
                warnings.simplefilter("error", MatrixRankWarning)
                dz = np.asarray(spsolve(jacobian_scaled, rhs_scaled), dtype=np.float64)
            if dz.shape != (3 * qf.node_count,) or not np.isfinite(dz).all():
                raise FloatingPointError("nonfinite Newton direction")
            # 后向误差用**未缩放**的 J 与 F 计算：它度量的是“dz 到底解没解对原方程”，
            # 在缩放后的系统上算会把行缩放的量级也带进来。
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
                    # 下面两项不进 CSV（`iteration` 行只写 residual/update/alpha），
                    # 但线搜索失败时必须能看出“是 merit 不降还是步长已到底”。
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
            # 两个判据同时对才收：残差够小但步长还大，说明停在了一个平坦区域，
            # 下一步仍可能大幅漂移（与 Gummel 共用 `is_converged` 的口径）。
            return make_result(True, iteration, rn, dn, "SUCCESS", physical, evaluation, qf)

    return make_result(False, max_iter, rn, dn, "MAX_ITER")
