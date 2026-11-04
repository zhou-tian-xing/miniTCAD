"""Gummel 解耦迭代（QF 口径，第 8.2 节）。

每个外层迭代把耦合方程拆成三步，每步只解一个变量：

1. Poisson：固定 `eta_n, eta_p`，解 `u`——非线性标量方程，内部 Newton + 阻尼线搜索；
2. 电子连续性：固定 `u` 与旧 `eta_p`，解 `eta_n`——线性化后的三对角方程；
3. 空穴连续性：固定 `u` 与新的 `eta_n`，解 `eta_p`——同 2。

三步之后得到 `star`，再对外层做一次线搜索：`F(star)` 的 scaled residual 必须下降才接受。
收敛要同时满足 scaled residual 与 scaled update（`solvers.convergence`），只看残差会把
停滞的迭代误判成收敛。

**为什么未知量取 `(u, eta_n, eta_p)` 而不是 `(u, n, p)`**：连续性方程的残差含 `np - 1`，
在 `np ≈ 1` 处两个大数相减会丢光有效位；改成 QF 变量后 `np - 1 = expm1(ref_n - ref_p + eta_n - eta_p)`，
直接对小的指数差求值。SG 通量同理，用 `expm1(Δeta)` 避开相邻结点浓渡差分的抵消。

**三对角从哪来**：控制体离散后，第 `i` 行只含 `i-1, i, i+1` 三个未知量——
边通量 `j_{i±1/2}` 各只依赖一侧相邻结点，复合项 `r_i` 只依赖本结点。
所以每个子步的 Jacobian 是 `(M-2)×(M-2)` 三对角矩阵，用 `scipy.sparse` 的 CSC 三对角装配后直接求解。
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve

from minitcad.core.contact import contact_targets_qf
from minitcad.core.device import NormalizedDevice
from minitcad.core.state import PhysicsDiagnostics, QFState, SolveResult, State
from minitcad.discretization import bernoulli, poisson_residual_nd
from minitcad.discretization.residual import ResidualEvaluation, evaluate_qf
from minitcad.discretization.sg import edge_flux_qf
from minitcad.solvers.convergence import is_converged, scaled_residual_qf, scaled_update

# 求解过程中 trial 点溢出是预期的：外层线搜索会拒绝它。
# 屏蔽 numpy 的 overflow 告警，但非有限值仍由 State/Fluxes 的构造检查与
# scaled_residual_qf 的有限性检查拦下。
_EXPECTED_FLOAT_ERRORS = (ValueError, FloatingPointError, OverflowError, ArithmeticError)


def _solve_tridiagonal(diagonal, left, right, rhs) -> NDArray[np.float64]:
    """解 `A x = b`，`A` 的主对角为 `diagonal`、下次对角为 `left`、上次对角为 `right`。

    `left[0]` 与 `right[-1]` 落在矩阵外，调用方不必初始化它们的有效值。
    """
    matrix = diags((left[1:], diagonal, right[:-1]), offsets=(-1, 0, 1), format="csc")
    x = np.asarray(spsolve(matrix, rhs), dtype=np.float64)
    if x.shape != rhs.shape or not np.isfinite(x).all():
        raise FloatingPointError("LINEAR_SOLVE_FAILED: invalid linear solution")
    return x


def _poisson_step(
    device: NormalizedDevice,
    u0: NDArray[np.float64],
    phi_n: NDArray[np.float64],
    phi_p: NDArray[np.float64],
    targets: tuple,
    *,
    max_iter: int = 35,
) -> NDArray[np.float64]:
    """固定 QF 势 `phi_n = ref_n + eta_n`、`phi_p = ref_p + eta_p`，解 Poisson 的 `u`。

    Newton 的 Jacobian 就是 Poisson 行对 `u` 的导数：
    `dF_i/du_i = -(1/h_i + 1/h_{i-1})/w_i - n_i - p_i`（`dp/du = -p`、`dn/du = n`），
    `dF_i/du_{i±1} = 1/(w_i h_{i±1/2})`。两端是 Dirichlet，不在未知量里。
    """
    h = device.mesh.h_edges_nd
    w = device.mesh.w_interior_nd
    m = device.mesh.node_count
    u = u0.copy()
    u[0], u[-1] = targets[0], targets[3]

    def evaluate(potential: NDArray[np.float64]):
        with np.errstate(over="ignore", invalid="ignore"):
            n = np.exp(potential + phi_n)
            p = np.exp(-potential - phi_p)
            state = State(u=potential, n=n, p=p)
        f = poisson_residual_nd(
            device.mesh.x_nodes_nd, state.u, state.n, state.p, device.cbar_nd, targets[0], targets[3]
        )
        scale = np.maximum.reduce((np.ones(m - 2), np.abs(device.cbar_nd), state.n[1:-1], state.p[1:-1]))
        return state, f, float(np.max(np.abs(f[1:-1]) / scale))

    for _ in range(max_iter):
        state, f, residual = evaluate(u)
        if residual <= 1e-10:
            return u
        diagonal = -(1.0 / h[:-1] + 1.0 / h[1:]) / w - state.n[1:-1] - state.p[1:-1]
        left = np.empty(m - 2, dtype=np.float64)
        right = np.empty(m - 2, dtype=np.float64)
        left[1:] = 1.0 / (w[1:] * h[1:-1])
        right[:-1] = 1.0 / (w[:-1] * h[1:-1])
        step = _solve_tridiagonal(diagonal, left, right, -f[1:-1])

        alpha = 1.0
        accepted = False
        while alpha >= 1e-6:
            trial = u.copy()
            trial[1:-1] += alpha * step
            try:
                _, _, trial_residual = evaluate(trial)
                if trial_residual < residual or trial_residual <= 1e-10:
                    u = trial
                    accepted = True
                    break
            except _EXPECTED_FLOAT_ERRORS:
                pass
            alpha *= 0.5
        if not accepted:
            raise ArithmeticError("POISSON_LINE_SEARCH_FAILED")
    raise ArithmeticError("POISSON_MAX_INNER")


def _qf_carrier_step(
    device: NormalizedDevice,
    u: NDArray[np.float64],
    eta_n: NDArray[np.float64],
    eta_p: NDArray[np.float64],
    *,
    electron: bool,
    ref_n: float = 0.0,
    ref_p: float = 0.0,
    max_iter: int = 30,
    tol: float = 1e-11,
) -> NDArray[np.float64]:
    """固定 `u` 与另一个载流子的 QF，解 `eta_n`（`electron=True`）或 `eta_p`。

    SRH 分母 `tau_p(n+1) + tau_n(p+1)` 在本子步开始时冻结：Gummel 的线性化只保留
    主未知量的一阶依赖，把分母当作常数，才能保证每条边只连两个未知量、Jacobian 保持三对角。
    分母冻结带来的误差由外层迭代消除。

    内层同样是阻尼 Newton：`delta` 由三对角求解得到，`alpha` 从 1 折半直到 scaled residual 下降。
    """
    h = device.mesh.h_edges_nd
    w = device.mesh.w_interior_nd
    du = np.diff(u)
    bp = bernoulli(du)
    bn = bernoulli(-du)

    with np.errstate(over="ignore", invalid="ignore"):
        n_lag = np.exp((u + ref_n) + eta_n)
        p_lag = np.exp((-u - ref_p) - eta_p)
    denominator = device.tau_n_nd * (p_lag[1:-1] + 1.0) + device.tau_p_nd * (n_lag[1:-1] + 1.0)

    def evaluate(eta: NDArray[np.float64]):
        """返回 `(residual, diag, left, right, scaled_residual)`，供 Newton 与线搜索共用。"""
        phi_n = eta if electron else eta_n
        phi_p = eta_p if electron else eta
        with np.errstate(over="ignore", invalid="ignore"):
            n = np.exp((u + ref_n) + phi_n)
            p = np.exp((-u - ref_p) - phi_p)
            flux = edge_flux_qf(
                u, phi_n, phi_p, h, device.mu_n_nd, device.mu_p_nd, ref_n=ref_n, ref_p=ref_p
            )
        # np - 1 用 expm1 直接算指数差，避免 np≈1 处的大数相减。
        rate = np.expm1((ref_n - ref_p) + (phi_n[1:-1] - phi_p[1:-1])) / denominator

        if electron:
            j = flux.jn_edges
            residual = np.diff(j) / w - rate
            mobility = device.mu_n_nd
            # SG 通量对 n(eta_n) 的导数：j 在边 i 上只含 n_i 与 n_{i+1}。
            left = -(mobility * bn[:-1] * n[:-2]) / (w * h[:-1])
            right = -(mobility * bp[1:] * n[2:]) / (w * h[1:])
            diagonal = (
                mobility * n[1:-1] * (bn[1:] / h[1:] + bp[:-1] / h[:-1]) / w
                + n[1:-1] * p[1:-1] / denominator
            )
        else:
            j = flux.jp_edges
            residual = np.diff(j) / w + rate
            mobility = device.mu_p_nd
            left = -(mobility * bp[:-1] * p[:-2]) / (w * h[:-1])
            right = -(mobility * bn[1:] * p[2:]) / (w * h[1:])
            diagonal = (
                mobility * p[1:-1] * (bp[1:] / h[1:] + bn[:-1] / h[:-1]) / w
                + n[1:-1] * p[1:-1] / denominator
            )

        scale = np.maximum.reduce((np.ones_like(residual), np.abs(j[1:]) / w, np.abs(j[:-1]) / w, np.abs(rate)))
        return residual, diagonal, left, right, float(np.max(np.abs(residual) / scale))

    current = (eta_n if electron else eta_p).copy()
    for _ in range(max_iter):
        residual, diagonal, left, right, rn = evaluate(current)
        if rn < tol:
            return current
        delta = _solve_tridiagonal(diagonal, left, right, residual)
        alpha = 1.0
        accepted = False
        for _ in range(30):
            trial = current.copy()
            trial[1:-1] += alpha * delta
            try:
                *_, trial_rn = evaluate(trial)
                if trial_rn < rn or trial_rn < tol:
                    current = trial
                    accepted = True
                    break
            except _EXPECTED_FLOAT_ERRORS:
                pass
            alpha *= 0.5
        if not accepted:
            name = "ELECTRON" if electron else "HOLE"
            raise ArithmeticError(f"{name}_INNER_LINE_SEARCH_FAILED at {rn:.3g}")
    name = "ELECTRON" if electron else "HOLE"
    raise ArithmeticError(f"{name}_INNER_MAX_ITER rn={rn:.3g}")


def _project_contacts(qf: QFState, targets: tuple) -> QFState:
    """把 6 个 Dirichlet 值写回两端结点；内部结点保持原值。"""
    u = qf.u.copy()
    eta_n = qf.eta_n.copy()
    eta_p = qf.eta_p.copy()
    u[0], u[-1] = targets[0], targets[3]
    eta_n[0], eta_n[-1] = targets[1], targets[4]
    eta_p[0], eta_p[-1] = targets[2], targets[5]
    return QFState(u=u, eta_n=eta_n, eta_p=eta_p, ref_n=qf.ref_n, ref_p=qf.ref_p)


def _interpolate_qf(old: QFState, target: QFState, alpha: float) -> QFState:
    """外层线搜索的试探点：`old + alpha * (target - old)`。

    `ref_n/ref_p` 是整次求解固定的锚点，插值中改变会让 `n = exp(u+ref_n+eta_n)` 不再对应
    同一个物理状态，所以这里显式检查它们没有变。
    """
    if old.ref_n != target.ref_n or old.ref_p != target.ref_p:
        raise ValueError("QF reference changed during iteration")
    return QFState(
        u=old.u + alpha * (target.u - old.u),
        eta_n=old.eta_n + alpha * (target.eta_n - old.eta_n),
        eta_p=old.eta_p + alpha * (target.eta_p - old.eta_p),
        ref_n=old.ref_n,
        ref_p=old.ref_p,
    )


def _physics_from_evaluation(evaluation: ResidualEvaluation) -> PhysicsDiagnostics:
    """把最终残差评估裁成对外的诊断包；QF 迭代变量本身不暴露。"""
    return PhysicsDiagnostics(
        fluxes=evaluation.fluxes,
        recombination_nd=evaluation.recombination_nd,
        residual_qf_nd=evaluation.f,
    )


def solve_gummel(
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

    `diagnostics`（可选）在成功时写入内部 QF 状态与最终通量/残差，供回归测试
    核对“求解器内部量”与“公开 `evaluate_qf` 重算量”逐位一致；失败时不写入。
    """
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
            or (name != "bias_V" and value <= 0)
        ):
            raise ValueError(f"invalid {name}")
    if alpha_min > 1.0:
        raise ValueError("alpha_min must be <= 1")
    if isinstance(max_iter, bool) or not isinstance(max_iter, (int, np.integer)) or max_iter < 1:
        raise ValueError("max_iter must be a positive integer")
    if diagnostics is not None:
        diagnostics.clear()

    targets = device.contact_targets(bias_V)  # (uL, nL, pL, uR, nR, pR)
    # 锚点取右端电子与左端空穴：n_R = exp(u_R + ref_n)、p_L = exp(-u_L - ref_p)，
    # 使两端接触的 eta_n（或 eta_p）中至少有一端为 0。锚点本身不影响残差
    # （见 tests/test_oracle_qf_residual.py 的平移不变性测试），只影响 exp 的取整。
    ref_n = float(np.log(targets[4]) - targets[3])
    ref_p = float(-np.log(targets[2]) - targets[0])
    targets_qf = contact_targets_qf(device=device, bias_V=bias_V, ref_n=ref_n, ref_p=ref_p)

    history: list[dict] = []

    def make_result(
        ok: bool,
        iterations: int,
        rn: float,
        dn: float,
        reason: str,
        qf: QFState | None = None,
        evaluation: ResidualEvaluation | None = None,
        physical: State | None = None,
    ) -> SolveResult:
        physics = None
        if ok:
            if qf is None or evaluation is None or physical is None:
                raise RuntimeError("successful solve is missing the final state")
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
        qf = _project_contacts(qf, targets_qf)
        physical = qf.to_physical()
        evaluation = evaluate_qf(state=qf, device=device, bias_V=bias_V)
        rn = scaled_residual_qf(qf, device, bias_V, evaluation, physical)
    except _EXPECTED_FLOAT_ERRORS as exc:
        history.append({"iteration": 0, "status": "INVALID_INITIAL_STATE", "reason": str(exc)})
        return make_result(False, 0, np.inf, np.inf, "INVALID_INITIAL_STATE")

    if rn <= tol_residual:
        return make_result(True, 0, rn, 0.0, "SUCCESS", qf=qf, evaluation=evaluation, physical=physical)

    dn = np.inf
    for iteration in range(1, max_iter + 1):
        try:
            # Poisson 子步收的是总 QF 势，所以先加回锚点：phi_n = eta_n + ref_n。
            u_star = _poisson_step(
                device, qf.u, qf.eta_n + qf.ref_n, qf.eta_p + qf.ref_p, targets
            )
            eta_n_star = _qf_carrier_step(
                device,
                u_star,
                qf.eta_n,
                qf.eta_p,
                electron=True,
                ref_n=qf.ref_n,
                ref_p=qf.ref_p,
            )
            eta_p_star = _qf_carrier_step(
                device,
                u_star,
                eta_n_star,
                qf.eta_p,
                electron=False,
                ref_n=qf.ref_n,
                ref_p=qf.ref_p,
            )
            star = _project_contacts(
                QFState(
                    u=u_star,
                    eta_n=eta_n_star,
                    eta_p=eta_p_star,
                    ref_n=qf.ref_n,
                    ref_p=qf.ref_p,
                ),
                targets_qf,
            )
        except _EXPECTED_FLOAT_ERRORS as exc:
            history.append({"iteration": iteration, "status": "SUBSTEP_FAILED", "reason": str(exc)})
            return make_result(False, iteration, rn, dn, "SUBSTEP_FAILED: " + str(exc))

        alpha = 1.0
        accepted = False
        while alpha >= alpha_min:
            trial_rn = np.inf
            trial_dn = np.inf
            trial_error = None
            try:
                trial_qf = _interpolate_qf(qf, star, alpha)
                trial_physical = trial_qf.to_physical()
                trial_evaluation = evaluate_qf(trial_qf, device, bias_V)
                trial_rn = scaled_residual_qf(trial_qf, device, bias_V, trial_evaluation, trial_physical)
                trial_dn = scaled_update(physical, trial_physical, n_floor)
                good = trial_rn < rn or trial_rn <= tol_residual
            except _EXPECTED_FLOAT_ERRORS as exc:
                good = False
                trial_error = str(exc)
            history.append(
                {
                    "iteration": iteration,
                    "alpha": float(alpha),
                    "residual_norm": float(trial_rn),
                    "update_norm": float(trial_dn),
                    "status": "ACCEPTED" if good else "REJECTED",
                    "reason": trial_error,
                }
            )
            if good:
                qf = trial_qf
                physical = trial_physical
                evaluation = trial_evaluation
                rn = trial_rn
                dn = trial_dn
                accepted = True
                break
            alpha *= 0.5
        if not accepted:
            return make_result(False, iteration, rn, dn, "LINE_SEARCH_FAILED")
        if is_converged(rn, dn, tol_residual, tol_update):
            return make_result(
                True, iteration, rn, dn, "SUCCESS", qf=qf, evaluation=evaluation, physical=physical
            )
    return make_result(False, max_iter, rn, dn, "MAX_ITER")
