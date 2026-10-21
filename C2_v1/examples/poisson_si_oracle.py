"""独立 SI Poisson 教学验证器：内部 C/m³，Dirichlet 边界 V。

请勿直接将其残差范数用于生产无量纲求解。
"""
import numpy as np

def poisson_residual_si_mms(
    x_m: np.ndarray,
    phi_V: np.ndarray,
    rho_C_m3: np.ndarray,
    epsilon_F_m: float,
    phi_left_V: float,
    phi_right_V: float,
) -> np.ndarray:
    """
    独立 SI 有限体积验证器。
    内部残差单位:C/m^3
    边界残差单位:V
    不复用生产无量纲 Poisson 残差。
    """
    x = np.asarray(x_m, dtype=np.float64)
    phi = np.asarray(phi_V, dtype=np.float64)
    rho = np.asarray(rho_C_m3, dtype=np.float64)
    if (
        x.ndim != 1
        or x.size < 3
        or phi.shape != x.shape
        or rho.shape != x.shape
    ):
        raise ValueError(
            "x, phi, rho 必须是等长、至少 3 节点的一维数组"
        )
    if not (
        np.isfinite(x).all()
        and np.isfinite(phi).all()
        and np.isfinite(rho).all()
    ):
        raise ValueError("x, phi, rho 必须有限")
    if np.any(np.diff(x) <= 0.0):
        raise ValueError("x 必须严格递增")
    epsilon = float(epsilon_F_m)
    if not (np.isfinite(epsilon) and epsilon > 0.0):
        raise ValueError(
            "epsilon_F_m 必须是正且有限的绝对介电常数"
        )
    bc = np.asarray(
        [phi_left_V, phi_right_V],
        dtype=np.float64,
    )
    if not np.isfinite(bc).all():
        raise ValueError("边界电势必须有限")
    out = np.empty(x.size, dtype=np.float64)
    out[0] = phi[0] - bc[0]
    out[-1] = phi[-1] - bc[-1]
    for i in range(1, x.size - 1):
        h_left = x[i] - x[i - 1]
        h_right = x[i + 1] - x[i]
        cv_width = (h_left + h_right) / 2.0
        flux_left = (
            epsilon
            * (phi[i] - phi[i - 1])
            / h_left
        )
        flux_right = (
            epsilon
            * (phi[i + 1] - phi[i])
            / h_right
        )
        out[i] = (
            -(flux_right - flux_left) / cv_width
            - rho[i]
        )
    if not np.isfinite(out).all():
        raise FloatingPointError("SI Poisson 残差非有限")
    return out
