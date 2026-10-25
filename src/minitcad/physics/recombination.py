#中能级SRH复合（无量纲）,n1=p1=ni
#利用准fermi能级提高了精度,2026/10/21
import numpy as np
from numpy.typing import ArrayLike,NDArray
def srh_recombination_nd(
        n:ArrayLike,
        p:ArrayLike,
        tau_n_nd:float,
        tau_p_nd:float
) -> float | NDArray[np.float64]:
    electrons = np.asarray(n,dtype = np.float64)
    holes = np.asarray(p,dtype = np.float64)
    if not (np.isfinite(tau_n_nd) and tau_n_nd > 0.0 and
            np.isfinite(tau_p_nd) and tau_p_nd > 0.0):
        raise ValueError("normalized lifetimes must be finite and > 0")
    if not (np.isfinite(electrons).all() and np.isfinite(holes).all()):
        raise ValueError("normalized carrier concentrations must be finite")
    if np.any(electrons <= 0.0) or np.any(holes <= 0.0):
        raise ValueError("normalized n and p must be strictly positive")
    try:
        rate = (electrons * holes - 1.0) / (tau_p_nd * (electrons + 1.0) + tau_n_nd * (holes + 1.0))
    except (FloatingPointError,ValueError) as exc:
        raise ValueError("invalid/overflowing normalized SRH intermediate") from exc
    if not np.isfinite(rate).all():
        raise ValueError("SRH result must be finite")
    if rate.ndim == 0:
        return float(rate)
    else:
        return np.asarray(rate,dtype = np.float64)
def srh_recombination_qf(
    u,
    eta_n,
    eta_p,
    tau_n_nd,
    tau_p_nd,
    *,
    ref_n=0.0,
    ref_p=0.0
):
    """
    R = (np - 1) /
        [tau_p*(n+1) + tau_n*(p+1)]
    其中
    np - 1 = expm1(fn - fp)
    """
    if (
        not np.isfinite(
            [tau_n_nd, tau_p_nd, ref_n, ref_p]
        ).all()
        or tau_n_nd <= 0
        or tau_p_nd <= 0
    ):
        raise ValueError("invalid SRH parameters")
    u = np.asarray(u,dtype=np.float64)
    eta_n = np.asarray(eta_n,dtype=np.float64)
    eta_p = np.asarray(eta_p,dtype=np.float64)
    if (
        u.shape != eta_n.shape
        or u.shape != eta_p.shape
    ):
        raise ValueError("SRH quasi-Fermi shape mismatch")
    if not all(
        np.isfinite(a).all()
        for a in (u, eta_n, eta_p)
    ):
        raise ValueError(
            "nonfinite SRH quasi-Fermi inputs"
        )
    n = np.exp((u + ref_n) + eta_n)
    p = np.exp((-u - ref_p) - eta_p)
    denominator = (
            tau_p_nd * (n + 1.0) + tau_n_nd * (p + 1.0)
        )
    # fn - fp = (ref_n-ref_p) + (eta_n-eta_p)
    num = np.expm1(
        (ref_n - ref_p) + (eta_n - eta_p)
    )
    rate = num / denominator
    if not np.isfinite(rate).all():
        raise FloatingPointError(
            "nonfinite quasi-Fermi SRH rate"
        )
    return (
        float(rate)
        if rate.ndim == 0
        else rate
    )