#连续性方程
import numpy as np
from numpy.typing import NDArray,ArrayLike
from minitcad.core.state import Fluxes
def _as_real_1d(
        name:str,
        values:ArrayLike
) -> NDArray[np.float64]:
    raw = np.asarray(values)
    if (raw.ndim != 1 or raw.dtype.kind not in "iuf"):
        raise ValueError(f"{name} must be 1D numeric array")
    result = np.asarray(raw,dtype = np.float64)
    if (not np.isfinite(result).all()):
        raise ValueError(f"{name} must contains finite numbers only")
    return result
def continuity_residual_nd(
        fluxes:Fluxes,
        recombination_nd:ArrayLike,
        w_interior_nd:ArrayLike
) -> tuple[NDArray[np.float64],NDArray[np.float64]]:
    jn = _as_real_1d("jn_edges",fluxes.jn_edges)
    jp = _as_real_1d("jp_edges",fluxes.jp_edges)
    r = _as_real_1d("recombination",recombination_nd)
    w = _as_real_1d("w_interior",w_interior_nd)
    m = jn.size + 1
    if (m < 3 or jp.shape != (m - 1,) or r.shape != (m,) or w.shape != (m - 2,)):
        raise ValueError(            "shape mismatch: flux edges (M-1), "
            "CV widths (M-2), SRH nodes (M) "
            "must agree; M >= 3"
        )
    if np.any(w <= 0.0):
        raise ValueError("w_interior_nd must be strictly positive")
    try:
        fn = np.diff(jn) / w - r[1:-1]
        fp = np.diff(jp) / w + r[1:-1]
    except FloatingPointError as exc:
        raise FloatingPointError("continuity residual must be finite") from exc
    return (fn,fp)