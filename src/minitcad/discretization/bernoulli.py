#伯努利函数
import numpy as np
from numpy.typing import NDArray,ArrayLike
def bernoulli(z : ArrayLike) -> float | NDArray[np.float64] :
    x = np.asarray(z,dtype = np.float64)
    if not np.all(np.isfinite(x)):
        raise ValueError("Bernoulli input must be finite")
    out = np.empty_like(x,dtype = np.float64)
    small = np.abs(x) < 1e-3
    pos = x > 50.0
    neg = x < -50.0
    med = ~ (small | pos | neg)
    t1 = x[small]
    out[small] = 1 - t1 / 2 + t1 ** 2 / 12 - t1 ** 4 / 720 + t1 ** 6 / 30240
    t2 = x[pos]
    exp_half = np.exp(-t2 / 2)
    out[pos] = (t2 * exp_half) * exp_half / (1 - exp_half ** 2)
    t3 = x[neg]
    out[neg] = - t3 / (1 - np.exp(t3))
    t4 = x[med]
    out[med] = t4 / np.expm1(t4)
    if out.ndim == 0:
        return float(out)
    else:
        return out