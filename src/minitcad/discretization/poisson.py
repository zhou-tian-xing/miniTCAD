#Poisson方程的离散残差
import numpy as np
from numpy.typing import NDArray,ArrayLike
def _vector(
        name:str,
        data:ArrayLike,
        size:int | None = None
) -> NDArray[np.float64]:
    arr = np.asarray(data,dtype = np.float64)
    if (arr.ndim != 1 or (
        size is not None and arr.size != size
    )):
        raise ValueError(
            f"{name} must be 1D vector of length {size}"
        )
    if not np.isfinite(arr).all():
        raise ValueError(
            f"{name} must contain only finite values"
        )
    return arr
def poisson_residual_nd(
        xi:ArrayLike,
        u:ArrayLike,
        n:ArrayLike,
        p:ArrayLike,
        cbar:ArrayLike,#C=ND​−NA​
        u_left:float,
        u_right:float
) -> NDArray[np.float64]:
    x = _vector("xi",xi)
    m = x.size
    if m < 3:
        raise ValueError(
            "At least three mesh nodes are required"
        )
    potential = _vector("u",u,m)
    electrons = _vector("n",n,m)
    holes = _vector("p",p,m)
    doping = _vector("Cbar",cbar)
    if doping.size == m:
        doping_interior = doping[1:-1]
    elif doping.size == m - 2:
        doping_interior = doping
    else:
        raise ValueError("Cbar must have length M-2 (CV) or M (Week02 test format)")
    if (np.any(holes <= 0.0)) or (np.any(electrons <= 0.0)):
        raise ValueError(
            "n and p must be strictly positive"
        )
    h = np.diff(x)
    if np.any(h <= 0):
        raise ValueError(
            "xi must be strictly increasing"
        )
    w = (h[:-1] + h[1:]) / 2.0
    gradients = np.diff(potential) / h
    target = np.array(
        [u_left,u_right],
        dtype = np.float64
    )
    if not np.isfinite(target).all():
        raise ValueError(
            "Dirichlet target potentials must be finite"
        )
    out = np.empty(m,dtype = np.float64)
    out[0] = potential[0] - target[0]
    out[-1] = potential[-1] - target[-1]
    out[1:-1] = (
        np.diff(gradients) / w
        + doping_interior
        + holes[1:-1]
        - electrons[1:-1]
    )
    if not np.isfinite(out).all():
        raise FloatingPointError(
            "Poisson residual is nonfinite; "
            "check grid/scaling"
        )
    return out
poisson_discretization = poisson_residual_nd