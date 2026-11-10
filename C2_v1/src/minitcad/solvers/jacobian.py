#给牛顿法提供雅可比矩阵
#九色稀疏差分法，2026/10/24
import numpy as np
from scipy.sparse import coo_matrix, csc_matrix
from minitcad.core.state import QFState
from minitcad.core.device import NormalizedDevice
from minitcad.discretization.residual import residual_qf
def jacobian_qf(
    state : QFState,
    device : NormalizedDevice,
    bias_V : float,
    *,
    rel_step : float = 1e-5
) -> csc_matrix:
    m = state.node_count
    if m != device.mesh.node_count:
        raise ValueError("QFState/device node_count mismatch")
    if isinstance(rel_step, (bool, np.bool_)) or not np.isfinite(rel_step) or rel_step <= 0:
        raise ValueError("rel_step must be positive and finite")
    z = state.pack_qf()
    nvar = z.size
    f0 = residual_qf(state = state,device = device,bias_V = bias_V)
    if f0.shape != (nvar,) or not np.isfinite(f0).all():
        raise FloatingPointError("invalid baseline QF residual")
    steps = rel_step * np.maximum(1.0, np.abs(z))
    if not np.isfinite(steps).all() or np.any(steps <= 0):
        raise ValueError("invalid Jacobian steps")
    rows = []
    cols = []
    data = []
    def evaluate(vector):
        qf = QFState.unpack_qf(z = vector,node_count = m,ref_n = state.ref_n,ref_p = state.ref_p)
        out = residual_qf(state = qf,device = device,bias_V = bias_V)
        if out.shape != (nvar,) or not np.isfinite(out).all():
            raise FloatingPointError("nonfinite perturbed residual")
        return out
    def local_derivative(column,step):
        for _ in range(12):
            if z[column] + step == z[column] or z[column] - step == z[column]:
                break
            zp = z.copy()
            zm = z.copy()
            zp[column] += step
            zm[column] -= step
            plus = None
            minus = None
            try:
                plus = evaluate(vector = zp)
            except (ValueError,FloatingPointError,OverflowError,ArithmeticError):
                pass
            try:
                minus = evaluate(vector = zm)
            except (ValueError,FloatingPointError,OverflowError,ArithmeticError):
                pass
            if plus is not None and minus is not None:
                return (plus - minus) / (2.0 * step)
            if plus is not None:
                return (plus - f0) / step
            if minus is not None:
                return (f0 - minus) / step
            step *= 0.5
        raise FloatingPointError(f"Cannot differentiate QF coordinate {column}")
    for component in range(3):
        for node_color in range(3):
            columns = 3 * np.arange(node_color,m,3) + component
            if columns.size == 0:
                continue
            zp = z.copy()
            zm = z.copy()
            zp[columns] += steps[columns]
            zm[columns] -= steps[columns]
            try:
                fp = evaluate(zp)
                fm = evaluate(zm)
            except (ValueError,FloatingPointError,OverflowError,ArithmeticError):
                fp = None
                fm = None
            for col in columns:
                k = col // 3
                start = 3 * max(0,k - 1)
                end = 3 * min(m,k + 2)
                if fp is not None and fm is not None:
                    local = (fp[start:end] - fm[start:end]) / (2.0 * steps[col])
                else:
                    local = local_derivative(col,steps[col])[start:end]
                if not np.isfinite(local).all():
                    raise FloatingPointError(f"nonfinite Jacobian column {col}")
                for row,value in zip(
                    range(start,end),local
                ):
                    if value != 0.0:
                        rows.append(row)
                        cols.append(int(col))
                        data.append(float(value))
    jac = coo_matrix((data,(rows,cols)),shape = (nvar,nvar)).tocsc()
    if not np.isfinite(jac.data).all():
        raise FloatingPointError("nonfinite sparse QF Jacobian")
    return jac