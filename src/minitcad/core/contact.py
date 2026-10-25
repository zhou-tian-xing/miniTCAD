"""理想欧姆接触。

约定：左边是 P，右边是 N；V_app = V_left - V_right，正的 V_app 是正向；右端接地。
接触处局部平衡，所以 n = e^a、p = e^-a，其中 a = asinh(C/2)，C 是净掺杂（已除以 n_i）。

``contact_targets`` 返回物理六元组；A 的 QF 残差使用 ``contact_targets_qf``，
它把同一组边界换算成 ``(u, eta_n, eta_p)`` 的 Dirichlet 目标。
"""

from __future__ import annotations

import math

from minitcad.core.region import DopingProfile

INVALID_PN_ORIENTATION = "INVALID_PN_ORIENTATION"


class InvalidPnOrientation(ValueError):
    """PN 要求左端 C < 0（P 型）、右端 C > 0（N 型）。"""

    code = INVALID_PN_ORIENTATION


def contact_potential_nd(net_doping_nd: float) -> float:
    """接触处的无量纲电势 a(C) = asinh(C/2)。

    用 asinh 而不是 C + sqrt(C^2 + 4)，高掺杂时后者会因相减丢精度。
    """
    if not math.isfinite(net_doping_nd):
        raise ValueError("net_doping_nd must be finite")
    return math.asinh(net_doping_nd / 2.0)


def contact_concentrations_nd(net_doping_nd: float) -> tuple[float, float]:
    """接触处的无量纲载流子浓度 (n, p)。

    由局部电中性 C + p - n = 0 和质量作用 np = 1 解出：n = e^a，p = e^-a。
    """
    a = contact_potential_nd(net_doping_nd)
    return math.exp(a), math.exp(-a)


def validate_pn_orientation(profile: DopingProfile) -> None:
    """检查是不是左 P 右 N；不是就报错并打印两端的实际 C 值。"""
    if not profile.is_left_p_right_n:
        raise InvalidPnOrientation(
            f"{INVALID_PN_ORIENTATION}: pn requires C_left < 0 < C_right; got "
            f"C_left = {profile.net_doping_left_cm3:.6g} cm^-3 and "
            f"C_right = {profile.net_doping_right_cm3:.6g} cm^-3. "
            "P must be on the left (x=0) and N on the right (x=L)"
        )


def contact_targets(device: object, bias_V: float) -> tuple[float, float, float, float, float, float]:
    """算出一组接触边界值 (u_L, n_L, p_L, u_R, n_R, p_R)，都是无量纲。

    bias_V 单位是 V，内部换成 v = V_app / V_T。
    只有左端电势随偏压移动：u_L = a(C_L) + v；右端接地但 u_R = a(C_R) 不为 0。
    浓度不随偏压变化。
    """
    if not math.isfinite(bias_V):
        raise ValueError("bias_V must be finite")
    v = bias_V / float(device.thermal_voltage_v)

    n_l, p_l = contact_concentrations_nd(device.net_doping_left_nd)
    n_r, p_r = contact_concentrations_nd(device.net_doping_right_nd)
    u_l = contact_potential_nd(device.net_doping_left_nd) + v
    u_r = contact_potential_nd(device.net_doping_right_nd)
    return u_l, n_l, p_l, u_r, n_r, p_r


def contact_targets_qf(
    device: object,
    bias_V: float,
    *,
    ref_n: float = 0.0,
    ref_p: float = 0.0,
) -> tuple[float, float, float, float, float, float]:
    """QF 残差用的接触边界值 ``(u_L, eta_nL, eta_pL, u_R, eta_nR, eta_pR)``。

    由物理接触六元组按 ``n = exp(u + ref_n + eta_n)``、
    ``p = exp(-u - ref_p - eta_p)`` 反解，因此本函数给出的是 QF 版 Dirichlet 目标。
    边界局部平衡时 ``eta_n = eta_p = -v``（左端）和 ``0``（右端），``v = V_app/V_T``。
    """
    if not math.isfinite(ref_n) or not math.isfinite(ref_p):
        raise ValueError("QF reference values must be finite")
    u_l, n_l, p_l, u_r, n_r, p_r = contact_targets(device, bias_V)
    return (
        u_l,
        math.log(n_l) - u_l - ref_n,
        -math.log(p_l) - u_l - ref_p,
        u_r,
        math.log(n_r) - u_r - ref_n,
        -math.log(p_r) - u_r - ref_p,
    )
