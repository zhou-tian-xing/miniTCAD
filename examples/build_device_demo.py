"""端到端演示：YAML -> DeviceConfig -> NormalizedDevice -> 接触目标。

运行（在仓库根目录）：

    python examples/build_device_demo.py
"""

from __future__ import annotations

import pathlib
import sys

_SRC = pathlib.Path(__file__).resolve().parents[1] / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from minitcad.core.device import build_device  # noqa: E402
from minitcad.io.config import load_config  # noqa: E402

YAML = pathlib.Path(__file__).resolve().parent / "pn_left_p_right_n.yaml"


def main() -> None:
    config = load_config(YAML)
    device = build_device(config)

    print("PN orientation: P(left, x=0) -> N(right, x=L)")
    print("V_app = V_left - V_right; +V = forward bias; right = 0 V")
    print()
    print(f"input unit system : cm (declared in {YAML.name})")
    print(f"nodes             : {device.mesh.node_count}")
    print(f"T                 : {device.temperature_k} K")
    print(f"V_T               : {device.thermal_voltage_v:.10g} V")
    print(f"L_0               : {device.l0_cm:.10g} cm")
    print(f"D_0               : {device.d0_cm2_s:.10g} cm^2/s")
    print(f"t_0               : {device.scales.t0_s:.10g} s")
    print(f"J_0               : {device.scales.j0_a_cm2:.10g} A/cm^2")
    print(f"C_left / C_right  : {device.net_doping_left_nd:.6g} / "
          f"{device.net_doping_right_nd:.6g}  (dimensionless)")
    print(f"material policy   : {device.material_parameter_policy}")
    print()

    for bias in (0.0, 0.2, 0.6):
        u_l, _n_l, _p_l, u_r, _n_r, _p_r = device.contact_targets(bias)
        barrier = u_r - u_l
        print(f"V_app = {bias:+.2f} V -> contact barrier u_R - u_L = {barrier:.6f}"
              f"  (V_bi = {device.thermal_voltage_v * (u_r - (u_l - bias / device.thermal_voltage_v)):.4f} V)")


if __name__ == "__main__":
    main()
