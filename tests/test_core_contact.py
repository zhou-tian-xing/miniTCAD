"""W2-B3 / T03 / T04：接触浓度、接触势与极性（第 4.1、4.2、11.2 节）。"""

from __future__ import annotations

import math

import pytest

from minitcad.core.contact import (
    INVALID_PN_ORIENTATION,
    InvalidPnOrientation,
    contact_concentrations_nd,
    contact_potential_nd,
    contact_targets,
    validate_pn_orientation,
)
from minitcad.core.device import (
    DeviceConfig,
    MeshSpec,
    NormalizationConfig,
    build_device,
)
from minitcad.core.material import SILICON
from minitcad.core.region import DopingProfile, Region

LENGTH_CM = 1.0e-3
JUNCTION_CM = 2.5e-4


class _DeviceStub:
    """独立的最小接触被测对象：不经过 B 的 DeviceConfig 构建路径。"""

    def __init__(self, c_left_nd: float, c_right_nd: float, vt_v: float) -> None:
        self.net_doping_left_nd = c_left_nd
        self.net_doping_right_nd = c_right_nd
        self.thermal_voltage_v = vt_v


def _example_device():
    config = DeviceConfig(
        length_cm=LENGTH_CM,
        area_cm2=1.0e-4,
        temperature_k=300.0,
        material=SILICON,
        regions=(
            Region(0.0, JUNCTION_CM, "silicon", acceptor_cm3=1.0e16, donor_cm3=0.0),
            Region(JUNCTION_CM, LENGTH_CM, "silicon", acceptor_cm3=0.0, donor_cm3=2.0e15),
        ),
        mesh=MeshSpec(kind="uniform", nodes=51),
        normalization=NormalizationConfig(reference_diffusivity_cm2_s=10.0),
    )
    return build_device(config)


@pytest.mark.parametrize("c_nd", [-1.0e6, -100.0, 100.0, 2.0e5])
def test_mass_action_and_neutrality(c_nd):
    # T03：np = 1 且 p - n + C = 0，误差 <= 1e-10。
    n, p = contact_concentrations_nd(c_nd)
    assert abs(n * p - 1.0) <= 1.0e-10
    assert abs(p - n + c_nd) <= 1.0e-10 * max(1.0, abs(c_nd))


def test_contact_potential_signs_follow_doping():
    # 左 P 要求 a(C_L) < 0，右 N 要求 a(C_R) > 0。
    assert contact_potential_nd(-1.0e6) < 0.0
    assert contact_potential_nd(2.0e5) > 0.0
    assert contact_potential_nd(0.0) == 0.0


def test_contact_targets_layout_and_zero_bias():
    stub = _DeviceStub(-1.0e6, 2.0e5, 0.0258519998)
    u_l, n_l, p_l, u_r, n_r, p_r = contact_targets(stub, 0.0)
    assert u_l == pytest.approx(contact_potential_nd(-1.0e6))
    assert u_r == pytest.approx(contact_potential_nd(2.0e5))
    # 右端为 0 V 参考并不意味着 u(l) = 0（第 4.2 节）。
    assert u_r != 0.0
    assert n_l == pytest.approx(math.exp(u_l))
    assert p_r == pytest.approx(math.exp(-u_r))


def test_forward_bias_lowers_the_barrier_by_exactly_v():
    # T04：V > 0 时 u_R - u_L 减少 v，偏差 <= 1e-10。
    vt = 0.0258519998
    stub = _DeviceStub(-1.0e6, 2.0e5, vt)
    barrier_0 = contact_targets(stub, 0.0)
    barrier_f = contact_targets(stub, 0.2)
    delta = (barrier_0[3] - barrier_0[0]) - (barrier_f[3] - barrier_f[0])
    assert delta == pytest.approx(0.2 / vt, rel=1e-12)


def test_bias_only_shifts_left_potential_never_concentrations():
    stub = _DeviceStub(-1.0e6, 2.0e5, 0.0258519998)
    at_zero = contact_targets(stub, 0.0)
    at_forward = contact_targets(stub, 0.3)
    assert at_forward[0] - at_zero[0] == pytest.approx(0.3 / 0.0258519998)
    assert at_forward[1:4] == pytest.approx(at_zero[1:4])
    assert at_forward[4:] == pytest.approx(at_zero[4:])


def test_reversed_pn_orientation_is_rejected_with_values():
    profile = DopingProfile(
        (
            Region(0.0, JUNCTION_CM, "silicon", acceptor_cm3=0.0, donor_cm3=1.0e16),
            Region(JUNCTION_CM, LENGTH_CM, "silicon", acceptor_cm3=2.0e15, donor_cm3=0.0),
        ),
        LENGTH_CM,
    )
    with pytest.raises(InvalidPnOrientation) as excinfo:
        validate_pn_orientation(profile)
    assert excinfo.value.code == INVALID_PN_ORIENTATION
    assert "C_left" in str(excinfo.value) and "C_right" in str(excinfo.value)


def test_pn_config_with_reversed_regions_is_rejected():
    with pytest.raises(InvalidPnOrientation):
        DeviceConfig(
            length_cm=LENGTH_CM,
            area_cm2=1.0e-4,
            temperature_k=300.0,
            material=SILICON,
            regions=(
                Region(0.0, JUNCTION_CM, "silicon", acceptor_cm3=0.0, donor_cm3=1.0e16),
                Region(JUNCTION_CM, LENGTH_CM, "silicon", acceptor_cm3=2.0e15, donor_cm3=0.0),
            ),
            mesh=MeshSpec(kind="uniform", nodes=11),
            normalization=NormalizationConfig(reference_diffusivity_cm2_s=10.0),
        )


def test_device_contact_targets_match_stub():
    device = _example_device()
    from_device = device.contact_targets(0.15)
    stub = _DeviceStub(
        device.net_doping_left_nd, device.net_doping_right_nd, device.thermal_voltage_v
    )
    assert from_device == pytest.approx(contact_targets(stub, 0.15))
