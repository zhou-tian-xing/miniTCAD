"""W2-B2：分段掺杂与控制体平均 ``Cbar_i``（第 5.1 节）。"""

from __future__ import annotations

import numpy as np
import pytest

from minitcad.core.mesh import uniform
from minitcad.core.region import DopingProfile, Region

L0_CM = 4.08845452e-3


def _pn_profile(length_cm=1.0e-3, junction_cm=2.5e-4) -> DopingProfile:
    return DopingProfile(
        (
            Region(0.0, junction_cm, "silicon", acceptor_cm3=1.0e16, donor_cm3=0.0),
            Region(junction_cm, length_cm, "silicon", acceptor_cm3=0.0, donor_cm3=2.0e15),
        ),
        length_cm,
    )


def test_net_doping_sign_and_values():
    profile = _pn_profile()
    assert profile.net_doping_left_cm3 == pytest.approx(-1.0e16)
    assert profile.net_doping_right_cm3 == pytest.approx(2.0e15)
    assert profile.is_left_p_right_n


def test_cv_average_equals_hand_computed_step():
    # 结点恰在结区上：左半控制体全在 P 区，右半全在 N 区。
    length_cm, junction = 1.0e-3, 5.0e-4
    profile = _pn_profile(length_cm, junction)
    mesh = uniform(length_cm, 3, L0_CM)
    node = mesh.x_nodes_cm[1]  # 0.5e-3 = 结区
    assert node == pytest.approx(junction)
    cbar = profile.cv_averaged_net_doping_cm3(mesh.x_nodes_cm, mesh.h_edges_cm)
    # (C_left*w + C_right*w)/2w = (C_left + C_right)/2
    expected = 0.5 * (profile.net_doping_left_cm3 + profile.net_doping_right_cm3)
    assert cbar[0] == pytest.approx(expected)


def test_cv_average_over_partially_covering_control_volume():
    length_cm, junction = 1.0e-3, 2.5e-4
    profile = _pn_profile(length_cm, junction)
    x = np.array([0.0, 2.0e-4, 4.0e-4, length_cm])
    h = np.diff(x)
    cbar = profile.cv_averaged_net_doping_cm3(x, h)
    # 结点 2 (x=4e-4) 的控制体是 [3e-4, 7e-4]，右侧全在 N 区，落在 N 区应为 +2e15。
    assert cbar[1] == pytest.approx(2.0e15)
    # 结点 1 (x=2e-4) 的控制体是 [1e-4, 3e-4]，跨结区：0.15e-3*C_left + 0.05e-3*C_right 除以 0.2e-3
    expected = (1.5e-4 * (-1.0e16) + 0.5e-4 * (2.0e15)) / 2.0e-4
    assert cbar[0] == pytest.approx(expected)


def test_gap_in_regions_is_rejected():
    with pytest.raises(ValueError):
        DopingProfile(
            (
                Region(0.0, 2.0e-4, "silicon", acceptor_cm3=1.0e16, donor_cm3=0.0),
                Region(3.0e-4, 1.0e-3, "silicon", acceptor_cm3=0.0, donor_cm3=2.0e15),
            ),
            1.0e-3,
        )


def test_overlap_in_regions_is_rejected():
    with pytest.raises(ValueError):
        DopingProfile(
            (
                Region(0.0, 5.0e-4, "silicon", acceptor_cm3=1.0e16, donor_cm3=0.0),
                Region(4.0e-4, 1.0e-3, "silicon", acceptor_cm3=0.0, donor_cm3=2.0e15),
            ),
            1.0e-3,
        )


def test_uncovered_device_is_rejected():
    with pytest.raises(ValueError):
        DopingProfile(
            (Region(0.0, 9.0e-4, "silicon", acceptor_cm3=1.0e16, donor_cm3=0.0),),
            1.0e-3,
        )


def test_negative_doping_is_rejected():
    with pytest.raises(ValueError):
        Region(0.0, 1.0e-4, "silicon", acceptor_cm3=-1.0, donor_cm3=0.0)


def test_degenerate_region_bounds_are_rejected():
    with pytest.raises(ValueError):
        Region(1.0e-4, 1.0e-4, "silicon", acceptor_cm3=1.0, donor_cm3=0.0)


def test_reversed_orientation_is_detected():
    reversed_profile = DopingProfile(
        (
            Region(0.0, 5.0e-4, "silicon", acceptor_cm3=0.0, donor_cm3=1.0e16),
            Region(5.0e-4, 1.0e-3, "silicon", acceptor_cm3=2.0e15, donor_cm3=0.0),
        ),
        1.0e-3,
    )
    assert not reversed_profile.is_left_p_right_n
