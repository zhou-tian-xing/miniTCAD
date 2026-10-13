"""W2-B1：``Mesh1D`` 与两个生成器（第 5.1、5.3 节）。"""

from __future__ import annotations

import numpy as np
import pytest

from minitcad.core.mesh import Mesh1D, segmented, uniform

L0_CM = 4.08845452e-3  # 附录 C


def test_uniform_shapes_and_strict_increase():
    mesh = uniform(1.0e-3, 201, L0_CM)
    assert mesh.node_count == 201
    assert mesh.h_edges_cm.size == 200
    assert mesh.w_interior_cm.size == 199
    assert np.all(np.diff(mesh.x_nodes_cm) > 0.0)
    assert mesh.x_nodes_cm[0] == 0.0
    assert mesh.x_nodes_cm[-1] == pytest.approx(1.0e-3)


def test_interior_control_volume_width_matches_two_sided_spacing():
    # 非均匀网格：w_i 必须取 (h_{i-1}+h_i)/2，不能错取下一条边。
    mesh = segmented(1.5, (0.0, 0.25, 0.75, 1.5), (2, 3, 3), L0_CM)
    h = mesh.h_edges_cm
    expected = 0.5 * (h[:-1] + h[1:])
    assert np.allclose(mesh.w_interior_cm, expected, rtol=0.0, atol=0.0)


def test_segmented_keeps_junction_on_a_node():
    # 结区 0.25 必须作为结点出现（第 5.3 节）。
    mesh = segmented(1.5, (0.0, 0.25, 0.75, 1.5), (2, 3, 3), L0_CM)
    assert np.isclose(mesh.x_nodes_cm, 0.25).any()


def test_cm_and_nd_share_one_scale():
    mesh = uniform(1.0e-3, 51, L0_CM)
    assert np.allclose(mesh.x_nodes_nd, mesh.x_nodes_cm / L0_CM, rtol=0.0, atol=0.0)
    assert mesh.normalization_scale() == pytest.approx(L0_CM, rel=1e-12)


def test_nd_round_trip_recovers_cm_grid():
    mesh = uniform(1.0e-3, 33, L0_CM)
    assert np.allclose(mesh.x_nodes_nd * L0_CM, mesh.x_nodes_cm, rtol=1e-12, atol=0.0)


def test_too_few_nodes_is_rejected():
    with pytest.raises(ValueError):
        uniform(1.0e-3, 2, L0_CM)


def test_nonpositive_length_is_rejected():
    with pytest.raises(ValueError):
        uniform(0.0, 10, L0_CM)


def test_breakpoints_must_cover_the_device():
    with pytest.raises(ValueError):
        segmented(1.5, (0.0, 0.5), (3,), L0_CM)  # 未覆盖到 1.5
    with pytest.raises(ValueError):
        segmented(1.5, (0.0, 0.5, 1.0), (3, 3), L0_CM)  # 末点不是 1.5


def test_non_monotone_breakpoints_are_rejected():
    with pytest.raises(ValueError):
        segmented(1.5, (0.0, 0.9, 0.5, 1.5), (2, 2, 2), L0_CM)


def test_mesh_arrays_are_read_only():
    mesh = uniform(1.0e-3, 11, L0_CM)
    with pytest.raises(ValueError):
        mesh.x_nodes_cm[0] = 1.0  # 冻结数组，禁止原地改写（第 5.3 节）


def test_inconsistent_nd_mesh_is_rejected():
    mesh = uniform(1.0e-3, 11, L0_CM)
    with pytest.raises(ValueError):
        Mesh1D(
            x_nodes_cm=mesh.x_nodes_cm,
            h_edges_cm=mesh.h_edges_cm,
            w_interior_cm=mesh.w_interior_cm,
            x_nodes_nd=mesh.x_nodes_nd * 2.0,  # 与 cm 网格不成单一尺度
            h_edges_nd=mesh.h_edges_nd,
            w_interior_nd=mesh.w_interior_nd,
        )
