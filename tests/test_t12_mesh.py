"""T12：结区加密网格与网格收敛表（§15 实验 C）。

两条线分开测，因为它们回答不同的问题：

* **加密序列**（uniform 201/401/801）：结点分布不变、只加密，误差必须单调下降，
  最细两级非零 `I` 的相对差 `≤1%`。
* **分布对比**（同为 201 个结点，uniform vs refined）：总结点数固定、只换分布，
  比较 `I`、势垒与峰值电场的误差和耗时。
"""

from __future__ import annotations

import dataclasses
import pathlib

import numpy as np
import pytest

from minitcad.core.device import MeshSpec, build_device
from minitcad.core.mesh import MIN_NODES, refined, uniform
from minitcad.io.config import load_config
from minitcad.verification.reference import (
    DEFAULT_TOLERANCE,
    convergence_table,
    is_monotone,
    relative_difference,
    sample,
)

EXAMPLES = pathlib.Path(__file__).resolve().parents[1] / "examples"
BASELINE = EXAMPLES / "pn_baseline_0_to_030.yaml"

LENGTH_CM = 1.0e-3
JUNCTION_CM = 2.5e-4
L0_CM = 1.0e-3


@pytest.fixture(scope="module")
def config():
    return load_config(BASELINE)


# --- 生成器本身 -------------------------------------------------------------


@pytest.mark.parametrize("nodes", [7, 101, 201, 401])
def test_refined_has_exactly_the_requested_node_count(nodes):
    """`nodes` 是总数，不是每段的数：T12 要在同结点数下比较分布。"""
    mesh = refined(LENGTH_CM, JUNCTION_CM, nodes, L0_CM, refine=20.0)
    assert mesh.node_count == nodes
    assert mesh.h_edges_cm.sum() == pytest.approx(LENGTH_CM, rel=0, abs=1e-18)


def test_refined_puts_the_junction_on_a_node():
    """结区不在结点上就谈不上“在结区加密”。"""
    for nodes in (201, 401):
        mesh = refined(LENGTH_CM, JUNCTION_CM, nodes, L0_CM, refine=20.0)
        assert np.any(mesh.x_nodes_cm == JUNCTION_CM)


def test_refined_spacing_grows_away_from_the_junction():
    """间距从结区向两端单调变疏，且最疏/最密之比等于 `refine`。"""
    mesh = refined(LENGTH_CM, JUNCTION_CM, 201, L0_CM, refine=20.0)
    junction_index = int(np.argwhere(mesh.x_nodes_cm == JUNCTION_CM)[0, 0])
    left = mesh.h_edges_cm[:junction_index]
    right = mesh.h_edges_cm[junction_index:]
    assert np.all(np.diff(left) < 0.0), "左侧应向 x=0 变疏（数组从 0 到结区递减）"
    assert np.all(np.diff(right) > 0.0), "右侧应向 x=L 变疏"
    assert mesh.h_edges_cm.max() / mesh.h_edges_cm.min() == pytest.approx(20.0, rel=0.05)


def test_refine_one_is_bit_identical_to_uniform():
    """`refine=1` 就是均匀网格；两者必须逐位相同，否则对比组不干净。"""
    graded = refined(LENGTH_CM, JUNCTION_CM, 201, L0_CM, refine=1.0)
    plain = uniform(LENGTH_CM, 201, L0_CM)
    assert np.array_equal(graded.x_nodes_cm, plain.x_nodes_cm)
    assert np.array_equal(graded.h_edges_cm, plain.h_edges_cm)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"refine": 0.5},
        {"refine": float("nan")},
        {"nodes": MIN_NODES - 1},
        {"nodes": 2.5},
        {"junction_cm": 0.0},
        {"junction_cm": LENGTH_CM},
    ],
)
def test_refined_rejects_invalid_parameters(kwargs):
    """非法参数显式失败，不静默修正（T12 的比较组一旦被悄悄改就不是比较组）。"""
    arguments = {
        "length_cm": LENGTH_CM,
        "junction_cm": JUNCTION_CM,
        "nodes": 201,
        "l0_cm": L0_CM,
        "refine": 20.0,
    }
    arguments.update(kwargs)
    with pytest.raises(ValueError):
        refined(**arguments)


def test_refined_mesh_reaches_the_device(config):
    """YAML 侧同样能用：`MeshSpec(kind="refined")` 自动取掺杂剖面的结区。"""
    graded = MeshSpec(kind="refined", nodes=201, refine=20.0)
    device = build_device(dataclasses.replace(config, mesh=graded))
    assert device.mesh.node_count == 201
    assert np.any(device.mesh.x_nodes_cm == JUNCTION_CM)


def test_refined_mesh_spec_needs_a_single_junction(config):
    """多结器件目前不支持：显式失败，而不是把第二个结悄悄忽略。"""
    with pytest.raises(ValueError):
        MeshSpec(kind="refined", nodes=201).build(
            LENGTH_CM, L0_CM, (JUNCTION_CM, 0.75e-3)
        )


# --- T12：加密序列与分布对比 ------------------------------------------------


@pytest.fixture(scope="module")
def ladder(config):
    """uniform 201/401/801，最细一级作参考。"""
    samples = [
        sample(config, MeshSpec(kind="uniform", nodes=nodes), label=f"uniform {nodes}")
        for nodes in (201, 401, 801)
    ]
    assert all(item.converged for item in samples), [item.label for item in samples]
    return samples


def test_current_error_decreases_monotonically_with_refinement(ladder):
    """T12 的核心：只加密、不改分布，误差必须单调下降。"""
    reference = ladder[-1]
    rows = convergence_table(ladder[:-1], reference)
    errors = [row["current_rel_error"] for row in rows]
    assert all(error > 0.0 for error in errors)
    assert is_monotone(errors), errors


def test_two_finest_levels_agree_within_one_percent(ladder):
    """T12 判据：最细两级非零 `I` 的相对差目标 `≤1%`。"""
    coarser, finest = ladder[-2], ladder[-1]
    assert coarser.current_density_A_cm2 != 0.0
    error = relative_difference(coarser.current_density_A_cm2, finest.current_density_A_cm2)
    assert error <= DEFAULT_TOLERANCE, f"{error:.3e} > {DEFAULT_TOLERANCE}"


def test_barrier_is_insensitive_to_the_mesh(ladder):
    """势垒由接触边界决定，网格变化不应改变它——变了就说明边界条件写错了。"""
    barriers = [item.barrier_V for item in ladder]
    assert max(barriers) - min(barriers) <= 1e-6 * abs(barriers[-1])


def test_refinement_at_fixed_node_count_beats_uniform(config, ladder):
    """实验 C：结点数相同，结区加密的电场误差必须显著小于均匀网格。"""
    plain = ladder[0]
    graded = sample(config, MeshSpec(kind="refined", nodes=201, refine=20.0), label="refined 201")
    reference = ladder[-1]
    assert graded.converged

    plain_error = relative_difference(plain.current_density_A_cm2, reference.current_density_A_cm2)
    graded_error = relative_difference(
        graded.current_density_A_cm2, reference.current_density_A_cm2
    )
    assert graded_error < plain_error, (graded_error, plain_error)

    plain_field = relative_difference(plain.peak_field_V_cm, reference.peak_field_V_cm)
    graded_field = relative_difference(graded.peak_field_V_cm, reference.peak_field_V_cm)
    assert graded_field < plain_field, (graded_field, plain_field)


def test_refined_and_uniform_agree_on_the_physics(config, ladder):
    """两种分布的解必须落在同一物理量级上（加密不该改变器件行为）。"""
    graded = sample(config, MeshSpec(kind="refined", nodes=201, refine=20.0), label="refined 201")
    plain = ladder[0]
    assert graded.current_density_A_cm2 == pytest.approx(plain.current_density_A_cm2, rel=1e-2)
    assert graded.barrier_V == pytest.approx(plain.barrier_V, rel=1e-6)
