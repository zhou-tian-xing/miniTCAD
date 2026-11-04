"""W2/W3 公共 API 冻结测试（§10.2、§14.7、I-01）。

B 维护的接口在 W2 末冻结、W3 起不得无通知改名。这里把“名字与形状”变成断言：
包级 ``__all__`` 能全部解析、关键函数的关键字参数名不变、CSV 表头顺序不变。
只检查契约，不重复检查物理正确性（那由 T01–T07 的 oracle 负责）。
"""

from __future__ import annotations

import dataclasses
import inspect

import pytest

import minitcad
import minitcad.core as core
import minitcad.io as io
from minitcad.discretization import ResidualEvaluation, evaluate_qf, residual_qf
from minitcad.simulation import Attempt, RunResult, solve_equilibrium_scalar


@pytest.mark.parametrize("module", [minitcad, core, io])
def test_dunder_all_resolves(module):
    missing = [name for name in module.__all__ if not hasattr(module, name)]
    assert missing == []


def test_contact_entry_points_keep_their_parameter_names():
    assert list(inspect.signature(core.contact_targets).parameters) == ["device", "bias_V"]
    signature = inspect.signature(core.contact_targets_qf)
    assert list(signature.parameters) == ["device", "bias_V", "ref_n", "ref_p"]
    for name in ("ref_n", "ref_p"):
        assert signature.parameters[name].kind is inspect.Parameter.KEYWORD_ONLY
        assert signature.parameters[name].default == 0.0


def test_device_and_state_shape_contracts():
    mesh_fields = [field.name for field in dataclasses.fields(core.Mesh1D)]
    assert mesh_fields == [
        "x_nodes_cm",
        "h_edges_cm",
        "w_interior_cm",
        "x_nodes_nd",
        "h_edges_nd",
        "w_interior_nd",
    ]
    assert [field.name for field in dataclasses.fields(core.State)] == ["u", "n", "p"]
    qf_fields = [field.name for field in dataclasses.fields(core.QFState)]
    assert qf_fields == ["u", "eta_n", "eta_p", "ref_n", "ref_p"]
    assert [field.name for field in dataclasses.fields(core.Fluxes)] == ["jn_edges", "jp_edges"]
    assert [field.name for field in dataclasses.fields(core.PhysicsDiagnostics)] == [
        "fluxes",
        "recombination_nd",
        "residual_qf_nd",
    ]


def test_normalized_device_exposes_the_frozen_scales():
    fields = set(core.NormalizedDevice.__dataclass_fields__)
    properties = {
        "thermal_voltage_v", "l0_cm", "d0_cm2_s", "ni_cm3",
        "net_doping_left_nd", "net_doping_right_nd", "contact_targets",
    }
    for name in ("mesh", "scales", "net_doping_nd", "cbar_nd", "mu_n_nd", "mu_p_nd", "tau_n_nd", "tau_p_nd"):
        assert name in fields, f"NormalizedDevice 缺少字段 {name}"
    for name in properties:
        assert hasattr(core.NormalizedDevice, name), f"NormalizedDevice 缺少属性 {name}"


def test_qf_residual_entry_points_are_stable():
    assert list(inspect.signature(residual_qf).parameters) == ["state", "device", "bias_V"]
    assert list(inspect.signature(evaluate_qf).parameters) == ["state", "device", "bias_V"]
    assert [field.name for field in dataclasses.fields(ResidualEvaluation)] == [
        "f",
        "fluxes",
        "recombination_nd",
    ]


def test_csv_header_order_is_frozen():
    assert io.CSV_COLUMNS == (
        "run_id", "config_sha256", "solver", "bias_V", "record_type", "quantity",
        "location_kind", "index", "position_cm", "iteration", "unit", "numeric_value",
        "theory_value", "theory_kind", "theory_abs_error", "theory_rel_error",
        "reference_value", "reference_kind", "reference_abs_error", "reference_rel_error",
        "status", "notes",
    )
    assert io.header_line() == ",".join(io.CSV_COLUMNS)
    assert io.check_header(io.CSV_COLUMNS) == []


def test_run_result_contract():
    assert [field.name for field in dataclasses.fields(Attempt)] == [
        "bias_V", "status", "runtime_s", "result",
    ]
    assert [field.name for field in dataclasses.fields(RunResult)] == [
        "config", "device", "attempts", "successful", "decision", "failure_reason",
        "run_id", "config_sha256",
    ]
    empty = RunResult(config=None, device=None)
    assert empty.accepted is False and empty.decision == "REJECTED"


def test_equilibrium_initialisation_is_exported():
    assert callable(solve_equilibrium_scalar)


def test_run_entry_points_are_exported_at_top_level():
    """§10.2：GUI 与 CLI 只通过 `run`/`write_required_outputs` 发起正式计算。"""
    assert callable(minitcad.run) and callable(minitcad.write_required_outputs)
    signature = inspect.signature(minitcad.run)
    assert list(signature.parameters) == ["config", "progress", "should_stop"]
    assert signature.parameters["progress"].default is None
    # 取消回调的默认值必须可调用且返回 False（§10.2 的 `lambda: False` 语义）。
    assert callable(signature.parameters["should_stop"].default)
    assert signature.parameters["should_stop"].default() is False

    writer = inspect.signature(minitcad.write_required_outputs)
    assert list(writer.parameters) == ["run_result", "path", "overwrite"]
    assert writer.parameters["overwrite"].default is False
