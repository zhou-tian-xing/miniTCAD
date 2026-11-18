"""w07-core-freeze：公开 API、YAML schema 与字段表的冻结契约。

`w07-core-freeze` 之后只允许 GUI、联调与阻塞性 bug 修复，不允许破坏性改名。本文件把
「不能改名/不能改顺序」的部分钉死；改这里就必须同时改 `docs/api_reference.md` 与
`docs/gui_api_contract.md`，并走接口会审的变更记录。
"""

from __future__ import annotations

import dataclasses
import inspect
import pathlib

import pytest

import minitcad
from minitcad.api import run, write_required_outputs
from minitcad.core.device import MeshSpec, SweepConfig
from minitcad.io.config import (
    config_fingerprint,
    load_config,
    parse_config_dict,
    save_config,
)
from minitcad.simulation.run_result import Attempt, RunResult

EXAMPLES = pathlib.Path(__file__).resolve().parents[1] / "examples"
BASELINE = EXAMPLES / "pn_baseline_0_to_030.yaml"
SI = EXAMPLES / "pn_left_p_right_n_si.yaml"
NEWTON = EXAMPLES / "pn_newton_adaptive_0_to_030.yaml"


def test_package_exports_are_frozen():
    assert minitcad.__all__ == [
        "run",
        "write_required_outputs",
        "load_config",
        "parse_config_dict",
        "save_config",
        "config_fingerprint",
    ]


@pytest.mark.parametrize(
    "function, expected",
    [
        (run, ["config", "progress", "should_stop"]),
        (write_required_outputs, ["run_result", "path", "overwrite"]),
        (load_config, ["path"]),
        (parse_config_dict, ["data"]),
        (save_config, ["config", "path"]),
        (config_fingerprint, ["config"]),
    ],
)
def test_public_signatures_are_frozen(function, expected):
    assert list(inspect.signature(function).parameters) == expected


def test_result_fields_are_frozen():
    assert [field.name for field in dataclasses.fields(RunResult)] == [
        "config", "device", "attempts", "successful", "decision", "failure_reason",
        "run_id", "config_sha256",
    ]
    assert [field.name for field in dataclasses.fields(Attempt)] == [
        "bias_V", "status", "runtime_s", "result", "note",
    ]
    assert [field.name for field in dataclasses.fields(SweepConfig)] == [
        "start_V", "stop_V", "step_V", "adaptive", "min_step_V", "max_step_V",
        "grow_factor", "shrink_factor",
    ]
    assert [field.name for field in dataclasses.fields(MeshSpec)] == [
        "kind", "nodes", "breakpoints_cm", "per_segment_nodes", "refine",
    ]


@pytest.mark.parametrize("example", [BASELINE, NEWTON])
def test_save_config_round_trips_to_the_same_fingerprint(example, tmp_path):
    """GUI 的“保存”必须能被 `load_config` 读回同一份规范化配置（T16 的前置）。"""
    config = load_config(example)
    target = save_config(config, tmp_path / "saved.yaml")
    again = load_config(target)
    assert config_fingerprint(again) == config_fingerprint(config)


def test_save_config_normalises_si_input_to_cm(tmp_path):
    """SI 声明的配置保存后按 cm 写：物理不变，`input_unit_system` 变成 cm。

    指纹随之改变是有意的——文件里声明的单位本来就该进入指纹。这里断言的是
    “归一化后是同一个器件”，而不是“同一个指纹”。
    """
    config = load_config(SI)
    assert config.input_unit_system == "si"
    again = load_config(save_config(config, tmp_path / "saved.yaml"))
    assert again.input_unit_system == "cm"
    assert again.length_cm == pytest.approx(config.length_cm, rel=0, abs=0)
    assert again.material.intrinsic_density_cm3 == pytest.approx(
        config.material.intrinsic_density_cm3
    )
    for left, right in zip(config.regions, again.regions):
        assert right.acceptor_cm3 == pytest.approx(left.acceptor_cm3)
        assert right.donor_cm3 == pytest.approx(left.donor_cm3)


def test_save_config_keeps_a_refined_mesh(tmp_path):
    config = dataclasses.replace(
        load_config(BASELINE), mesh=MeshSpec(kind="refined", nodes=201, refine=20.0)
    )
    again = load_config(save_config(config, tmp_path / "refined.yaml"))
    assert again.mesh.kind == "refined"
    assert again.mesh.refine == pytest.approx(20.0)
    assert config_fingerprint(again) == config_fingerprint(config)


def test_fingerprint_covers_the_mesh_and_sweep_options():
    """指纹必须覆盖全部求解相关的配置，否则两次不同的运行会共用同一个 SHA。"""
    base = load_config(BASELINE)
    variants = [
        dataclasses.replace(base, mesh=dataclasses.replace(base.mesh, refine=10.0)),
        dataclasses.replace(
            base, sweep=dataclasses.replace(base.sweep, adaptive=True)
        ),
        dataclasses.replace(
            base, sweep=dataclasses.replace(base.sweep, min_step_V=1e-4)
        ),
        dataclasses.replace(
            base, solver=dataclasses.replace(base.solver, method="newton")
        ),
    ]
    fingerprints = {config_fingerprint(variant) for variant in variants}
    assert len(fingerprints) == len(variants)
    assert config_fingerprint(base) not in fingerprints


def test_yaml_schema_rejects_unknown_fields():
    """未知字段必须报错：GUI 写错名字时应当当场失败，而不是被静默忽略。"""
    from minitcad.io.units import UnitError

    data = {
        "units": {"system": "cm"},
        "device": {
            "type": "pn",
            "orientation": "p_left_n_right",
            "length_cm": 1.0e-3,
            "area_cm2": 1.0e-4,
            "temperature_K": 300.0,
        },
        "material": {
            "name": "silicon",
            "epsilon_relative": 11.7,
            "intrinsic_density_cm3": 1.0e10,
            "mobility_n_cm2_per_Vs": 1350.0,
            "mobility_p_cm2_per_Vs": 480.0,
            "lifetime_n_s": 1.0e-6,
            "lifetime_p_s": 2.0e-6,
            "bandgap_eV": 1.12,
        },
        "regions": [
            {"start_cm": 0.0, "end_cm": 1.0e-4, "acceptor_cm3": 1.0e16, "donor_cm3": 0.0},
            {"start_cm": 1.0e-4, "end_cm": 1.0e-3, "acceptor_cm3": 0.0, "donor_cm3": 2.0e15},
        ],
        "mesh": {"type": "uniform", "nodes": 51},
        "normalization": {"reference_diffusivity_cm2_s": 10.0},
        "contacts": {"right_voltage_V": 0.0},
        "sweep": {"stop_V": 0.1, "adaptive": True},
    }
    assert parse_config_dict(data).sweep.adaptive is True
    data["sweep"] = {"stop_V": 0.1, "adaptiv": True}
    with pytest.raises(UnitError):
        parse_config_dict(data)


def test_cli_exit_codes_are_frozen():
    """退出码是 GUI/CI 的契约：0 接受、2 拒绝、3 未实现。"""
    from minitcad import cli

    assert (cli.EXIT_OK, cli.EXIT_REJECTED, cli.EXIT_NOT_WIRED) == (0, 2, 3)
    assert "P(left, x=0) -> N(right, x=L)" in cli.BANNER
