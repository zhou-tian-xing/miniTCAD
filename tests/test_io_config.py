"""W2：YAML 配置校验，以及「同一物理器件 cm 输入与 SI 输入得到同一无量纲状态」（第 3.3、13 节）。"""

from __future__ import annotations

import pathlib

import pytest

from minitcad.core.device import build_device
from minitcad.io.config import load_config, parse_config_dict
from minitcad.io.units import UnitError

EXAMPLES = pathlib.Path(__file__).resolve().parents[1] / "examples"
CM_YAML = EXAMPLES / "pn_left_p_right_n.yaml"
SI_YAML = EXAMPLES / "pn_left_p_right_n_si.yaml"


def test_example_yaml_loads():
    config = load_config(CM_YAML)
    assert config.length_cm == pytest.approx(1.0e-3)
    assert config.device_type == "pn"
    assert config.regions[0].acceptor_cm3 == pytest.approx(1.0e16)


def test_cm_and_si_inputs_give_the_same_normalized_device():
    cm_device = build_device(load_config(CM_YAML))
    si_device = build_device(load_config(SI_YAML))
    assert si_device.l0_cm == pytest.approx(cm_device.l0_cm, rel=1e-12)
    assert si_device.thermal_voltage_v == pytest.approx(cm_device.thermal_voltage_v, rel=1e-15)
    assert si_device.area_cm2 == pytest.approx(cm_device.area_cm2, rel=1e-12)
    assert si_device.mesh.node_count == cm_device.mesh.node_count
    assert si_device.mesh.x_nodes_cm == pytest.approx(cm_device.mesh.x_nodes_cm, rel=1e-12)
    assert si_device.net_doping_nd == pytest.approx(cm_device.net_doping_nd, rel=1e-12)
    assert si_device.cbar_nd == pytest.approx(cm_device.cbar_nd, rel=1e-12)
    assert si_device.contact_targets(0.3) == pytest.approx(cm_device.contact_targets(0.3))


def _minimal(**device_overrides) -> dict:
    return {
        "units": {"system": "cm"},
        "device": {
            "type": "pn",
            "orientation": "p_left_n_right",
            "length_cm": 1.0e-3,
            "area_cm2": 1.0e-4,
            "temperature_K": 300.0,
            **device_overrides,
        },
        "material": {"name": "silicon"},
        "regions": [
            {"start_cm": 0.0, "end_cm": 2.5e-4, "material": "silicon", "acceptor_cm3": 1.0e16, "donor_cm3": 0.0},
            {"start_cm": 2.5e-4, "end_cm": 1.0e-3, "material": "silicon", "acceptor_cm3": 0.0, "donor_cm3": 2.0e15},
        ],
        "mesh": {"type": "uniform", "nodes": 21},
        "normalization": {"reference_diffusivity_cm2_s": 10.0},
    }


def test_missing_units_system_is_rejected():
    data = _minimal()
    del data["units"]
    with pytest.raises(UnitError):
        parse_config_dict(data)


def test_unsuffixed_physical_field_is_rejected():
    data = _minimal()
    data["device"]["length"] = 1.0e-3  # 无单位后缀
    del data["device"]["length_cm"]
    with pytest.raises(UnitError):
        parse_config_dict(data)


def test_mixed_unit_systems_in_one_config_are_rejected():
    data = _minimal()
    data["device"]["length_m"] = 1.0e-5  # 声明 cm 却给了 SI 后缀
    del data["device"]["length_cm"]
    with pytest.raises(UnitError):
        parse_config_dict(data)


def test_negative_temperature_is_rejected():
    with pytest.raises(ValueError):
        parse_config_dict(_minimal(temperature_K=-1.0))


def test_reversed_pn_regions_are_rejected():
    data = _minimal()
    data["regions"] = [
        {"start_cm": 0.0, "end_cm": 2.5e-4, "material": "silicon", "donor_cm3": 1.0e16},
        {"start_cm": 2.5e-4, "end_cm": 1.0e-3, "material": "silicon", "acceptor_cm3": 2.0e15},
    ]
    with pytest.raises(ValueError):
        parse_config_dict(data)


def test_both_normalization_references_are_rejected():
    data = _minimal()
    data["normalization"] = {
        "reference_diffusivity_cm2_s": 10.0,
        "reference_mobility_cm2_per_Vs": 1000.0,
    }
    with pytest.raises(ValueError):
        parse_config_dict(data)


def test_unknown_solver_method_is_rejected():
    data = _minimal()
    data["solver"] = {"method": "relaxation"}
    with pytest.raises(UnitError):
        parse_config_dict(data)


def test_unknown_material_name_is_rejected():
    data = _minimal()
    data["material"] = {"name": "germanium"}
    with pytest.raises(ValueError):
        parse_config_dict(data)
