"""读取 YAML/字典配置并校验，返回 DeviceConfig。

用 units.system 声明输入是 cm 还是 si；si 在这里换算成 cm 再往下走。
字段后缀必须和声明一致，出现不认识的字段（例如没有单位的 density）就报错。
"""

from __future__ import annotations

from pathlib import Path

import yaml

from minitcad.core.device import (
    DEVICE_GENERIC_1D,
    DEVICE_PN,
    ORIENTATION_P_LEFT_N_RIGHT,
    DeviceConfig,
    MeshSpec,
    NormalizationConfig,
)
from minitcad.core.material import Material, get_material
from minitcad.core.region import Region
from minitcad.io.units import CM, SI, UnitError, convert, require_suffixed_field

# 每个物理量：(cm 字段名, si 字段名, 量类别)
_FIELDS: dict[str, tuple[str, str, str]] = {
    "length": ("length_cm", "length_m", "length"),
    "area": ("area_cm2", "area_m2", "area"),
    "intrinsic_density": ("intrinsic_density_cm3", "intrinsic_density_m3", "density"),
    "mobility_n": ("mobility_n_cm2_per_Vs", "mobility_n_m2_per_Vs", "mobility"),
    "mobility_p": ("mobility_p_cm2_per_Vs", "mobility_p_m2_per_Vs", "mobility"),
    "diffusivity": ("reference_diffusivity_cm2_s", "reference_diffusivity_m2_s", "diffusivity"),
    "mobility_ref": ("reference_mobility_cm2_per_Vs", "reference_mobility_m2_per_Vs", "mobility"),
    "region_start": ("start_cm", "start_m", "length"),
    "region_end": ("end_cm", "end_m", "length"),
    "acceptor": ("acceptor_cm3", "acceptor_m3", "density"),
    "donor": ("donor_cm3", "donor_m3", "density"),
}

# 每个块允许出现的键
_ALLOWED: dict[str, frozenset[str]] = {
    "root": frozenset(
        {"units", "device", "material", "regions", "mesh", "normalization", "contacts", "solver", "sweep", "output"}
    ),
    "device": frozenset(
        {"type", "orientation", "length_cm", "length_m", "area_cm2", "area_m2", "temperature_K"}
    ),
    "material": frozenset(
        {
            "name", "epsilon_relative", "temperature_K", "bandgap_eV",
            "lifetime_n_s", "lifetime_p_s",
            "intrinsic_density_cm3", "intrinsic_density_m3",
            "mobility_n_cm2_per_Vs", "mobility_n_m2_per_Vs",
            "mobility_p_cm2_per_Vs", "mobility_p_m2_per_Vs",
        }
    ),
    "region": frozenset(
        {
            "start_cm", "start_m", "end_cm", "end_m", "material",
            "acceptor_cm3", "acceptor_m3", "donor_cm3", "donor_m3",
        }
    ),
    "mesh": frozenset({"type", "nodes", "breakpoints", "per_segment_nodes"}),
    "normalization": frozenset(
        {
            "reference_diffusivity_cm2_s", "reference_diffusivity_m2_s",
            "reference_mobility_cm2_per_Vs", "reference_mobility_m2_per_Vs",
        }
    ),
    "contacts": frozenset(
        {"type", "left_terminal", "right_terminal", "right_voltage_V", "bias_definition"}
    ),
}

_VALID_SOLVER_METHODS = ("newton", "gummel")


def _check_keys(block: dict, block_name: str, where: str) -> None:
    """块里出现不认识的键就报错。"""
    unknown = sorted(set(block) - _ALLOWED[block_name])
    if unknown:
        raise UnitError(f"{where}: unknown field(s) {unknown}; declare units explicitly")


def _as_mapping(value: object, where: str) -> dict:
    if not isinstance(value, dict):
        raise UnitError(f"{where} must be a mapping, got {type(value).__name__}")
    return value


def _read(block: dict, logical: str, system: str) -> float | None:
    """读出一个物理量并换算成 cm；字段不存在返回 None。

    单位后缀必须和声明的单位制一致。
    """
    cm_key, si_key, kind = _FIELDS[logical]
    keys = [key for key in (cm_key, si_key) if key in block]
    if len(keys) > 1:
        raise UnitError(f"give only one of {cm_key!r}/{si_key!r}")
    if not keys:
        return None
    key = keys[0]
    require_suffixed_field(key, system)
    value = block[key]
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        raise UnitError(f"{key} must be a number, got {type(value).__name__}")
    return convert(kind, float(value), system, CM)


def _read_required(block: dict, logical: str, system: str) -> float:
    """和 _read 一样，但字段缺失时报错。"""
    value = _read(block, logical, system)
    if value is None:
        raise UnitError(f"missing required field {_FIELDS[logical][0]!r}")
    return value


def _float_or(value: object, fallback: float) -> float:
    """value 是数字就取它，否则用 fallback。"""
    if value is None:
        return fallback
    return float(value)


def _parse_material(block: dict, system: str) -> Material:
    """由 material 块构造 Material，没给的字段回落到参考材料。"""
    _check_keys(block, "material", "material")
    reference = get_material(str(block.get("name", "silicon")))
    ni = _read(block, "intrinsic_density", system)
    mu_n = _read(block, "mobility_n", system)
    mu_p = _read(block, "mobility_p", system)
    return Material(
        name=reference.name,
        symbol=reference.symbol,
        epsilon_relative=_float_or(block.get("epsilon_relative"), reference.epsilon_relative),
        reference_temperature_k=_float_or(
            block.get("temperature_K"), reference.reference_temperature_k
        ),
        intrinsic_density_cm3=_float_or(ni, reference.intrinsic_density_cm3),
        mobility_n_cm2_per_Vs=_float_or(mu_n, reference.mobility_n_cm2_per_Vs),
        mobility_p_cm2_per_Vs=_float_or(mu_p, reference.mobility_p_cm2_per_Vs),
        lifetime_n_s=_float_or(block.get("lifetime_n_s"), reference.lifetime_n_s),
        lifetime_p_s=_float_or(block.get("lifetime_p_s"), reference.lifetime_p_s),
        bandgap_eV=_float_or(block.get("bandgap_eV"), reference.bandgap_eV),
    )


def _parse_regions(blocks: object, system: str) -> tuple[Region, ...]:
    """由 regions 列表构造 Region 元组；重叠/缺口留给 DopingProfile 检查。"""
    if not isinstance(blocks, list) or not blocks:
        raise UnitError("regions must be a non-empty list")
    regions = []
    for index, raw in enumerate(blocks):
        where = f"regions[{index}]"
        block = _as_mapping(raw, where)
        _check_keys(block, "region", where)
        regions.append(
            Region(
                start_cm=_read_required(block, "region_start", system),
                end_cm=_read_required(block, "region_end", system),
                material=str(block.get("material", "silicon")),
                acceptor_cm3=_float_or(_read(block, "acceptor", system), 0.0),
                donor_cm3=_float_or(_read(block, "donor", system), 0.0),
            )
        )
    return tuple(regions)


def _parse_mesh(block: dict, system: str) -> MeshSpec:
    """由 mesh 块构造 MeshSpec。"""
    _check_keys(block, "mesh", "mesh")
    kind = str(block.get("type", "")).strip().lower()
    if kind == "uniform":
        return MeshSpec(kind="uniform", nodes=int(block["nodes"]))
    if kind == "segmented":
        breaks = tuple(_breakpoint(value, system) for value in block["breakpoints"])
        return MeshSpec(
            kind="segmented",
            breakpoints_cm=breaks,
            per_segment_nodes=tuple(int(n) for n in block["per_segment_nodes"]),
        )
    raise UnitError(f"mesh.type must be 'uniform' or 'segmented', got {kind!r}")


def _breakpoint(value: object, system: str) -> float:
    """分段点也是长度，按声明单位制换算成 cm。"""
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        raise UnitError("breakpoint values must be numbers")
    return convert("length", float(value), system, CM)


def _parse_normalization(block: dict, system: str) -> NormalizationConfig:
    """由 normalization 块构造 NormalizationConfig（D0 只能给一种）。"""
    _check_keys(block, "normalization", "normalization")
    return NormalizationConfig(
        reference_diffusivity_cm2_s=_read(block, "diffusivity", system),
        reference_mobility_cm2_per_Vs=_read(block, "mobility_ref", system),
    )


def _check_solver(block: object) -> None:
    """只检查 solver.method 的名字合法，W2 不执行求解。"""
    if not isinstance(block, dict):
        return
    method = block.get("method")
    if method is not None and str(method).strip().lower() not in _VALID_SOLVER_METHODS:
        raise UnitError(f"solver.method must be one of {_VALID_SOLVER_METHODS}, got {method!r}")


def parse_config_dict(data: dict) -> DeviceConfig:
    """由已解析的字典构造 DeviceConfig。"""
    root = _as_mapping(data, "config")
    _check_keys(root, "root", "config")

    units_block = root.get("units")
    if not isinstance(units_block, dict) or "system" not in units_block:
        raise UnitError("units.system must be declared explicitly ('cm' or 'si')")
    system = str(units_block["system"]).strip().lower()
    if system not in (CM, SI):
        raise UnitError(f"units.system must be 'cm' or 'si', got {system!r}")

    device_block = _as_mapping(root.get("device"), "device")
    _check_keys(device_block, "device", "device")
    material_block = _as_mapping(root.get("material"), "material")
    mesh_block = _as_mapping(root.get("mesh"), "mesh")
    normalization_block = _as_mapping(root.get("normalization"), "normalization")

    temperature_K = device_block.get("temperature_K")
    if not isinstance(temperature_K, (int, float)) or isinstance(temperature_K, bool):
        raise UnitError("device.temperature_K must be an explicit number")

    device_type = str(device_block.get("type", DEVICE_PN)).strip().lower()
    orientation = str(device_block.get("orientation", ORIENTATION_P_LEFT_N_RIGHT)).strip().lower()

    contacts = _as_mapping(root.get("contacts", {}), "contacts")
    if contacts:
        _check_keys(contacts, "contacts", "contacts")
    if device_type == DEVICE_PN and contacts:
        if str(contacts.get("left_terminal", "p")).lower() != "p" or str(
            contacts.get("right_terminal", "n")
        ).lower() != "n":
            raise UnitError("device.type pn requires left_terminal=p and right_terminal=n")
        if str(contacts.get("bias_definition", "left_minus_right")) != "left_minus_right":
            raise UnitError("contacts.bias_definition must be 'left_minus_right'")

    _check_solver(root.get("solver"))

    return DeviceConfig(
        length_cm=_read_required(device_block, "length", system),
        area_cm2=_read_required(device_block, "area", system),
        temperature_k=float(temperature_K),
        material=_parse_material(material_block, system),
        regions=_parse_regions(root.get("regions"), system),
        mesh=_parse_mesh(mesh_block, system),
        normalization=_parse_normalization(normalization_block, system),
        device_type=device_type,
        orientation=orientation,
        right_voltage_V=float(contacts.get("right_voltage_V", 0.0)),
    )


def load_config(path: str | Path) -> DeviceConfig:
    """从 YAML 文件读取配置。"""
    with open(path, "r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle)
    if not isinstance(data, dict):
        raise UnitError(f"{path}: YAML root must be a mapping")
    return parse_config_dict(data)
