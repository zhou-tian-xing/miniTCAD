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
    SOLVER_METHODS,
    OutputConfig,
    SolverConfig,
    SweepConfig,
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
        return MeshSpec(kind="uniform", nodes=_node_count(block["nodes"]))
    if kind == "segmented":
        breaks = tuple(_breakpoint(value, system) for value in block["breakpoints"])
        return MeshSpec(
            kind="segmented",
            breakpoints_cm=breaks,
            per_segment_nodes=tuple(_node_count(n) for n in block["per_segment_nodes"]),
        )
    raise UnitError(f"mesh.type must be 'uniform' or 'segmented', got {kind!r}")


def _node_count(value: object) -> int:
    """结点数必须是整数；浮点 201.7 不能被静默截断成 201（第 5.4 节）。"""
    if isinstance(value, bool) or not isinstance(value, int):
        raise UnitError(f"node counts must be integers, got {value!r}")
    return value


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


def _parse_solver(block: object) -> SolverConfig:
    """`solver` 块；缺省时用 dataclass 的默认值。"""
    if block is None:
        return SolverConfig()
    if not isinstance(block, dict):
        raise UnitError("solver must be a mapping")
    allowed = {"method", "jacobian", "tol_residual", "tol_update", "max_iter",
               "alpha_min", "positive_scale_floor", "n_floor"}
    unknown = sorted(set(block) - allowed)
    if unknown:
        raise UnitError(f"solver: unknown field(s) {unknown}")
    method = str(block.get("method", "gummel")).strip().lower()
    if method not in SOLVER_METHODS:
        raise UnitError(f"solver.method must be one of {SOLVER_METHODS}, got {method!r}")
    return SolverConfig(
        method=method,
        jacobian=str(block.get("jacobian", "sparse_fd")).strip().lower(),
        tol_residual=_positive_float(block.get("tol_residual", 1.0e-8), "solver.tol_residual"),
        tol_update=_positive_float(block.get("tol_update", 1.0e-6), "solver.tol_update"),
        max_iter=_integer(block.get("max_iter", 60), "solver.max_iter"),
        alpha_min=_positive_float(block.get("alpha_min", 1.0e-6), "solver.alpha_min"),
        n_floor=_positive_float(
            block.get("n_floor", block.get("positive_scale_floor", 1.0e-15)), "solver.n_floor"
        ),
    )


def _parse_sweep(block: object) -> SweepConfig:
    """`sweep` 块。`step_V` 缺省时取 `initial_step_V`（两种写法等价，W2 起就接受）。

    `adaptive` 打开 W5 的自适应继续法；`min_step_V/max_step_V/grow_factor/shrink_factor`
    只在自适应模式下起作用，但**始终解析并校验**——留着不校验的话，一个写错的自适应参数
    会在固定步进模式下静默通过，等切换模式时才炸。
    """
    if block is None:
        return SweepConfig()
    if not isinstance(block, dict):
        raise UnitError("sweep must be a mapping")
    allowed = {"start_V", "stop_V", "step_V", "initial_step_V", "min_step_V",
               "max_step_V", "grow_factor", "shrink_factor", "adaptive"}
    unknown = sorted(set(block) - allowed)
    if unknown:
        raise UnitError(f"sweep: unknown field(s) {unknown}")
    step = block.get("step_V", block.get("initial_step_V", 0.05))
    defaults = SweepConfig()
    adaptive = block.get("adaptive", False)
    if not isinstance(adaptive, bool):
        raise UnitError("sweep.adaptive must be a boolean")
    try:
        return SweepConfig(
            start_V=_float(block.get("start_V", 0.0), "sweep.start_V"),
            stop_V=_float(block.get("stop_V", 0.3), "sweep.stop_V"),
            step_V=_float(step, "sweep.step_V"),
            adaptive=adaptive,
            min_step_V=_float(block.get("min_step_V", defaults.min_step_V), "sweep.min_step_V"),
            max_step_V=_float(block.get("max_step_V", defaults.max_step_V), "sweep.max_step_V"),
            grow_factor=_float(block.get("grow_factor", defaults.grow_factor), "sweep.grow_factor"),
            shrink_factor=_float(
                block.get("shrink_factor", defaults.shrink_factor), "sweep.shrink_factor"
            ),
        )
    except ValueError as exc:
        raise UnitError(f"sweep: {exc}") from exc


def _parse_output(block: object) -> OutputConfig:
    """`output` 块：目录与文件单位制；单位制只影响 writer（§12.1）。"""
    if block is None:
        return OutputConfig()
    if not isinstance(block, dict):
        raise UnitError("output must be a mapping")
    allowed = {"directory", "unit_system"}
    unknown = sorted(set(block) - allowed)
    if unknown:
        raise UnitError(f"output: unknown field(s) {unknown}")
    return OutputConfig(
        directory=str(block.get("directory", "runs/pn_example")),
        unit_system=str(block.get("unit_system", "cm")).strip().lower(),
    )


def _float(value: object, where: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise UnitError(f"{where} must be a number, got {value!r}")
    return float(value)


def _positive_float(value: object, where: str) -> float:
    number = _float(value, where)
    if number <= 0.0:
        raise UnitError(f"{where} must be > 0, got {number!r}")
    return number


def _integer(value: object, where: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise UnitError(f"{where} must be an integer, got {value!r}")
    return value


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
        right_voltage_V = contacts.get("right_voltage_V", 0.0)
        if isinstance(right_voltage_V, bool) or not isinstance(right_voltage_V, (int, float)):
            raise UnitError("contacts.right_voltage_V must be a number")
        if float(right_voltage_V) != 0.0:
            raise UnitError(
                "contacts.right_voltage_V must be 0 V: the right terminal is the 0 V reference "
                "and the bias is applied at the left terminal (V_app = V_left - V_right)"
            )

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
        input_unit_system=system,
        solver=_parse_solver(root.get("solver")),
        sweep=_parse_sweep(root.get("sweep")),
        output=_parse_output(root.get("output")),
    )


def load_config(path: str | Path) -> DeviceConfig:
    """从 YAML 文件读取配置。"""
    with open(path, "r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle)
    if not isinstance(data, dict):
        raise UnitError(f"{path}: YAML root must be a mapping")
    return parse_config_dict(data)


def config_payload(config: DeviceConfig) -> dict:
    """规范化配置字典：字段全、顺序无关（`sort_keys`），供指纹与复现审计。"""
    return {
        "length_cm": config.length_cm,
        "area_cm2": config.area_cm2,
        "temperature_K": config.temperature_k,
        "input_unit_system": config.input_unit_system,
        "device": {
            "type": config.device_type,
            "orientation": config.orientation,
            "right_voltage_V": config.right_voltage_V,
        },
        "material": {
            "name": config.material.name,
            "epsilon_relative": config.material.epsilon_relative,
            "reference_temperature_K": config.material.reference_temperature_k,
            "intrinsic_density_cm3": config.material.intrinsic_density_cm3,
            "mobility_n_cm2_per_Vs": config.material.mobility_n_cm2_per_Vs,
            "mobility_p_cm2_per_Vs": config.material.mobility_p_cm2_per_Vs,
            "lifetime_n_s": config.material.lifetime_n_s,
            "lifetime_p_s": config.material.lifetime_p_s,
            "bandgap_eV": config.material.bandgap_eV,
        },
        "regions": [
            {
                "start_cm": region.start_cm,
                "end_cm": region.end_cm,
                "material": region.material,
                "acceptor_cm3": region.acceptor_cm3,
                "donor_cm3": region.donor_cm3,
            }
            for region in config.regions
        ],
        "mesh": {
            "kind": config.mesh.kind,
            "nodes": config.mesh.nodes,
            "breakpoints_cm": list(config.mesh.breakpoints_cm),
            "per_segment_nodes": list(config.mesh.per_segment_nodes),
        },
        "normalization": {
            "reference_diffusivity_cm2_s": config.normalization.reference_diffusivity_cm2_s,
            "reference_mobility_cm2_per_Vs": config.normalization.reference_mobility_cm2_per_Vs,
        },
        "solver": {
            "method": config.solver.method,
            "jacobian": config.solver.jacobian,
            "tol_residual": config.solver.tol_residual,
            "tol_update": config.solver.tol_update,
            "max_iter": config.solver.max_iter,
            "alpha_min": config.solver.alpha_min,
            "n_floor": config.solver.n_floor,
        },
        "sweep": {
            "start_V": config.sweep.start_V,
            "stop_V": config.sweep.stop_V,
            "step_V": config.sweep.step_V,
        },
        "output": {
            "directory": config.output.directory,
            "unit_system": config.output.unit_system,
        },
    }


def config_fingerprint(config: DeviceConfig) -> str:
    """规范化配置的 SHA-256，写进 CSV 的 `config_sha256` 与报告的 Run summary。"""
    import hashlib
    import json

    payload = json.dumps(config_payload(config), sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()
