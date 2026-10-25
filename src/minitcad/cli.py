"""命令行入口（W3 骨架）。

`validate` 已可用：校验 YAML 配置，或校验已有 comparison.csv 的表头。
`solve` 只报告“未接线”；在 W4-B4 接上 run() 之前不产生任何结果文件，
也不允许打印看似成功的 I-V。

用法（在仓库根目录）：

    PYTHONPATH=src python -m minitcad.cli validate examples/pn_left_p_right_n.yaml
    PYTHONPATH=src python -m minitcad.cli validate runs/pn_example/comparison.csv
    PYTHONPATH=src python -m minitcad.cli solve examples/pn_left_p_right_n.yaml
"""

from __future__ import annotations

import argparse
import pathlib
import sys

import yaml

from minitcad.core.device import build_device
from minitcad.io.config import load_config
from minitcad.io.results_csv import check_header, read_header

# 第 14.7 节要求 CLI 首屏写死方向与极性（ASCII，避免 Windows 控制台编码问题）。
BANNER = (
    "PN orientation: P(left, x=0) -> N(right, x=L) | "
    "V_app = V_left - V_right | +V = forward bias | right = 0 V"
)

EXIT_OK = 0
EXIT_REJECTED = 2
EXIT_NOT_WIRED = 3


def _declared_units(path: pathlib.Path) -> tuple[str, str]:
    """读回 YAML 里的输入/输出单位声明，仅用于首屏展示。"""
    try:
        raw = yaml.safe_load(path.read_text(encoding="utf-8"))
    except Exception:
        return "<unreadable>", "<unreadable>"
    if not isinstance(raw, dict):
        return "<unreadable>", "<unreadable>"
    unit_input = str((raw.get("units") or {}).get("system", "<undeclared>"))
    unit_output = str((raw.get("output") or {}).get("unit_system", "cm"))
    return unit_input, unit_output


def _validate_yaml(path: pathlib.Path) -> int:
    """校验 YAML 配置并打印启动摘要；失败时给出可读原因，不打印堆栈。"""
    try:
        config = load_config(path)
        device = build_device(config)
    except Exception as exc:  # 配置/材料/网格/区域错误都在这里统一拒绝
        print(f"REJECTED: {type(exc).__name__}: {exc}", file=sys.stderr)
        return EXIT_REJECTED

    unit_input, unit_output = _declared_units(path)
    print(BANNER)
    print(f"config       : {path}")
    print(f"input units  : {unit_input}")
    print("core units   : dimensionless")
    print(f"output units : {unit_output}")
    print(f"nodes        : {device.mesh.node_count}")
    print(f"T / V_T      : {device.temperature_k} K / {device.thermal_voltage_v:.8g} V")
    print(f"L_0 / D_0    : {device.l0_cm:.8g} cm / {device.d0_cm2_s:.8g} cm^2/s")
    print(
        f"C_L / C_R    : {device.net_doping_left_nd:.6g} / "
        f"{device.net_doping_right_nd:.6g} (normalized)"
    )
    print(f"material     : {device.material_parameter_policy} @ {device.temperature_k} K")
    return EXIT_OK


def _validate_csv(path: pathlib.Path) -> int:
    """只检查表头是否符合第 12.1 节冻结的 long-form schema，不重算数据。"""
    try:
        header = read_header(path)
    except OSError as exc:
        print(f"REJECTED: cannot read {path}: {exc}", file=sys.stderr)
        return EXIT_REJECTED
    problems = check_header(header)
    if problems:
        for problem in problems:
            print(f"REJECTED: {problem}", file=sys.stderr)
        return EXIT_REJECTED
    print(f"header OK: {len(header)} columns match the frozen long-form schema")
    return EXIT_OK


def _validate(args: argparse.Namespace) -> int:
    path = pathlib.Path(args.path)
    if not path.exists():
        print(f"REJECTED: no such file: {path}", file=sys.stderr)
        return EXIT_REJECTED
    if path.suffix.lower() == ".csv":
        return _validate_csv(path)
    return _validate_yaml(path)


def _solve(args: argparse.Namespace) -> int:
    """W3 只立骨架：明确未接线，绝不伪造结果。"""
    print(
        "solve: not wired yet (W4-B4 connects run()); no result files were written",
        file=sys.stderr,
    )
    print(f"solve: requested output directory ignored: {args.out}", file=sys.stderr)
    return EXIT_NOT_WIRED


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="minitcad",
        description="1D silicon PN drift-diffusion (dimensionless core, cm I/O by default)",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    validate = sub.add_parser("validate", help="validate a YAML config or a comparison.csv header")
    validate.add_argument("path", help="*.yaml/*.yml config, or an existing *.csv result")
    validate.set_defaults(func=_validate)

    solve = sub.add_parser("solve", help="solve a config (not wired until W4)")
    solve.add_argument("path", help="YAML config")
    solve.add_argument("--out", default="runs/pn_example", help="output directory")
    solve.set_defaults(func=_solve)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
