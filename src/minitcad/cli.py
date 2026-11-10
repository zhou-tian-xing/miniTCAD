"""命令行入口：`minitcad validate` 与 `minitcad solve`（§10.2、§14.7）。

`solve` 与 GUI 走同一个 `api.run()` + `api.write_required_outputs()`，
因此两端的偏压点、失败原因与三文件内容必然一致。

退出码：0 = ACCEPTED，2 = REJECTED/配置非法，3 = 求解器未落地。
W5 起 `method: newton` 已接线，3 只作为**保留码**：`run()` 若对某个尚未实现的能力
抛 `NotImplementedError`，CLI 必须报 `NOT_WIRED` 而不是伪装成求解失败（2）。
"""

from __future__ import annotations

import argparse
import pathlib
import sys

import yaml

from minitcad.api import run, write_required_outputs
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
    print(f"solver       : {config.solver.method} (max_iter {config.solver.max_iter})")
    print(
        f"sweep        : {config.sweep.start_V} -> {config.sweep.stop_V} V, "
        f"step {config.sweep.step_V} V"
    )
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
    """求解并写出三个固定文件；失败也照实打印原因，不产生假的成功摘要。"""
    path = pathlib.Path(args.path)
    try:
        config = load_config(path)
    except Exception as exc:
        print(f"REJECTED: {type(exc).__name__}: {exc}", file=sys.stderr)
        return EXIT_REJECTED

    directory = args.out or config.output.directory
    print(BANNER)
    print(f"config  : {path}")
    print(f"output  : {directory}")
    print(f"solver  : {config.solver.method}")

    def progress(bias_V: float, status: str) -> None:
        print(f"  [{status:4s}] V_app = {bias_V:.4f} V")

    try:
        result = run(config, progress=progress)
    except NotImplementedError as exc:
        print(f"NOT_WIRED: {exc}", file=sys.stderr)
        return EXIT_NOT_WIRED
    except Exception as exc:
        print(f"REJECTED: {type(exc).__name__}: {exc}", file=sys.stderr)
        return EXIT_REJECTED

    try:
        written = write_required_outputs(result, directory, overwrite=args.overwrite)
    except FileExistsError as exc:
        print(f"REJECTED: {exc}", file=sys.stderr)
        return EXIT_REJECTED

    print(f"decision: {result.decision}")
    if result.failure_reason is not None:
        print(f"reason  : {result.failure_reason}")
    print(f"run_id  : {result.run_id}")
    # 自适应模式下 `attempts` 可能包含同一区段的多次尝试（回滚），因此这里说的是
    # “尝试”而不是“偏压点”：把回滚次数如实暴露出来，别让 9/12 看起来像丢了 3 个偏压。
    print(f"accepted: {len(result.successful)} of {len(result.attempts)} attempts")
    if len(result.attempts) != len(result.successful):
        print(f"rollback: {len(result.attempts) - len(result.successful)} failed attempt(s) kept in comparison.csv")
    for name, written_path in zip(("comparison.csv", "comparison.png", "convergence.md"), written):
        print(f"  wrote {name:16s} {written_path}")
    return EXIT_OK if result.accepted else EXIT_REJECTED


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="minitcad",
        description="1D silicon PN drift-diffusion (dimensionless core, cm I/O by default)",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    validate = sub.add_parser("validate", help="validate a YAML config or a comparison.csv header")
    validate.add_argument("path", help="*.yaml/*.yml config, or an existing *.csv result")
    validate.set_defaults(func=_validate)

    solve = sub.add_parser("solve", help="solve a config and write the three required files")
    solve.add_argument("path", help="YAML config")
    solve.add_argument("--out", default=None, help="output directory (overrides output.directory)")
    solve.add_argument(
        "--overwrite",
        action="store_true",
        help="replace existing comparison.csv/comparison.png/convergence.md",
    )
    solve.set_defaults(func=_solve)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
