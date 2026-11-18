"""控制台输出必须是纯 ASCII（2026-11-18 CI 失败的回归护栏）。

`run_mesh_convergence.py` 的表头里曾有一个 `Δ`：本机 GBK 能编码、GitHub 的 Windows
runner 用 cp1252 不能，脚本直接抛 `UnicodeEncodeError`。规律是——**程序的 stdout 是
ASCII 通道**，非 ASCII 只应该出现在写进文件的内容里（报告、CSV、PNG 文字）。

三条检查：
* 静态：`src/minitcad/`、`scripts/` 里 `print(...)`、`sys.stdout.write(...)`、
  argparse 的 `help/description/epilog` 不得含非 ASCII 字面量；
* 结构：argparse 的 `description` 必须是字面量（用 `__doc__` 就绕过了上面一条）；
* 动态：把 stdout 强制成 cp1252 跑 CLI 与各脚本的 `--help`，退出码必须为 0 且输出是 ASCII。
"""

from __future__ import annotations

import ast
import os
import pathlib
import subprocess
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
PRINT_NAMES = {"print"}
SINK_SUFFIXES = (".write",)
TEXT_KWARGS = ("help", "description", "epilog")
SCAN_GLOBS = ("src/minitcad/**/*.py", "scripts/*.py")
SCRIPTS = sorted(p.name for p in (ROOT / "scripts").glob("*.py") if not p.name.startswith("_"))


def _call_name(node: ast.expr) -> str:
    if isinstance(node, ast.Name):
        return node.id
    if isinstance(node, ast.Attribute):
        return "." + node.attr
    return ""


def _strings_in(node: ast.expr):
    """取出这个表达式里真正会变成文本的字符串片段（含 f-string 的常量部分）。"""
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        yield node.value
    elif isinstance(node, ast.JoinedStr):
        for value in node.values:
            if isinstance(value, ast.Constant) and isinstance(value.value, str):
                yield value.value


def _offending_lines(path: pathlib.Path) -> list[tuple[int, str]]:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    bad: list[tuple[int, str]] = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        name = _call_name(node.func)
        is_sink = name in PRINT_NAMES or name.endswith(SINK_SUFFIXES)
        is_help = any(kw.arg in TEXT_KWARGS and kw.value is not None for kw in node.keywords)
        if not (is_sink or is_help):
            continue
        for arg in list(node.args) + [kw.value for kw in node.keywords]:
            if arg is None:
                continue
            for text in _strings_in(arg):
                if not text.isascii():
                    bad.append((node.lineno, text[:60]))
    return bad


def _dynamic_descriptions(path: pathlib.Path) -> list[int]:
    """`description=__doc__` 会把中文 docstring 送进 stdout，静态检查看不见。"""
    tree = ast.parse(path.read_text(encoding="utf-8"))
    bad: list[int] = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        if not _call_name(node.func).endswith("ArgumentParser"):
            continue
        for kw in node.keywords:
            if kw.arg != "description":
                continue
            if not (isinstance(kw.value, ast.Constant) and isinstance(kw.value.value, str)):
                bad.append(node.lineno)
    return bad


def test_console_strings_are_ascii():
    offenders: list[str] = []
    scanned = 0
    for pattern in SCAN_GLOBS:
        for path in sorted(ROOT.glob(pattern)):
            scanned += 1
            for lineno, text in _offending_lines(path):
                offenders.append(f"{path.relative_to(ROOT)}:{lineno}: {text!r}")
    assert scanned > 0
    assert offenders == [], "console output must be ASCII:\n" + "\n".join(offenders)


def test_argparse_description_is_a_literal():
    offenders: list[str] = []
    for pattern in SCAN_GLOBS:
        for path in sorted(ROOT.glob(pattern)):
            for lineno in _dynamic_descriptions(path):
                offenders.append(f"{path.relative_to(ROOT)}:{lineno}")
    assert offenders == [], (
        "argparse description must be an ASCII string literal, not __doc__: "
        + ", ".join(offenders)
    )


def _run_cp1252(argv: list[str], tmp_path: pathlib.Path) -> subprocess.CompletedProcess:
    env = dict(os.environ)
    env["PYTHONIOENCODING"] = "cp1252"
    env["PYTHONPATH"] = os.pathsep.join([str(ROOT / "src"), env.get("PYTHONPATH", "")]).rstrip(
        os.pathsep
    )
    return subprocess.run(argv, cwd=ROOT, env=env, capture_output=True)


@pytest.mark.parametrize(
    "argv",
    [
        ["validate", "examples/pn_baseline_0_to_030.yaml"],
        ["solve", "examples/pn_baseline_0_to_030.yaml"],
    ],
)
def test_cli_output_survives_cp1252_stdout(argv, tmp_path):
    """cp1252 是 GitHub Windows runner 的默认控制台编码：非 ASCII 会直接抛异常。"""
    extra = ["--out", str(tmp_path / "run")] if argv[0] == "solve" else []
    completed = _run_cp1252([sys.executable, "-m", "minitcad.cli", *argv, *extra], tmp_path)
    assert completed.returncode == 0, completed.stderr.decode(errors="replace")
    assert completed.stdout.isascii(), completed.stdout.decode("cp1252", errors="replace")
    assert completed.stderr.isascii()


@pytest.mark.parametrize("script", SCRIPTS)
def test_script_help_survives_cp1252_stdout(script):
    completed = _run_cp1252(
        [sys.executable, str(ROOT / "scripts" / script), "--help"],
        pathlib.Path(ROOT),
    )
    assert completed.returncode == 0, completed.stderr.decode(errors="replace")
    assert completed.stdout.isascii(), completed.stdout.decode("cp1252", errors="replace")
