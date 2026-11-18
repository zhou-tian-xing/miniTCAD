"""一条命令跑完 Gate 2 / Gate 1 的验收项，打印 PASS/FAIL 表（§11.2、W6-B2）。

```bash
python scripts/run_gates.py              # 全部
python scripts/run_gates.py --gate 2     # 只看 Gate 2
python scripts/run_gates.py --list       # 只列清单
```

每一项对应一个测试文件（必要时加 `-k`），通过就直接调 `pytest` 子进程跑，
失败时把该子进程的输出尾巴贴出来——不做二次解析，避免"脚本说通过了、人却不知道跑了什么"。

退出码：0 = 所有被选中的项通过；1 = 有失败。
"""

from __future__ import annotations

import argparse
import pathlib
import subprocess
import sys
from dataclasses import dataclass

ROOT = pathlib.Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class Item:
    """一个验收项：ID、所属 Gate、说明、pytest 参数。"""

    test_id: str
    gate: int
    title: str
    args: tuple[str, ...]


ITEMS: tuple[Item, ...] = (
    Item("T08", 1, "0 V equilibrium: mass action / zero current / contact neutrality",
         ("tests/test_week04_gate1.py",)),
    Item("T09", 1, "total current conserved along x", ("tests/test_week04_gate1.py",)),
    Item("T10", 2, "Jacobian vs independent finite differences", ("tests/test_t10_jacobian.py",)),
    Item("T11", 2, "Gummel vs Newton cross-check", ("tests/test_t11_dual_solver.py",)),
    Item("T13", 2, "adaptive rollback and failure records", ("tests/test_t13_rollback.py",)),
    Item("T12", 3, "mesh convergence ladder", ("tests/test_t12_mesh.py",)),
    Item("T14", 3, "three-file schema and independent acceptor", ("tests/test_t14_output.py",)),
    Item("T15", 3, "reproducibility", ("tests/test_t15_reproducible.py",)),
    Item("API", 3, "public API / YAML / CLI freeze",
         ("tests/test_w07_freeze.py", "tests/test_public_api.py")),
)


def _run(item: Item, pytest_args: list[str]) -> tuple[bool, str]:
    command = [sys.executable, "-m", "pytest", *item.args, "-q", *pytest_args]
    completed = subprocess.run(
        command, cwd=ROOT, capture_output=True, text=True, encoding="utf-8", errors="replace"
    )
    output = (completed.stdout or "") + (completed.stderr or "")
    return completed.returncode == 0, output.strip().splitlines()[-1] if output.strip() else ""


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run the acceptance items and print a table.")
    parser.add_argument("--gate", type=int, choices=(1, 2, 3), default=None)
    parser.add_argument("--list", action="store_true", help="only print the item list")
    parser.add_argument("-k", dest="keyword", default=None, help="extra -k filter for pytest")
    args = parser.parse_args(argv)

    selected = [item for item in ITEMS if args.gate is None or item.gate == args.gate]
    if args.list:
        for item in selected:
            print(f"  {item.test_id:4s} Gate {item.gate}  {item.title}")
        return 0

    extra = ["-k", args.keyword] if args.keyword else []
    print(f"running {len(selected)} acceptance item(s) from {ROOT}")
    print(f"{'ID':5s} {'gate':4s} {'result':7s} summary")
    failures: list[str] = []
    for item in selected:
        passed, summary = _run(item, extra)
        print(f"{item.test_id:5s} {item.gate:<4d} {'PASS' if passed else 'FAIL':7s} {summary}")
        if not passed:
            failures.append(item.test_id)

    print()
    if failures:
        print(f"REJECTED: {len(failures)} item(s) failed -> {failures}")
        print("rerun one item with: python -m pytest " + (selected[0].args[0] if selected else "tests/") + " -q")
        return 1
    print(f"ACCEPTED: {len(selected)} item(s) passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
