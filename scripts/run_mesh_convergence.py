"""生成 T12 的网格收敛表（§15 实验 C）。

```bash
python scripts/run_mesh_convergence.py                          # 默认基准配置
python scripts/run_mesh_convergence.py --config examples/pn_left_p_right_n.yaml
```

输出两张 markdown 表：加密序列（uniform 201/401/801）与同结点数的分布对比
（uniform vs junction-refined）。最细的 uniform 网格作为参考量。

`reports/week07_mesh_convergence.md` 里的数字就是本脚本的输出；表内每行都能重跑。
"""

from __future__ import annotations

import argparse
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from minitcad.core.device import MeshSpec  # noqa: E402
from minitcad.io.config import load_config  # noqa: E402
from minitcad.verification.reference import convergence_table, is_monotone, sample  # noqa: E402

DEFAULT_CONFIG = "examples/pn_baseline_0_to_030.yaml"
LADDER = (201, 401, 801)
REFINE_FACTORS = (1.0, 10.0, 20.0, 50.0)


def _row(cells: list[str]) -> str:
    return "| " + " | ".join(cells) + " |"


def _format_table(rows: list[dict]) -> list[str]:
    lines = [
        _row(["mesh", "M", "refine", "J [A/cm^2]", "|ΔJ|/J", "E_peak [V/cm]", "|ΔE|/E", "iters", "t [s]"]),
        _row(["---"] * 9),
    ]
    for row in rows:
        lines.append(
            _row(
                [
                    row["kind"],
                    str(row["nodes"]),
                    "-" if row["refine"] is None else f"{row['refine']:g}",
                    f"{row['current_density_A_cm2']:.9e}",
                    f"{row['current_rel_error']:.3e}",
                    f"{row['peak_field_V_cm']:.4e}",
                    f"{row['peak_field_rel_error']:.3e}",
                    f"{row['iterations']}",
                    f"{row['runtime_s']:.2f}",
                ]
            )
        )
    return lines


def _tolerate_unencodable_stdout() -> None:
    """表头含 `Δ`；cp1252 控制台（如 GitHub 的 Windows runner）没有这个字符。

    只在写不出去时退化成 `?`，不改变能正常编码的控制台（GBK/UTF-8）的输出。
    """
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(errors="replace")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=DEFAULT_CONFIG)
    args = parser.parse_args(argv)
    _tolerate_unencodable_stdout()

    config = load_config(pathlib.Path(args.config))
    reference = sample(config, MeshSpec(kind="uniform", nodes=LADDER[-1]), label="reference")
    if not reference.converged:
        print("reference mesh did not converge", file=sys.stderr)
        return 1

    print(f"# T12 mesh convergence — {args.config}")
    print(f"\nreference: uniform {LADDER[-1]} nodes, J = {reference.current_density_A_cm2:.9e} A/cm^2, "
          f"bias = {reference.bias_V:g} V, solver = {config.solver.method}")

    ladder = [
        sample(config, MeshSpec(kind="uniform", nodes=nodes), label=f"uniform {nodes}")
        for nodes in LADDER
    ]
    print("\n## A. Refinement ladder (uniform)\n")
    ladder_rows = convergence_table(ladder, reference)
    print("\n".join(_format_table(ladder_rows)))
    errors = [row["current_rel_error"] for row in ladder_rows]
    print(f"\nmonotone decreasing: {is_monotone(errors)}; "
          f"two finest relative difference: {errors[-2]:.3e}")

    print(f"\n## B. Distribution at fixed node count ({LADDER[0]} nodes)\n")
    distribution = [
        sample(config, MeshSpec(kind="uniform", nodes=LADDER[0]), label="uniform")
    ] + [
        sample(
            config,
            MeshSpec(kind="refined", nodes=LADDER[0], refine=refine),
            label=f"refined refine={refine:g}",
        )
        for refine in REFINE_FACTORS
    ]
    print("\n".join(_format_table(convergence_table(distribution, reference))))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
