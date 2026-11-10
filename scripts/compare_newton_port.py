"""核对根包并入的 Newton 与 `C2_v1` 里 A 的实现是否逐位一致。

两个包同名（`minitcad`），不能在同一进程里同时导入，所以分两步跑、再用 JSON 对比：

```bash
python scripts/compare_newton_port.py run --src src       --out .tmp/root.json
python scripts/compare_newton_port.py run --src C2_v1/src --out .tmp/c2v1.json
python scripts/compare_newton_port.py diff .tmp/root.json .tmp/c2v1.json
```

`diff` 比较每个偏压的 `converged / iterations / residual_norm / update_norm / reason`，
以及最终 `u/n/p` 三个数组的 SHA-256。任何一项不同就以退出码 1 结束。

为什么值得单独留一个脚本：并入 A 的求解器时，B 改了文档、异常类别和参数表，
“数值没动”是必须能被复核的断言，而不是一句说明。W6-B2 的
`scripts/compare_solvers.py` 会复用这里的 JSON 形状。
"""

from __future__ import annotations

import argparse
import hashlib
import json
import pathlib
import sys

# 序列解析、状态构造与 JSON 形状必须与两个包都兼容，因此只用两边共有的 API：
# core.device 的 MeshSpec/NormalizationConfig/DeviceConfig/build_device、
# simulation.equilibrium 的标量初值、solvers.newton 的 solve_newton。
_BIASES = (0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3)


def _build_device():
    from minitcad.core.device import (
        DeviceConfig,
        MeshSpec,
        NormalizationConfig,
        build_device,
    )
    from minitcad.core.material import Material
    from minitcad.core.region import Region

    # 与 examples/pn_baseline_0_to_030.yaml 同一器件，但直接构造，
    # 避免依赖各包 YAML 解析器是否支持 solver/sweep 块。
    material = Material(
        name="silicon",
        symbol="Si",
        reference_temperature_k=300.0,
        epsilon_relative=11.7,
        intrinsic_density_cm3=1.0e10,
        mobility_n_cm2_per_Vs=1350.0,
        mobility_p_cm2_per_Vs=480.0,
        lifetime_n_s=1.0e-6,
        lifetime_p_s=2.0e-6,
        bandgap_eV=1.12,
    )
    config = DeviceConfig(
        length_cm=1.0e-3,
        area_cm2=1.0e-4,
        temperature_k=300.0,
        material=material,
        regions=(
            Region(start_cm=0.0, end_cm=2.5e-4, material=material, acceptor_cm3=1.0e16, donor_cm3=0.0),
            Region(start_cm=2.5e-4, end_cm=1.0e-3, material=material, acceptor_cm3=0.0, donor_cm3=2.0e15),
        ),
        mesh=MeshSpec(kind="uniform", nodes=201),
        normalization=NormalizationConfig(reference_diffusivity_cm2_s=10.0),
    )
    return build_device(config)


def _sha(array) -> str:
    import numpy as np

    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()


def _run(src: str, out_path: pathlib.Path) -> None:
    sys.path.insert(0, src)
    from minitcad.simulation.equilibrium import solve_equilibrium_scalar
    from minitcad.solvers.newton import solve_newton

    device = _build_device()
    targets = device.contact_targets(0.0)
    state = solve_equilibrium_scalar(
        device.mesh.x_nodes_nd, device.cbar_nd, targets[0], targets[3]
    ).state

    records = []
    for bias in _BIASES:
        result = solve_newton(device, bias, state)
        record = {
            "bias": bias,
            "converged": bool(result.converged),
            "iterations": int(result.iterations),
            "residual_norm": float(result.residual_norm),
            "update_norm": float(result.update_norm),
            "reason": result.reason,
        }
        if result.state is not None:
            record["sha_u"] = _sha(result.state.u)
            record["sha_n"] = _sha(result.state.n)
            record["sha_p"] = _sha(result.state.p)
            # 只沿已收敛点推进，模拟 sweep 的热启动链。
            state = result.state
        records.append(record)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(records, indent=2), encoding="utf-8")
    print(f"wrote {len(records)} records to {out_path}")


def _diff(left_path: pathlib.Path, right_path: pathlib.Path) -> int:
    left = json.loads(left_path.read_text(encoding="utf-8"))
    right = json.loads(right_path.read_text(encoding="utf-8"))
    keys = ("converged", "iterations", "residual_norm", "update_norm", "reason",
            "sha_u", "sha_n", "sha_p")
    problems = 0
    for a, b in zip(left, right):
        mismatched = [key for key in keys if a.get(key) != b.get(key)]
        flag = "OK  " if not mismatched else "DIFF"
        print(f"{flag} {a['bias']:5.2f} V  {a['reason']}")
        if mismatched:
            problems += 1
            for key in mismatched:
                print(f"       {key}: {a.get(key)!r} != {b.get(key)!r}")
    if len(left) != len(right):
        problems += 1
        print(f"record count differs: {len(left)} vs {len(right)}")
    print("identical" if problems == 0 else f"{problems} bias point(s) differ")
    return 0 if problems == 0 else 1


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    run = sub.add_parser("run", help="run the Newton chain against one package tree")
    run.add_argument("--src", required=True, help="directory holding the `minitcad` package")
    run.add_argument("--out", required=True, type=pathlib.Path)

    diff = sub.add_parser("diff", help="compare two run outputs")
    diff.add_argument("left", type=pathlib.Path)
    diff.add_argument("right", type=pathlib.Path)

    args = parser.parse_args(argv)
    if args.command == "run":
        _run(args.src, args.out)
        return 0
    return _diff(args.left, args.right)


if __name__ == "__main__":
    raise SystemExit(main())
