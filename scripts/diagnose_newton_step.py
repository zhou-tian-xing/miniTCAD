"""诊断 Newton 的阶跃上限：为什么从 0.05 V 直接跳到 0.10 V 会失败。

```bash
python scripts/diagnose_newton_step.py              # 默认 0.05 -> 0.10 V
python scripts/diagnose_newton_step.py --from 0.05 --to 0.08
```

打印目标偏压上的真实求解结果，并复现 Newton 的第一步（按目标偏压投影接触边界），
看残差落在哪些结点。[5] 是对照组：同样的边界台阶也出现在能收敛的跃迁里。

输出全 ASCII，避免 Windows 控制台编码问题。
"""

from __future__ import annotations

import argparse
import pathlib
import sys


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--from", dest="src_V", type=float, default=0.05)
    parser.add_argument("--to", dest="dst_V", type=float, default=0.10)
    parser.add_argument("--config", default="examples/pn_baseline_0_to_030.yaml")
    args = parser.parse_args(argv)

    import numpy as np

    from minitcad.core.contact import contact_targets_qf
    from minitcad.core.device import build_device
    from minitcad.core.state import QFState
    from minitcad.discretization.residual import evaluate_qf
    from minitcad.io.config import load_config
    from minitcad.simulation.sweep import equilibrium_start_state
    from minitcad.solvers import solve_gummel, solve_newton
    from minitcad.solvers.convergence import residual_scales_qf
    from minitcad.solvers.newton import _project_contacts

    config = load_config(pathlib.Path(args.config))
    device = build_device(config)
    equilibrium = equilibrium_start_state(device)

    # 源偏压用一个能收敛的解；这里用 Gummel，避免把 Newton 自己的失败混进初值。
    source = solve_gummel(device, args.src_V, equilibrium)
    if not source.converged:
        print(f"source bias {args.src_V} V did not converge ({source.reason}); pick another --from")
        return 1
    print(f"[1] source solution at {args.src_V:g} V via gummel:")
    print(f"    max |jp_edges| = {np.abs(source.physics.fluxes.jp_edges).max():.3e}")

    # 目标偏压：Newton 从同一个源解出发
    result = solve_newton(device, args.dst_V, source.state)
    print(f"[2] solve_newton({args.dst_V:g} V) from that state:")
    print(f"    converged={result.converged} iterations={result.iterations} reason={result.reason}")
    print(f"    residual_norm={result.residual_norm:.3e} update_norm={result.update_norm:.3e}")

    # Newton 的第一步：按目标偏压的接触边界投影，检查残差落在哪些结点
    targets = device.contact_targets(args.dst_V)
    ref_n = float(np.log(targets[4]) - targets[3])
    ref_p = float(-np.log(targets[2]) - targets[0])
    qf = _project_contacts(
        QFState.from_physical(source.state, ref_n=ref_n, ref_p=ref_p),
        contact_targets_qf(device, args.dst_V, ref_n=ref_n, ref_p=ref_p),
    )
    physical = qf.to_physical()
    evaluation = evaluate_qf(state=qf, device=device, bias_V=args.dst_V)
    scales = residual_scales_qf(qf, device, args.dst_V, evaluation, physical).reshape(-1, 3)
    rows = evaluation.f.reshape(-1, 3)
    scaled = np.abs(rows) / scales

    print(f"[3] residual right after the contact projection to {args.dst_V:g} V:")
    print(f"    max |jp_edges| = {np.abs(evaluation.fluxes.jp_edges).max():.3e}   <-- was 1e-3 above")
    print(f"    node 1 raw rows  (u, eta_n, eta_p) = {np.array2string(rows[1], precision=4)}")
    print(f"    node 1 scales                      = {np.array2string(scales[1], precision=4)}")
    print(f"    node 1 scaled rows                 = {np.array2string(scaled[1], precision=4)}")
    interior = np.concatenate([scaled[2:-1].ravel(), scaled[0], scaled[-1]])
    print(f"    max scaled residual elsewhere      = {interior.max():.3e}")
    print("    -> the misbehaviour is localized at node 1, and there |f| == scale,")
    print("       because the row scale is computed from the spurious flux itself.")

    print("[4] same target, smaller jump (bisect the step):")
    for fraction in (0.5, 0.25, 0.125):
        middle = args.src_V + fraction * (args.dst_V - args.src_V)
        step = solve_newton(device, middle, source.state)
        print(f"    {args.src_V:g} -> {middle:.6g} V : converged={step.converged} reason={step.reason}")

    print("[5] control: the SAME signature also shows up in a jump that converges")
    print("    (so the boundary layer alone is not what makes it fail)")

    def signature(start_state, target_V):
        tgt = device.contact_targets(target_V)
        n_ref = float(np.log(tgt[4]) - tgt[3])
        p_ref = float(-np.log(tgt[2]) - tgt[0])
        trial = _project_contacts(
            QFState.from_physical(start_state, ref_n=n_ref, ref_p=p_ref),
            contact_targets_qf(device, target_V, ref_n=n_ref, ref_p=p_ref),
        )
        body = trial.to_physical()
        evaluation = evaluate_qf(state=trial, device=device, bias_V=target_V)
        row_scales = residual_scales_qf(
            trial, device, target_V, evaluation, body
        ).reshape(-1, 3)
        return (
            float(np.abs(evaluation.fluxes.jp_edges).max()),
            np.abs(evaluation.f).reshape(-1, 3)[1] / row_scales[1],
        )

    for label, start, target_V in (
        ("0.00 -> 0.05 V (converges)", equilibrium, 0.05),
        (f"{args.src_V:g} -> {args.dst_V:g} V (fails)", source.state, args.dst_V),
    ):
        peak, node1 = signature(start, target_V)
        print(f"    {label:28s} max |jp| = {peak:.2e}  node 1 scaled = {np.array2string(node1, precision=3)}")
    print("    -> identical boundary signature, different outcome: the layer is not the criterion.")

    print("[6] control: a 0.05 V ladder where every jump has the SAME delta v")
    print("    (so 'the step is too big' cannot explain the failures either)")
    from minitcad.core.device import SolverConfig
    from minitcad.simulation.dc import solve_single_bias

    newton = SolverConfig(method="newton")
    gummel = SolverConfig(method="gummel")
    state = equilibrium
    for step in range(1, 7):
        target = round(step * 0.05, 12)
        tried, _ = solve_single_bias(device, target, state, newton)
        peak, node1 = signature(state, target)
        print(
            f"    target {target:5.2f} V  jump {0.05:.2f}  "
            f"{'PASS' if tried.converged else 'FAIL':4s} it={tried.iterations:2d}  "
            f"max |jp| = {peak:.2e}"
        )
        accepted, _ = solve_single_bias(device, target, state, gummel)
        state = accepted.state
    print("    -> same delta v, same signature: 0.10/0.15 fail while 0.20/0.25/0.30 pass.")
    print("       The cause is NOT the step size and NOT the boundary layer. Still open.")
    return 0


if __name__ == "__main__":
    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))
    raise SystemExit(main())
