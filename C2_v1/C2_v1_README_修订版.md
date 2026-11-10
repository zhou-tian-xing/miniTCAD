# MiniTCAD 一维 PN — 数值求解开发版

> 文档状态：与 `C2_v1.zip` 当前源代码快照核对（2026-10-27）。**已有求解能力不等于通过项目 Gate 验收**。  
> 方向和极性：**左 P、右 N；右端接地；`V_app = V_left - V_right`；正向偏压 `V_app > 0`**。

## 当前已实现

| 模块 | 实现路径 | 状态 |
|---|---|---|
| 器件/材料/非均匀网格/掺杂控制体平均 | `src/minitcad/core/` | 已有代码及测试 |
| cm 归一化、SI 基础单位换算 | `physics/normalization.py`, `io/units.py` | 已有底层能力；**无完整 YAML 输入适配** |
| Poisson、SG、SRH、continuity | `src/minitcad/discretization/`, `physics/recombination.py` | 已实现 |
| W3 物理变量残差 `residual(State,...)` | `discretization/residual.py` | 保留的旧兼容接口，**不是当前求解器主路径** |
| 准费米残差与诊断 `evaluate_qf`, `residual_qf` | `discretization/residual.py` | 当前主路径 |
| 零偏压平衡态标量求解 | `simulation/equilibrium.py` | 已实现 |
| Gummel / Newton | `solvers/gummel.py`, `solvers/newton.py` | 已实现 |
| 九色分组稀疏差分 Jacobian | `solvers/jacobian.py` | 已实现 `jacobian_qf` |
| 反向偏压自适应扫描 | `simulation/reverse_sweep.py` | 已实现**反向限定版**，非通用 DC sweep |

在当前快照运行 `python -m pytest -q`：**195 passed**（只代表现有自动测试结果）。

## 安装与验证

在仓库根目录运行（Python 3.10+）：

```bash
python -m pip install -e ".[dev]"
python -m pytest -q
python -m examples.verify_poisson
```

VS Code 可选项目虚拟环境的 Python 解释器，使用 **MiniTCAD: pytest** 或 **MiniTCAD: Poisson SI MMS** 任务。

## 可运行示例：求零偏压平衡态并用 Newton 复核

```python
from minitcad.core.device import (
    DeviceConfig, MeshSpec, NormalizationConfig, build_device,
)
from minitcad.core.material import SILICON
from minitcad.core.region import Region
from minitcad.simulation.equilibrium import solve_equilibrium_scalar
from minitcad.solvers.newton import solve_newton

cfg = DeviceConfig(
    length_cm=1e-3,
    area_cm2=1e-4,
    temperature_k=300.0,
    material=SILICON,
    regions=(
        Region(0, 2.5e-4, "silicon", 1e16, 0),
        Region(2.5e-4, 1e-3, "silicon", 0, 2e15),
    ),
    mesh=MeshSpec(kind="uniform", nodes=201),
    normalization=NormalizationConfig(reference_diffusivity_cm2_s=10.0),
)
device = build_device(cfg)
u_left, _, _, u_right, _, _ = device.contact_targets(0.0)
eq = solve_equilibrium_scalar(
    device.mesh.x_nodes_nd,
    device.cbar_nd,
    u_left,
    u_right,
)
assert eq.converged and eq.state is not None, eq.reason

result = solve_newton(device=device, bias_V=0.0, initial_state=eq.state)
assert result.converged and result.state is not None, result.reason
print(result.reason, result.iterations, result.residual_norm)
```

`State(u,n,p)` 是公开的**物理状态（但数值已无量纲化）**。当前求解器内部使用 `QFState(u,eta_n,eta_p,ref_n,ref_p)`，并由 `evaluate_qf` 产生完整残差、边电流与复合率。不要把 `residual(State,...)` 误认为 Newton 当前的 Jacobian 对象。

材料温度要求：内置 `SILICON` 只给出与 **300 K** 对应的参数；运行时温度不匹配会报错。`NormalizationConfig` 必须在参考扩散系数、参考迁移率两种参数中**恰好选择一项**。

## 尚未完成 / 不能宣称通过

- 统一 YAML 解析和校验、完整 SI YAML 适配。
- `minitcad solve` / `validate` CLI、统一公共 `run()` / `RunResult` 工作流。
- 通用正反向 fixed/adaptive DC sweep（现有 `solve_reverse_sweep()` 只支持负向目标电压）。
- 必需的 `comparison.csv`、`comparison.png`、`convergence.md` 三文件 writer。
- PySide6 GUI、Windows 打包与正式发布。
- Gate 1、Gate 2、Gate 3 的**独立验收报告及签署证据**（当前压缩包未提供，不得从单测通过推断 Gate 已通过）。

## 文档导航

- [当前求解器接口](docs/current_solver_interfaces.md)：**开发时先看这一份**。
- [物理与归一化约定](physics.md)：物理理论说明 + 当前实现约束。
- [八周实施验收规范](MiniTCAD_1D_PN_8周实施验收规范_cm工程单位版_VSCode.md)：**项目目标/验收基线，不是已经完成的功能表**。
- [W2→W3 历史接口记录](docs/archive/a_b_interface_week02_history.md)：保留历史，不作为当前求解器接口。
- [本次文档变更清单](DOCS_CHANGELOG.md)。

**注意：** 原 README 引用的 `reports/W3D2_修改报告.md` 不在此压缩包中，因此已移除该失效链接。
