# MiniTCAD 一维 PN — W3 物理离散接口版

本包是规范所定义 **W3 物理/离散残差阶段**，不是已完成的 DC 求解器。唯一规范是根目录随附的《MiniTCAD_1D_PN_8周实施验收规范_cm工程单位版_VSCode.md》。

## VS Code / 命令行

在仓库根目录用 Python 3.10+（含 NumPy、SciPy）运行：

```bash
python -m pip install -e ".[dev]"
python -m pytest -q
python -m examples.verify_poisson
```

VS Code 选择项目虚拟环境的 Python 解释器后，可使用自带的 **MiniTCAD: pytest** / **MiniTCAD: Poisson SI MMS** 任务。

Python 中调用物理离散核心（结果 `F` 为无量纲 `3M`，**不是求解后的收敛解**）：

```python
import numpy as np
from minitcad.core.device import DeviceConfig, MeshSpec, NormalizationConfig, build_device
from minitcad.core.material import SILICON
from minitcad.core.region import Region
from minitcad.core.state import State
from minitcad.discretization.residual import residual

cfg = DeviceConfig(
    length_cm=1e-3, area_cm2=1e-4, temperature_k=300.0,
    material=SILICON,
    regions=(Region(0, 2.5e-4, "silicon", 1e16, 0),
             Region(2.5e-4, 1e-3, "silicon", 0, 2e15)),
    mesh=MeshSpec(kind="uniform", nodes=201),
    normalization=NormalizationConfig(reference_diffusivity_cm2_s=10.0),
)
device = build_device(cfg)
u = np.linspace(*[device.contact_targets(0.0)[i] for i in (0, 3)],
                device.mesh.node_count)
state = State(u, np.exp(u), np.exp(-u))
F = residual(state, device, 0.0)
print(F.shape, np.isfinite(F).all())
```

参数和方向：**左 P，右 N；`V_app = V_left - V_right`；右端接地；`V_app>0` 为正向**。cm 工程量在 `build_device` 处归一化，求解核心不接收 SI 残差或物理电压电势数组。

## 交付边界

W4 Equilibrium/Gummel，W5 Newton/Jacobian，以及 B 的配置读写、CLI、完整 SI YAML 适配、DC sweep、固定三文件 writer 尚未实现。当前只可执行 W2/W3 模块验证，不能运行 `minitcad solve` 或宣称 Gate 1/Gate 2 通过。

详细逐文件记录：`reports/W3D2_修改报告.md`；A/B 接口确认表：`docs/a_b_interface_week02.md`。
