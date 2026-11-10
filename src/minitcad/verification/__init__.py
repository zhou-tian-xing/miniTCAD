"""独立验证工具：不从生产代码里借用任何“被验证对象”。

`verification` 与 `solvers` 的依赖方向是单向的：这里可以 import `discretization.residual`
（被微分的残差），但**不允许** import `solvers.jacobian`——否则 T10 就变成“用 A 的实现
检验 A 的实现”。`tests/test_t10_jacobian.py` 里有一条断言直接检查这一点（读源码文本），
因为 import 图本身挡不住后人加一行 import。
"""
