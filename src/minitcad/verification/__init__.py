"""独立验证工具：不借用被验证对象。

本包可以 import `discretization.residual`，但**不允许** import `solvers.jacobian`，
否则 T10 就变成“用实现检验实现”。`tests/test_t10_jacobian.py` 直接读源码检查这一点。
"""
