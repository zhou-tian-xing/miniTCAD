"""miniTCAD：一维硅 PN 结漂移扩散模拟。

内部计算全部用无量纲量；cm 是默认工程单位，SI 输入输出都经 `minitcad.io.units` 换算。
公开入口在 `w07-core-freeze` 冻结，见 `docs/api_reference.md`。
"""

from minitcad.api import run, write_required_outputs
from minitcad.io.config import config_fingerprint, load_config, parse_config_dict, save_config

__all__ = [
    "run",
    "write_required_outputs",
    "load_config",
    "parse_config_dict",
    "save_config",
    "config_fingerprint",
]
