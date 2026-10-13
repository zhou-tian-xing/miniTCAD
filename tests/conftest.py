"""让 ``pytest`` 在未执行 ``pip install -e .`` 时也能导入 ``minitcad``（src 布局）。"""

from __future__ import annotations

import pathlib
import sys

_SRC = pathlib.Path(__file__).resolve().parents[1] / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))
