"""一次 run 的结果容器（W3 最小版）。

字段依据第 10.2 节：attempts / successful / decision / failure_reason / config / device，
单点尝试保留 bias_V / status / result / runtime_s。W4 起由 sweep 填充。
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite

DECISION_ACCEPTED = "ACCEPTED"
DECISION_REJECTED = "REJECTED"
DECISION_CANCELLED = "CANCELLED"
DECISIONS = (DECISION_ACCEPTED, DECISION_REJECTED, DECISION_CANCELLED)

STATUS_PASS = "PASS"
STATUS_FAIL = "FAIL"
STATUS_CANCELLED = "CANCELLED"
STATUSES = (STATUS_PASS, STATUS_FAIL, STATUS_CANCELLED)


@dataclass(frozen=True, slots=True)
class Attempt:
    """一个已尝试的偏压点。result 为 SolveResult；W4 前可以为 None。"""

    bias_V: float
    status: str
    runtime_s: float = 0.0
    result: object | None = None

    def __post_init__(self) -> None:
        if not isfinite(self.bias_V):
            raise ValueError("bias_V must be finite")
        if self.status not in STATUSES:
            raise ValueError(f"status must be one of {STATUSES}, got {self.status!r}")
        if not isfinite(self.runtime_s) or self.runtime_s < 0.0:
            raise ValueError("runtime_s must be finite and >= 0")


@dataclass(frozen=True, slots=True)
class RunResult:
    """一次完整扫描的结果；GUI/CLI/writer 都只消费这个对象。"""

    config: object
    device: object
    attempts: tuple[Attempt, ...] = ()
    successful: tuple[Attempt, ...] = ()
    decision: str = DECISION_REJECTED
    failure_reason: str | None = None

    def __post_init__(self) -> None:
        if self.decision not in DECISIONS:
            raise ValueError(f"decision must be one of {DECISIONS}, got {self.decision!r}")
        if any(attempt.status != STATUS_PASS for attempt in self.successful):
            raise ValueError("successful may only contain PASS attempts")
        if not set(map(id, self.successful)).issubset(map(id, self.attempts)):
            raise ValueError("successful must be a subset of attempts")

    @property
    def accepted(self) -> bool:
        return self.decision == DECISION_ACCEPTED
