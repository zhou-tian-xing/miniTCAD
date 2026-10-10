#SRH复合
from dataclasses import dataclass
from math import exp, isfinite
from constants import kbv
@dataclass(frozen=True, slots=True)
class SRHParameters:
    tau_n_s: float
    tau_p_s: float
    trap_offset_ev: float  # E_t - E_i, in eV
    def __post_init__(self) -> None:
        if not isfinite(self.tau_n_s) or self.tau_n_s <= 0.0:
            raise ValueError("tau_n_s must be finite and > 0")
        if not isfinite(self.tau_p_s) or self.tau_p_s <= 0.0:
            raise ValueError("tau_p_s must be finite and > 0")
        if not isfinite(self.trap_offset_ev):
            raise ValueError("trap_offset_ev must be finite")
def srh_recombination_rate(
    n_per_cm3: float,
    p_per_cm3: float,
    ni_per_cm3: float,
    temperature_k: float,
    params: SRHParameters,
) -> float:
    if not isfinite(n_per_cm3) or n_per_cm3 < 0.0:
        raise ValueError("n_per_cm3 must be finite and >= 0")
    if not isfinite(p_per_cm3) or p_per_cm3 < 0.0:
        raise ValueError("p_per_cm3 must be finite and >= 0")
    if not isfinite(ni_per_cm3) or ni_per_cm3 <= 0.0:
        raise ValueError("ni_per_cm3 must be finite and > 0")
    if not isfinite(temperature_k) or temperature_k <= 0.0:
        raise ValueError("temperature_k must be finite and > 0")
    thermal_energy_ev = kbv * temperature_k
    if thermal_energy_ev == 0.0:
        raise ValueError("temperature_k is too small for this calculation")
    eta_t = params.trap_offset_ev / thermal_energy_ev
    try:
        n1 = ni_per_cm3 * exp(eta_t)
        p1 = ni_per_cm3 * exp(-eta_t)
    except OverflowError as exc:
        raise ValueError("trap energy causes exponential overflow") from exc
    numerator = n_per_cm3 * p_per_cm3 - ni_per_cm3**2
    denominator = (
        params.tau_p_s * (n_per_cm3 + n1)
        + params.tau_n_s * (p_per_cm3 + p1)
    )
    if not isfinite(numerator) or not isfinite(denominator) or denominator <= 0.0:
        raise ValueError("SRH intermediate result is not finite or valid")
    rate = numerator / denominator
    if not isfinite(rate):
        raise ValueError("SRH rate is not finite")
    return rate