"""W2-B4：cm <-> SI 双向适配与字段后缀契约（第 3.3 节）。"""

from __future__ import annotations

import pytest

from minitcad.io.units import (
    CM,
    SI,
    UnitError,
    UnitSystem,
    cm_to_si,
    convert,
    require_suffixed_field,
    si_to_cm,
)

SAMPLES = {
    "length": 1.0e-3,
    "area": 1.0e-4,
    "density": 1.0e16,
    "mobility": 1350.0,
    "diffusivity": 10.0,
    "electric_field": 1.0e4,
    "current_density": 1.0e-2,
    "recombination_rate": 1.0e21,
    "charge_density": 1.0e-3,
    "permittivity": 1.03594e-12,
    "potential": 0.6,
    "current": 1.0e-3,
    "lifetime": 1.0e-6,
    "temperature": 300.0,
}


@pytest.mark.parametrize("quantity,value", sorted(SAMPLES.items()))
def test_cm_to_si_to_cm_round_trip(quantity, value):
    assert si_to_cm(quantity, cm_to_si(quantity, value)) == pytest.approx(value, rel=1e-12)


@pytest.mark.parametrize("quantity,value", sorted(SAMPLES.items()))
def test_si_to_cm_to_si_round_trip(quantity, value):
    assert cm_to_si(quantity, si_to_cm(quantity, value)) == pytest.approx(value, rel=1e-12)


def test_mobility_factor_is_the_documented_one():
    assert cm_to_si("mobility", 1350.0) == pytest.approx(0.135)
    assert cm_to_si("density", 1.0e16) == pytest.approx(1.0e22)
    assert cm_to_si("diffusivity", 10.0) == pytest.approx(1.0e-3)


def test_same_system_conversion_is_identity():
    assert convert("density", 1.0e16, CM, CM) == 1.0e16


def test_unknown_quantity_and_non_finite_are_rejected():
    with pytest.raises(UnitError):
        cm_to_si("mass", 1.0)
    with pytest.raises(UnitError):
        cm_to_si("density", float("inf"))


def test_units_system_parsing():
    assert UnitSystem.parse("CM").name == CM
    assert UnitSystem.parse(" si ").name == SI
    with pytest.raises(UnitError):
        UnitSystem.parse("mks")


def test_suffixed_field_accepts_matching_system():
    assert require_suffixed_field("length_cm", CM) == "_cm"
    assert require_suffixed_field("intrinsic_density_m3", SI) == "_m3"


def test_suffixed_field_rejects_mixed_systems():
    with pytest.raises(UnitError):
        require_suffixed_field("density_m3", CM)
    with pytest.raises(UnitError):
        require_suffixed_field("length_cm", SI)


def test_unsuffixed_field_is_rejected():
    # 无单位的 density=1e16 一律拒绝，不得按数值大小猜测单位（第 3.3 节）。
    with pytest.raises(UnitError):
        require_suffixed_field("density", CM)
    with pytest.raises(UnitError):
        require_suffixed_field("density", SI)
