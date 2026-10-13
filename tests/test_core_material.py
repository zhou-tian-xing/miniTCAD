"""W2-B2：cm-first 材料与温度策略（第 2.1、2.4 节）。"""

from __future__ import annotations

import pytest

from minitcad.core.material import (
    POLICY_EXPLICIT_AT_T,
    SILICON,
    Material,
    MaterialTemperatureMismatch,
    get_material,
    material_parameter_policy,
)


def test_silicon_absolute_permittivity_matches_appendix_d2():
    # 附录 D.2：ε_s ≈ 1.0359399740976e-12 F/cm
    assert SILICON.epsilon_f_per_cm == pytest.approx(1.0359399740976e-12, rel=1e-9)


def test_get_material_is_case_insensitive():
    assert get_material("Si") is SILICON
    assert get_material(" silicon ") is SILICON


def test_unknown_material_is_rejected():
    with pytest.raises(ValueError):
        get_material("germanium")


def test_nonpositive_parameter_is_rejected():
    with pytest.raises(ValueError):
        Material(
            name="bad",
            symbol="X",
            epsilon_relative=11.7,
            reference_temperature_k=300.0,
            intrinsic_density_cm3=0.0,  # 非正 n_i
            mobility_n_cm2_per_Vs=1400.0,
            mobility_p_cm2_per_Vs=450.0,
            lifetime_n_s=1e-6,
            lifetime_p_s=2e-6,
            bandgap_eV=1.12,
        )


def test_temperature_mismatch_is_rejected_not_silently_fixed():
    # 第 2.4 节：标定温度 300 K 的参数不得用于 400 K 的运行。
    with pytest.raises(MaterialTemperatureMismatch):
        material_parameter_policy(SILICON, 400.0)


def test_temperature_match_returns_policy_tag():
    assert material_parameter_policy(SILICON, 300.0) == POLICY_EXPLICIT_AT_T


def test_nonpositive_run_temperature_is_rejected():
    with pytest.raises(ValueError):
        material_parameter_policy(SILICON, 0.0)
