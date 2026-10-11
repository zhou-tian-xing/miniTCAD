#物理量测试
#python -m pytest
import pytest
import materials
import PhysicsBase
#测试热电压
def test_thermal_voltage_300k():
    vt1 = PhysicsBase.thermal_voltage(300)
    vt2 = PhysicsBase.thermal_voltage(600)
    assert vt1 == pytest.approx(
        0.025852,
        rel = 1e-6,
        abs = 0
    )
    assert vt2 == pytest.approx(
        2 * vt1,
        rel = 1e-6,
        abs = 0
    )
#热电压的温度异常值输入
@pytest.mark.parametrize(
        "Bad_Temperature",
        [0,-10.0,-100.0]
)
def test_nonpositive_temperature(Bad_Temperature):
    with pytest.raises(ValueError):
        PhysicsBase.thermal_voltage(Bad_Temperature)
#硅的介电常数
def test_silicon_permittivity():
    epsilon = materials.SILICON.epsilon_s
    assert epsilon == pytest.approx(
        1.0359399740976e-12,
        rel = 1e-8,
        abs = 0
    )
#参考ni模型的温度异常输入
def test_reference_ni_rejects_other_temperature():
    model = PhysicsBase.ReferenceIntrinsicCarrierModel()
    with pytest.raises(ValueError):
        model.intrinsic_concentration(materials.SILICON,400)
#Einstein关系
def test_einstein_relation():
    mobility = PhysicsBase.ConstantMobilityModel().evaluate(
        materials.SILICON,300.0
    )
    diffusion = PhysicsBase.einstein_diffusion(mobility,300.00)
    vt = PhysicsBase.thermal_voltage(300)
    assert diffusion.electron_cm2_per_s == pytest.approx(
        mobility.electron_cm2_per_vs * vt,
        rel = 1e-12,
        abs = 0.0
    )
    assert diffusion.hole_cm2_per_s == pytest.approx(
        mobility.hole_cm2_per_vs * vt,
        rel = 1e-12,
        abs = 0.0
    )
    assert diffusion.electron_cm2_per_s == pytest.approx(
        36.1928, rel = 1e-5, abs = 0.0
    )
    assert diffusion.hole_cm2_per_s == pytest.approx(
        11.6334, rel = 1e-5, abs = 0.0
    )
#本征Debye长度
def test_intrinsic_debye_length():
    ni_model = PhysicsBase.ReferenceIntrinsicCarrierModel()
    scales = PhysicsBase.build_normalization_scales(
        materials.SILICON,300,ni_model
    )
    assert scales.l0_cm == pytest.approx(
        4.08845e-3,
        rel = 1e-5,
        abs = 0.0,
    )
#SRH复合:平衡
def test_srh_equilibrium():
    ni = 1.0e10
    parameters = PhysicsBase.SRHParameters(
        tau_n_s = 1e-6,
        tau_p_s = 2e-6,
        trap_offset_ev = 0.0
    )
    rate = PhysicsBase.srh_recombination_rate(
        n_per_cm3 = 2.0 * ni,
        p_per_cm3 = 0.5 * ni,
        ni_per_cm3 = ni,
        temperature_k = 300.0,
        params = parameters
    )
    assert rate == pytest.approx(
        0,0,abs = 1e-12
    )
#SRH复合：非平衡
def test_srh_nonequilibrium_sign():
    ni = 1.0e10
    parameters = PhysicsBase.SRHParameters(
        tau_n_s = 1e-6,
        tau_p_s = 2e-6,
        trap_offset_ev = 0.0
    )
    recombination = PhysicsBase.srh_recombination_rate(
        2.0 * ni, ni, ni, 300.0, parameters
    )
    generation = PhysicsBase.srh_recombination_rate(
        0.5 * ni, ni, ni, 300.0, parameters
    )
    assert recombination == pytest.approx(
        1.25e15, rel=1e-12, abs=0.0
    )
    assert generation == pytest.approx(
        -1.0e15, rel=1e-12, abs=0.0
    )
#SRH的异常输入
@pytest.mark.parametrize(
    "n,p,ni",
    [
        (-1.0, 1e10, 1e10),
        (1e10, -1.0, 1e10),
        (1e10, 1e10, 0.0),
    ]
)
def test_srh_rejects_invalid_density(n,p,ni):
    parameters = PhysicsBase.SRHParameters(
        tau_n_s = 1e-6,
        tau_p_s = 2e-6,
        trap_offset_ev = 0.0,
    )
    with pytest.raises(ValueError):
        PhysicsBase.srh_recombination_rate(
            n,p,ni,300.0,parameters
        )