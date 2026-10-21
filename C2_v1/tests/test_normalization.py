#第2周 cm→无量纲（D0显式）；不同T只接受该T的材料参数
from dataclasses import replace
import numpy as np
import pytest

from minitcad.core.constants import q,k_b_j_per_k
from minitcad.core.material import SILICON, MaterialTemperatureMismatch
from minitcad.io.units import convert_units
from minitcad.physics.normalization import (
    build_normalization_scales,normalize_cm,denormalize_cm,
)


def test_cm_default_scales_and_D0_policy():
    scale=build_normalization_scales(SILICON,300.0,reference_diffusivity_cm2_s=10.0)
    vt=k_b_j_per_k*300.0/q
    assert scale.v0_v==pytest.approx(vt,rel=1e-14)
    assert scale.l0_cm==pytest.approx(4.08845e-3,rel=1e-5)
    assert scale.d0_cm2_s==10.0
    assert scale.mu0_cm2_per_vs==pytest.approx(10.0/vt,rel=1e-14)
    assert scale.j0_a_cm2==pytest.approx(q*10.0*SILICON.intrinsic_density_cm3/scale.l0_cm)
    assert scale.r0_cm3_s==pytest.approx(10.0*SILICON.intrinsic_density_cm3/scale.l0_cm**2)
    m0=1000.0
    scales_mu=build_normalization_scales(SILICON,300.,reference_mobility_cm2_per_Vs=m0)
    assert scales_mu.d0_cm2_s==pytest.approx(m0*vt,rel=1e-14)
    assert scales_mu.d0_cm2_s!=10.0


def test_reference_T_and_invalid_D0():
    material350=replace(SILICON,reference_temperature_k=350.0,
                        intrinsic_density_cm3=2e10,mobility_n_cm2_per_Vs=1200.0)
    s=build_normalization_scales(material350,350,reference_diffusivity_cm2_s=10.0)
    assert s.v0_v==pytest.approx(k_b_j_per_k*350/q)
    assert s.n0_cm3==2e10
    with pytest.raises(MaterialTemperatureMismatch):
        build_normalization_scales(SILICON,350,reference_diffusivity_cm2_s=10.0)
    with pytest.raises(ValueError):
        build_normalization_scales(SILICON,300)
    with pytest.raises(ValueError):
        build_normalization_scales(SILICON,300,reference_diffusivity_cm2_s=10.0,
                                   reference_mobility_cm2_per_Vs=1000.0)
    with pytest.raises(ValueError):
        build_normalization_scales(SILICON,300,reference_diffusivity_cm2_s=-10.0)


@pytest.mark.parametrize("name,value",[
    ("x_cm",[1e-5,2e-4]),("phi_V",[-0.5,0.0,0.3]),
    ("density_cm3",[1e10,2e16]),("net_doping_cm3",[-1e16,2e15]),
    ("mobility_cm2_per_Vs",[1350.,480.]),("diffusivity_cm2_s",[1.,10.]),
    ("E_V_cm",[0.,1e4]),("J_A_cm2",[-1e-2,2e-2]),
    ("R_cm3_s",[0.,1e20]),("tau_s",[1e-6,2e-6]),
])
def test_cm_dimensionless_roundtrip(name,value):
    scales=build_normalization_scales(SILICON,300,reference_diffusivity_cm2_s=10.0)
    nd=normalize_cm(value,name,scales)
    recovered=denormalize_cm(nd,name,scales)
    np.testing.assert_allclose(recovered,value,rtol=3e-15,atol=0)


@pytest.mark.parametrize("kind,x",[
    ("length",[1e-3,1e-5]),("area",[1e-4,1.5]),
    ("density",[1e10,1e16]),("mobility",[1350.,480.]),
    ("diffusivity",[10.,100.]),("electric_field",[1e3,1e6]),
    ("current_density",[0.,0.03]),("recombination",[3e8,1e21]),
    ("charge_density",[1e-9,1e-4]),("permittivity",[1e-12,2e-12]),
    ("voltage",[-1.,0.,1.]),("current",[1e-4,2e-5]),
    ("lifetime",[1e-6,2e-6]),("temperature",[300.,350.]),
])
def test_cm_si_roundtrip_both_directions(kind,x):
    cm=np.array(x)
    si=convert_units(cm,kind,'cm','si')
    np.testing.assert_allclose(convert_units(si,kind,'si','cm'),cm,rtol=3e-15,atol=0)
    np.testing.assert_allclose(convert_units(convert_units(si,kind,'si','cm'),kind,'cm','si'),si,rtol=3e-15,atol=0)


def test_SI_adapt_cm_normalize_restore_SI():
    scales=build_normalization_scales(SILICON,300,reference_diffusivity_cm2_s=10.0)
    cases=[("length","x_cm",0.003),("density","density_cm3",1e16),
           ("mobility","mobility_cm2_per_Vs",1350.),
           ("diffusivity","diffusivity_cm2_s",10.),
           ("electric_field","E_V_cm",300.),
           ("current_density","J_A_cm2",0.001)]
    for kind,qname,value in cases:
        si=convert_units(value,kind,'cm','si')
        from_si=convert_units(si,kind,'si','cm')
        nd=normalize_cm(from_si,qname,scales)
        back_cm=denormalize_cm(nd,qname,scales)
        np.testing.assert_allclose(convert_units(back_cm,kind,'cm','si'),si,rtol=4e-15)
    with pytest.raises(ValueError):
        convert_units(1.,'density','unknown','cm')
    with pytest.raises(ValueError):
        convert_units(1.,'nonsense','cm','si')


def test_einstein_relation_same_T_and_single_D0():
    scales=build_normalization_scales(SILICON,300,reference_diffusivity_cm2_s=10.0)
    mn=SILICON.mobility_n_cm2_per_Vs/scales.mu0_cm2_per_vs
    mp=SILICON.mobility_p_cm2_per_Vs/scales.mu0_cm2_per_vs
    assert mn*scales.d0_cm2_s==pytest.approx(SILICON.mobility_n_cm2_per_Vs*scales.v0_v)
    assert mp*scales.d0_cm2_s==pytest.approx(SILICON.mobility_p_cm2_per_Vs*scales.v0_v)
