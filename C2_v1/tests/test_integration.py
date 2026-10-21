"""W2→W3 A/B 契约：仅以 core.Mesh1D/State/NormalizedDevice 串联完整离散,vibecoding"""
from dataclasses import replace
from types import SimpleNamespace
import numpy as np
import pytest
from minitcad.core.contact import contact_targets,InvalidPnOrientation
from minitcad.core.device import (
    DeviceConfig,MeshSpec,NormalizationConfig,build_device,
)
from minitcad.core.mesh import Mesh1D,uniform,segmented
from minitcad.core.material import SILICON
from minitcad.core.region import Region,DopingProfile
from minitcad.core.state import State,Fluxes
from minitcad.discretization.sg import edge_flux
from minitcad.discretization.continuity import continuity_residual_nd
from minitcad.discretization.residual import residual
from minitcad.physics.recombination import srh_recombination_nd
def sample_config(mesh=None):
    return DeviceConfig(
        length_cm=1e-3,area_cm2=1e-4,temperature_k=300.0,
        material=SILICON,
        regions=(
            Region(0.0,2.5e-4,"silicon",1e16,0.0),
            Region(2.5e-4,1e-3,"silicon",0.0,2e15),
        ),
        mesh=mesh or MeshSpec(kind="uniform",nodes=201),
        normalization=NormalizationConfig(reference_diffusivity_cm2_s=10.0),
    )
def test_B_device_and_contact_polarity():
    d=build_device(sample_config())
    assert d.mesh.node_count==201
    assert d.cbar_nd.shape==(199,)
    assert d.net_doping_nd.shape==(201,)
    assert d.net_doping_left_nd<0 and d.net_doping_right_nd>0
    targets0=contact_targets(d,0.)
    targets1=contact_targets(d,0.1)
    assert targets1[0]-targets0[0]==pytest.approx(0.1/d.thermal_voltage_v,abs=1e-12)
    assert targets1[3:]==targets0[3:]
    assert targets1[1:3]==targets0[1:3]
    assert (targets0[1]*targets0[2])==pytest.approx(1.,rel=1e-10)
    assert (targets0[4]*targets0[5])==pytest.approx(1.,rel=1e-10)
    assert targets0[0]<targets0[3]
def test_reversed_orientation_REJECTED():
    cfg=sample_config()
    with pytest.raises(InvalidPnOrientation,match="INVALID_PN_ORIENTATION"):
        replace(cfg,regions=tuple(replace(r,acceptor_cm3=r.donor_cm3,
                                         donor_cm3=r.acceptor_cm3) for r in cfg.regions))
    with pytest.raises(ValueError,match="grounded"):
        replace(cfg,right_voltage_V=0.1)
    with pytest.raises(ValueError,match="generic_1d"):
        replace(cfg,device_type="generic_1d")
def test_CV_average_at_step_not_point_sampling():
    profile=sample_config().doping_profile
    mesh=uniform(1e-3,201,1e-4)
    cv=profile.cv_averaged_net_doping_cm3(mesh.x_nodes_cm,mesh.h_edges_cm)
    junction_index=50
    assert mesh.x_nodes_cm[junction_index]==pytest.approx(2.5e-4)
    assert cv[junction_index-1]==pytest.approx((-1e16+2e15)/2)
    assert profile.net_doping_cm3([2.5e-4])[0]==2e15
    # 结点恰好在跃变处，其控制体平均不等于任一侧的浓度。
def test_nonuniform_segment_mesh_and_cv_invariants():
    mesh_spec=MeshSpec('segmented',breakpoints_cm=(0.,2.5e-4,3e-4,1e-3),
                       per_segment_nodes=(26,51,42))
    device=build_device(sample_config(mesh_spec))
    h=device.mesh.h_edges_nd
    np.testing.assert_allclose(device.mesh.w_interior_nd,0.5*(h[:-1]+h[1:]),rtol=1e-13)
    assert device.mesh.node_count==117
    with pytest.raises(ValueError,match="include junction"):
        build_device(sample_config(MeshSpec('segmented',breakpoints_cm=(0.,3e-4,1e-3),
                                           per_segment_nodes=(10,10))))
    with pytest.raises(ValueError,match="integer"):
        uniform(1,3.9,0.1)
    with pytest.raises(ValueError,match="integer"):
        segmented(1,[0,1],[3.8],0.1)
def test_state_pack_copy_and_only_readonly_arrays():
    u=np.array([-1.,0.,1.]);n=np.array([0.2,1.,5.]);p=np.array([5.,1.,0.2])
    state=State(u,n,p)
    x=state.pack()
    assert x.tolist()==[-1.,0.2,5.,0.,1.,1.,1.,5.,0.2]
    repacked=State.unpack(x,3)
    np.testing.assert_array_equal(repacked.pack(),x)
    u[0]=3.5;n[0]=99.0;x[1]=-1.0
    assert state.u[0]==-1.0 and state.n[0]==0.2
    with pytest.raises(ValueError):
        state.u.setflags(write=True)
    with pytest.raises(ValueError):
        State(u,[1.,0.,1.],p)
    with pytest.raises(ValueError):
        State.unpack(np.ones(8),3)
def test_W3_F_layout_and_all_six_BC_rows():
    device=build_device(sample_config(MeshSpec(kind='uniform',nodes=11)))
    u=np.linspace(*[contact_targets(device,0.)[i] for i in (0,3)],11)
    state=State(u=np.asarray(u),n=np.exp(u),p=np.exp(-u))
    f=residual(state,device,0.0)
    assert f.shape==(3*11,) and np.isfinite(f).all()
    expected=np.array(contact_targets(device,0.))
    np.testing.assert_allclose(f[[0,1,2,30,31,32]],
                               [state.u[0]-expected[0],state.n[0]-expected[1],
                                state.p[0]-expected[2],state.u[-1]-expected[3],
                                state.n[-1]-expected[4],state.p[-1]-expected[5]],
                               atol=1e-8,rtol=0)
    # 这是构造的Boltzmann初态，并非已收敛Poisson解；绝不宣称Gate 1完成。
    flux=edge_flux(state,device)
    assert np.max(abs(flux.jn_edges))<1e-4
    assert np.max(abs(flux.jp_edges))<1e-4
    assert np.max(abs(f[4:-3:3]))<1e-4
    assert np.max(abs(f[5:-3:3]))<1e-4
    new=residual(state,device,0.1)
    assert new[0]==pytest.approx(-0.1/device.thermal_voltage_v,abs=1e-12)
    np.testing.assert_allclose(new[1:3],f[1:3],atol=0,rtol=0)
    np.testing.assert_allclose(new[-3:],f[-3:],atol=0,rtol=0)
def test_exact_discrete_current_identity_with_nonzero_R():
    device=build_device(sample_config(MeshSpec(kind='uniform',nodes=9)))
    m=device.mesh.node_count
    state=State(np.linspace(-2.,2.,m),np.linspace(0.1,2,m),np.linspace(3,0.2,m))
    flux=edge_flux(state,device)
    rate=srh_recombination_nd(state.n,state.p,device.tau_n_nd,device.tau_p_nd)
    fn,fp=continuity_residual_nd(flux,rate,device.mesh.w_interior_nd)
    full=residual(state,device,0.)
    np.testing.assert_allclose(full[4:-3:3],fn,rtol=0,atol=0)
    np.testing.assert_allclose(full[5:-3:3],fp,rtol=0,atol=0)
    current_div=np.diff(flux.jn_edges+flux.jp_edges)/device.mesh.w_interior_nd
    np.testing.assert_allclose(fn+fp,current_div,rtol=1e-12,atol=1e-12)
    bad=State(np.ones(4),np.ones(4),np.ones(4))
    with pytest.raises(ValueError,match="node_count"):
        residual(bad,device,0.)
    with pytest.raises(ValueError,match="shape"):
        continuity_residual_nd(flux,np.ones(m-2),device.mesh.w_interior_nd)
def test_T03_contact_neutrality_and_mass_action():
    d=build_device(sample_config())
    for doping,n,p in ((d.net_doping_left_nd,*contact_targets(d,0.)[1:3]),
                       (d.net_doping_right_nd,*contact_targets(d,0.)[4:6])):
        assert abs((doping+p-n)/max(1.,abs(doping)))<=1e-10
        assert abs(n*p-1.)<=1e-10
    with pytest.raises(ValueError):
        d.cbar_nd.setflags(write=True)
    with pytest.raises(ValueError):
        d.mesh.x_nodes_nd.setflags(write=True)
def test_three_node_manual_coupled_residual_with_nonzero_SRH():
    device=build_device(sample_config(MeshSpec(kind='uniform',nodes=3)))
    h=device.mesh.h_edges_nd
    w=device.mesh.w_interior_nd[0]
    state=State(u=[0.,0.,0.],n=[2.,4.,7.],p=[2.,3.,5.])
    taup=device.tau_p_nd; taun=device.tau_n_nd
    source=(4.*3.-1.)/(taup*(4.+1.)+taun*(3.+1.))
    jn_left=device.mu_n_nd*(4.-2.)/h[0]
    jn_right=device.mu_n_nd*(7.-4.)/h[1]
    jp_left=-device.mu_p_nd*(3.-2.)/h[0]
    jp_right=-device.mu_p_nd*(5.-3.)/h[1]
    Fu=device.cbar_nd[0]+state.p[1]-state.n[1]
    Fn=(jn_right-jn_left)/w-source
    Fp=(jp_right-jp_left)/w+source
    got=residual(state,device,0.)
    np.testing.assert_allclose(got[3:6],[Fu,Fn,Fp],rtol=1e-12,atol=1e-12)
    assert got.shape==(9,)