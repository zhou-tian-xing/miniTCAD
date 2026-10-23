"""QFState and production residual interface tests."""
#vibecoding

import numpy as np
import pytest

from minitcad.core.state import State, QFState

from minitcad.core.contact import (
    contact_targets_qf,
)

from minitcad.core.device import (
    build_device,
    DeviceConfig,
    MeshSpec,
    NormalizationConfig,
)

from minitcad.core.material import SILICON
from minitcad.core.region import Region

from minitcad.discretization.residual import (
    evaluate_qf,
    residual_qf,
    ResidualEvaluation,
)


def make_device():
    return build_device(
        DeviceConfig(
            length_cm=1e-3,
            area_cm2=1e-4,
            temperature_k=300.0,
            material=SILICON,

            # Region arguments:
            # start, end, material, acceptor, donor
            regions=(
                Region(
                    0.0, 5e-4,
                    "silicon",
                    1e11, 0.0,
                ),
                Region(
                    5e-4, 1e-3,
                    "silicon",
                    0.0, 2e11,
                ),
            ),

            mesh=MeshSpec(
                kind="uniform",
                nodes=31,
            ),

            normalization=NormalizationConfig(
                reference_diffusivity_cm2_s=10.0
            ),
        )
    )


def test_qf_roundtrip_and_pack():
    state = State(
        u=np.array([-10.0, 0.0, 10.0]),
        n=np.array([1e-5, 1.0, 1e5]),
        p=np.array([1e5, 1.0, 1e-5]),
    )

    qf = QFState.from_physical(
        state,
        ref_n=0.3,
        ref_p=-0.2,
    )

    restored = qf.to_physical()

    np.testing.assert_allclose(
        restored.u,
        state.u,
        rtol=0,
        atol=0,
    )

    np.testing.assert_allclose(
        restored.n,
        state.n,
        rtol=1e-14,
    )

    np.testing.assert_allclose(
        restored.p,
        state.p,
        rtol=1e-14,
    )

    qf2 = QFState.unpack_qf(
        qf.pack_qf(),
        3,
        ref_n=0.3,
        ref_p=-0.2,
    )

    np.testing.assert_array_equal(
        qf2.pack_qf(),
        qf.pack_qf(),
    )

    # Internal state arrays must be read-only
    assert not qf.u.flags.writeable
    assert not qf.eta_n.flags.writeable
    assert not qf.eta_p.flags.writeable


@pytest.mark.parametrize(
    "bias",
    [-0.05, 0.0, 0.05],
)
def test_qf_contact_sign(bias):
    device = make_device()

    v = bias / device.thermal_voltage_v

    qft = contact_targets_qf(
        device,
        bias,
    )

    # Actual QF levels:
    # left = -v, right = 0
    np.testing.assert_allclose(
        (
            qft[1],
            qft[2],
            qft[4],
            qft[5],
        ),
        (-v, -v, 0.0, 0.0),
        rtol=0,
        atol=1e-13,
    )

    # Recover the six physical contact values
    t = device.contact_targets(bias)

    reconstructed = [
        np.exp(qft[0] + qft[1]),
        np.exp(-qft[0] - qft[2]),
        np.exp(qft[3] + qft[4]),
        np.exp(-qft[3] - qft[5]),
    ]

    np.testing.assert_allclose(
        reconstructed,
        [t[1], t[2], t[4], t[5]],
        rtol=1e-14,
        atol=1e-14,
    )


def test_residual_evaluation_contract():
    device = make_device()
    m = device.mesh.node_count

    t = device.contact_targets(0.02)

    u = np.linspace(
        t[0],
        t[3],
        m,
    )

    eta_n = np.linspace(
        -0.1,
        0.02,
        m,
    )

    eta_p = np.linspace(
        -0.05,
        0.03,
        m,
    )

    qf = QFState(
        u=u,
        eta_n=eta_n,
        eta_p=eta_p,
    )

    ev = evaluate_qf(
        qf,
        device,
        0.02,
    )

    assert isinstance(ev, ResidualEvaluation)

    assert ev.f.shape == (3 * m,)
    assert ev.f.dtype == np.float64

    assert ev.fluxes.jn_edges.shape == (m - 1,)
    assert ev.fluxes.jp_edges.shape == (m - 1,)

    assert ev.recombination_nd.shape == (m,)

    np.testing.assert_array_equal(
        residual_qf(qf, device, 0.02),
        ev.f,
    )

    # Discrete total-current identity
    j_total = (
        ev.fluxes.jn_edges
        + ev.fluxes.jp_edges
    )

    f = ev.f.reshape(m, 3)

    np.testing.assert_allclose(
        f[1:-1, 1] + f[1:-1, 2],
        np.diff(j_total)
        / device.mesh.w_interior_nd,
        rtol=1e-12,
        atol=1e-10,
    )


def test_qf_reference_shift_invariance():
    device = make_device()
    m = device.mesh.node_count

    u = np.linspace(-1.5, 1.6, m)
    en = np.linspace(-0.3, 0.2, m)
    ep = np.linspace(0.1, -0.4, m)

    original = QFState(
        u=u,
        eta_n=en,
        eta_p=ep,
    )

    # Change references while keeping actual
    # quasi-Fermi levels unchanged.
    shifted = QFState(
        u=u,
        eta_n=en - 0.7,
        eta_p=ep + 0.4,
        ref_n=0.7,
        ref_p=-0.4,
    )

    ev0 = evaluate_qf(
        original,
        device,
        0.0,
    )

    ev1 = evaluate_qf(
        shifted,
        device,
        0.0,
    )

    np.testing.assert_allclose(
        ev1.f,
        ev0.f,
        rtol=1e-12,
        atol=1e-12,
    )

    np.testing.assert_allclose(
        ev1.fluxes.jn_edges,
        ev0.fluxes.jn_edges,
        rtol=1e-12,
        atol=1e-12,
    )

    np.testing.assert_allclose(
        ev1.fluxes.jp_edges,
        ev0.fluxes.jp_edges,
        rtol=1e-12,
        atol=1e-12,
    )

    np.testing.assert_allclose(
        ev1.recombination_nd,
        ev0.recombination_nd,
        rtol=1e-12,
        atol=1e-12,
    )


def test_invalid_qf_inputs():
    with pytest.raises(ValueError):
        QFState(
            u=np.zeros(2),
            eta_n=np.zeros(2),
            eta_p=np.zeros(2),
        )

    with pytest.raises(ValueError):
        QFState(
            u=np.zeros(3),
            eta_n=np.zeros(2),
            eta_p=np.zeros(3),
        )

    with pytest.raises(ValueError):
        QFState(
            u=np.zeros(3),
            eta_n=np.array(
                [0.0, float("nan"), 0.0]
            ),
            eta_p=np.zeros(3),
        )

    with pytest.raises(ValueError):
        QFState.unpack_qf(
            np.zeros(8),
            3,
        )

    with pytest.raises(
        (ValueError, FloatingPointError)
    ):
        QFState(
            u=np.full(3, 800.0),
            eta_n=np.zeros(3),
            eta_p=np.zeros(3),
        ).to_physical()

    device = make_device()

    small_qf = QFState(
        u=np.zeros(3),
        eta_n=np.zeros(3),
        eta_p=np.zeros(3),
    )

    with pytest.raises(
        ValueError,
        match="node_count",
    ):
        evaluate_qf(
            small_qf,
            device,
            0.0,
        )