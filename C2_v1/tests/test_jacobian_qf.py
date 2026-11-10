
"""Week 5 Day 3: QF Jacobian smoke tests + directional FD check."""

import numpy as np
import pytest
from scipy.sparse import issparse
from scipy.sparse.linalg import spsolve

from minitcad.core.device import (
    DeviceConfig, MeshSpec, NormalizationConfig, build_device,
)
from minitcad.core.material import SILICON
from minitcad.core.region import Region
from minitcad.core.state import QFState
from minitcad.discretization.residual import residual_qf
from minitcad.solvers.jacobian import jacobian_qf


def make_case(nodes=31):
    device = build_device(DeviceConfig(
        length_cm=1e-3,
        area_cm2=1e-4,
        temperature_k=300.0,
        material=SILICON,
        regions=(
            Region(0.0, 5e-4, "silicon", 1e11, 0.0),
            Region(5e-4, 1e-3, "silicon", 0.0, 2e11),
        ),
        mesh=MeshSpec(kind="uniform", nodes=nodes),
        normalization=NormalizationConfig(
            reference_diffusivity_cm2_s=10.0,
        ),
    ))
    bias = 0.0
    target = device.contact_targets(bias)
    x = np.linspace(0.0, 1.0, nodes)

    qf = QFState(
        u=np.linspace(target[0], target[3], nodes)
          + 0.02 * np.sin(np.pi * x),
        eta_n=0.01 * np.sin(2.0 * np.pi * x) - 0.3,
        eta_p=0.01 * np.cos(2.0 * np.pi * x) + 0.2,
        ref_n=0.3,
        ref_p=-0.2,
    )
    return device, bias, qf


def test_sparse_matrix_shape_and_locality():
    device, bias, qf = make_case(nodes=7)
    m = qf.node_count
    jac = jacobian_qf(qf, device, bias)

    assert issparse(jac)
    assert jac.shape == (3*m, 3*m)
    assert jac.nnz <= 27*m
    assert np.isfinite(jac.data).all()

    rr, cc = jac.nonzero()
    assert np.all(np.abs(rr // 3 - cc // 3) <= 1)


@pytest.mark.parametrize("nodes", [7, 31])
def test_boundary_qf_dirichlet_derivatives(nodes):
    device, bias, qf = make_case(nodes)
    jac = jacobian_qf(qf, device, bias).toarray()

    for row in (
        0, 1, 2,
        3*nodes-3, 3*nodes-2, 3*nodes-1
    ):
        expected = np.zeros(3*nodes)
        expected[row] = 1.0
        np.testing.assert_allclose(
            jac[row], expected, atol=1e-8, rtol=0.0
        )


def test_directional_derivative_T10_initial():
    device, bias, qf = make_case(nodes=31)
    z = qf.pack_qf()
    jac = jacobian_qf(qf, device, bias)

    rng = np.random.default_rng(2026)
    direction = rng.standard_normal(z.size)
    direction /= np.max(np.abs(direction))
    eps = 1e-5

    def f_from_z(vec):
        trial = QFState.unpack_qf(
            vec,
            qf.node_count,
            ref_n=qf.ref_n,
            ref_p=qf.ref_p,
        )
        return residual_qf(trial, device, bias)

    independent_fd = (
        f_from_z(z + eps * direction)
        - f_from_z(z - eps * direction)
    ) / (2.0 * eps)

    actual = jac @ direction

    error = np.max(
        np.abs(actual - independent_fd)
    ) / max(
        1.0, np.max(np.abs(independent_fd))
    )

    print("T10 preliminary directional error:", error)
    assert error <= 1e-4


def test_one_linear_newton_step_is_finite():
    device, bias, qf = make_case(nodes=31)
    jac = jacobian_qf(qf, device, bias)
    f = residual_qf(qf, device, bias)

    dz = spsolve(jac, -f)

    assert dz.shape == (3 * qf.node_count,)
    assert np.isfinite(dz).all()


def test_reference_shift_same_jacobian():
    device, bias, qf = make_case(nodes=7)

    shifted = QFState(
        u=qf.u,
        eta_n=qf.eta_n - 0.4,
        eta_p=qf.eta_p + 0.6,
        ref_n=qf.ref_n + 0.4,
        ref_p=qf.ref_p - 0.6,
    )

    jac_a = jacobian_qf(qf, device, bias).toarray()
    jac_b = jacobian_qf(shifted, device, bias).toarray()

    np.testing.assert_allclose(
        jac_a, jac_b, rtol=2e-4, atol=1e-6
    )


def test_poisson_derivative_signs():
    device, bias, qf = make_case(nodes=7)
    jac = jacobian_qf(qf, device, bias)
    physical = qf.to_physical()

    h = device.mesh.h_edges_nd
    w = device.mesh.w_interior_nd
    i = 3
    row = 3 * i

    expected = -(1 / h[i - 1] + 1 / h[i]) / w[i - 1]
    expected -= physical.n[i] + physical.p[i]

    np.testing.assert_allclose(
        jac[row, row], expected, rtol=1e-7
    )
    np.testing.assert_allclose(
        jac[row, row + 1], -physical.n[i], rtol=1e-7
    )
    np.testing.assert_allclose(
        jac[row, row + 2], -physical.p[i], rtol=1e-7
    )


def test_standard_201_nodes_smoke():
    device, bias, qf = make_case(nodes=201)
    jac = jacobian_qf(qf, device, bias)

    assert jac.shape == (603, 603)
    assert jac.nnz < 0.05 * 603**2
    assert np.isfinite(jac.data).all()