
"""Week 6/7 A-side numerical regression tests."""

import numpy as np
import pytest

from minitcad.core.device import (
    DeviceConfig,
    MeshSpec,
    NormalizationConfig,
    build_device,
)
from minitcad.core.material import SILICON
from minitcad.core.region import Region
from minitcad.core.state import QFState
from minitcad.discretization.residual import (
    evaluate_qf,
    residual_qf,
)
from minitcad.simulation.equilibrium import (
    solve_equilibrium_scalar,
)
from minitcad.solvers import solve_gummel, solve_newton
from minitcad.solvers.jacobian import jacobian_qf


def make_case(doping=1e14, nodes=101, segmented=False):
    L = 2e-4

    if segmented:
        mesh = MeshSpec(
            kind="segmented",
            breakpoints_cm=(0., 8e-5, 1e-4, 1.2e-4, L),
            per_segment_nodes=(21, 41, 41, 21),
        )
    else:
        mesh = MeshSpec(kind="uniform", nodes=nodes)

    device = build_device(DeviceConfig(
        length_cm=L,
        area_cm2=1e-4,
        temperature_k=300.,
        material=SILICON,
        regions=(
            Region(0., L / 2, "silicon", doping, 0.),
            Region(L / 2, L, "silicon", 0., doping),
        ),
        mesh=mesh,
        normalization=NormalizationConfig(
            reference_diffusivity_cm2_s=10.,
        ),
    ))

    targets = device.contact_targets(0.)

    eq = solve_equilibrium_scalar(
        device.mesh.x_nodes_nd,
        device.cbar_nd,
        targets[0],
        targets[3],
    )

    assert eq.converged, eq.reason

    return device, eq.state


@pytest.mark.parametrize("doping", [1e12, 1e14, 1e16])
@pytest.mark.parametrize("bias", [0.05, 0.20, -0.10])
def test_gummel_pressure(doping, bias):
    device, start = make_case(doping=doping)

    result = solve_gummel(device, bias, start)

    assert result.converged, result.reason
    assert result.physics is not None
    assert result.residual_norm <= 1e-8
    assert result.update_norm <= 1e-6

    flux = result.physics.fluxes
    j = flux.jn_edges + flux.jp_edges

    nonconservation = (
        np.max(np.abs(j - np.median(j)))
        / max(1., abs(np.median(j)))
    )

    assert nonconservation <= 1e-3


def test_newton_continuation():
    device, state = make_case()

    for bias in (0.05, 0.10, 0.15, 0.20):
        result = solve_newton(device, bias, state)

        assert result.converged, (
            bias, result.reason
        )
        assert result.physics is not None

        # Only accepted results become initial states.
        state = result.state


def test_jacobian_directional_fd():
    device, physical = make_case(nodes=31)
    qf = QFState.from_physical(physical)

    z = qf.pack_qf()
    direction = np.random.default_rng(2026).normal(
        size=z.size
    )
    direction /= np.max(np.abs(direction))

    eps = 1e-5

    def evaluate_z(vector):
        candidate = QFState.unpack_qf(
            vector,
            qf.node_count,
            ref_n=qf.ref_n,
            ref_p=qf.ref_p,
        )
        return residual_qf(candidate, device, 0.)

    fd = (
        evaluate_z(z + eps * direction)
        - evaluate_z(z - eps * direction)
    ) / (2.0 * eps)

    jd = jacobian_qf(qf, device, 0.) @ direction

    error = (
        np.max(np.abs(fd - jd))
        / max(1., np.max(np.abs(fd)))
    )

    assert error <= 1e-4


def test_nonuniform_qf_diagnostics():
    device, start = make_case(segmented=True)

    # Existing internal Gummel test hook:
    debug = {}

    result = solve_gummel(
        device,
        0.05,
        start,
        diagnostics=debug,
    )

    assert result.converged
    assert result.physics is not None

    # Construct from exact accepted QF data.
    # Do not reconstruct QF from physical State here.
    qf = QFState(
        debug["u"],
        debug["eta_n"],
        debug["eta_p"],
        debug["ref_n"],
        debug["ref_p"],
    )

    ev = evaluate_qf(qf, device, 0.05)

    np.testing.assert_allclose(
        result.physics.residual_qf_nd,
        ev.f,
        rtol=1e-8,
        atol=1e-8,
    )
    np.testing.assert_allclose(
        result.physics.fluxes.jn_edges,
        ev.fluxes.jn_edges,
        rtol=1e-8,
        atol=1e-8,
    )
    np.testing.assert_allclose(
        result.physics.fluxes.jp_edges,
        ev.fluxes.jp_edges,
        rtol=1e-8,
        atol=1e-8,
    )


def test_mesh_current_convergence():
    currents = []

    for nodes in (101, 201, 401):
        device, start = make_case(nodes=nodes)
        result = solve_gummel(device, 0.05, start)

        assert result.converged

        flux = result.physics.fluxes
        current_nd = (
            flux.jn_edges[0] + flux.jp_edges[0]
        )

        currents.append(
            current_nd * device.scales.j0_a_cm2
        )

    scale = max(
        abs(currents[-1]),
        abs(currents[-2]),
    )

    # This benchmark has demonstrably nonzero J.
    assert scale > 0.

    error = (
        abs(currents[-1] - currents[-2])
        / scale
    )

    assert error <= 0.01