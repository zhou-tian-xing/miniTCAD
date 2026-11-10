
"""Reverse-bias numerical stability tests."""

import numpy as np
import pytest

from test_stability import make_case
from minitcad.simulation.reverse_sweep import (
    solve_reverse_sweep,
)


@pytest.mark.parametrize(
    "doping",
    [1e14, 1e16],
)
def test_uniform_reverse_sweep_reaches_minus_one_volt(
    doping,
):
    device, equilibrium = make_case(
        doping=doping,
    )

    result = solve_reverse_sweep(
        device,
        equilibrium,
        -1.0,
    )

    assert result.reached_target
    assert result.points[-1].bias_V == pytest.approx(
        -1.0
    )
    assert len(result.points) > 2

    biases = [
        row.bias_V for row in result.points
    ]
    assert np.all(np.diff(biases) < 0)

    for row in result.points:
        out = row.result

        assert out.converged
        assert out.state is not None
        assert out.physics is not None
        assert out.residual_norm <= 1e-8
        assert out.update_norm <= 1e-6


def test_segmented_reverse_bias_with_explicit_tolerance():
    device, equilibrium = make_case(
        doping=1e16,
        segmented=True,
    )

    result = solve_reverse_sweep(
        device,
        equilibrium,
        -0.4,
        initial_step_V=0.03,
        max_step_V=0.07,
        tol_residual=1e-7,
    )

    assert result.reached_target

    assert result.points[-1].bias_V == pytest.approx(
        -0.4
    )

    for item in result.points:
        assert item.result.residual_norm <= 1e-7


def test_invalid_reverse_target_is_rejected():
    device, equilibrium = make_case()

    with pytest.raises(ValueError):
        solve_reverse_sweep(
            device,
            equilibrium,
            0.2,
        )