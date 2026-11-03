"""W2-B3：``NormalizedDevice`` 的 cm 基准尺度（附录 C）与 D_0 二选一（第 3.5 节）。"""

from __future__ import annotations

import numpy as np
import pytest

from minitcad.core.device import (
    DeviceConfig,
    MeshSpec,
    NormalizationConfig,
    build_device,
)
from minitcad.core.material import SILICON, MaterialTemperatureMismatch
from minitcad.core.region import Region

LENGTH_CM = 1.0e-3
JUNCTION_CM = 2.5e-4


def _config(**overrides) -> DeviceConfig:
    params = dict(
        length_cm=LENGTH_CM,
        area_cm2=1.0e-4,
        temperature_k=300.0,
        material=SILICON,
        regions=(
            Region(0.0, JUNCTION_CM, "silicon", acceptor_cm3=1.0e16, donor_cm3=0.0),
            Region(JUNCTION_CM, LENGTH_CM, "silicon", acceptor_cm3=0.0, donor_cm3=2.0e15),
        ),
        mesh=MeshSpec(kind="uniform", nodes=21),
        normalization=NormalizationConfig(reference_diffusivity_cm2_s=10.0),
    )
    params.update(overrides)
    return DeviceConfig(**params)


def test_scales_match_appendix_c():
    device = build_device(_config())
    assert device.thermal_voltage_v == pytest.approx(0.0258519998, rel=1e-8)
    assert device.l0_cm == pytest.approx(4.08845452e-3, rel=1e-6)
    assert device.scales.t0_s == pytest.approx(1.67154603e-6, rel=1e-5)
    assert device.scales.j0_a_cm2 == pytest.approx(3.91878307e-6, rel=1e-5)
    assert device.ni_cm3 == 1.0e10


def test_reference_mobility_variant_gives_same_length_and_different_d0():
    by_d0 = build_device(_config())
    # mu_0 = 1000 cm^2/(V s) -> D_0 = mu_0 * V_T（第 3.5 节的备选参数化）
    by_mu0 = build_device(
        _config(normalization=NormalizationConfig(reference_mobility_cm2_per_Vs=1000.0))
    )
    assert by_mu0.l0_cm == pytest.approx(by_d0.l0_cm, rel=1e-12)
    assert by_mu0.d0_cm2_s == pytest.approx(1000.0 * by_d0.thermal_voltage_v, rel=1e-12)
    assert by_mu0.d0_cm2_s != pytest.approx(by_d0.d0_cm2_s)


def test_normalization_requires_exactly_one_reference():
    with pytest.raises(ValueError):
        NormalizationConfig()
    with pytest.raises(ValueError):
        NormalizationConfig(
            reference_diffusivity_cm2_s=10.0, reference_mobility_cm2_per_Vs=1000.0
        )


def test_normalized_arrays_have_expected_lengths():
    device = build_device(_config())
    m = device.mesh.node_count
    assert device.net_doping_nd.size == m
    assert device.cbar_nd.size == m - 2
    assert device.mesh.normalization_scale() == pytest.approx(device.l0_cm, rel=1e-12)


def test_dimensionless_mobilities_and_lifetimes_use_scales():
    device = build_device(_config())
    assert device.mu_n_nd == pytest.approx(1400.0 / device.scales.mu0_cm2_per_vs)
    assert device.tau_n_nd == pytest.approx(1.0e-6 / device.scales.t0_s)


def test_temperature_mismatch_blocks_build():
    with pytest.raises(MaterialTemperatureMismatch):
        build_device(_config(temperature_k=350.0))


def test_segmented_mesh_spec_builds():
    device = build_device(
        _config(
            mesh=MeshSpec(
                kind="segmented",
                breakpoints_cm=(0.0, JUNCTION_CM, LENGTH_CM),
                per_segment_nodes=(11, 31),
            )
        )
    )
    assert device.mesh.node_count == 41  # 11 + 31 - 1 共享结点


def test_segmented_mesh_must_include_junction_breakpoint():
    # 声明 segmented 就必须把结区坐标写进断点，否则“结区加密”没有兑现（第 5.1、5.3 节）。
    with pytest.raises(ValueError) as excinfo:
        build_device(
            _config(
                mesh=MeshSpec(
                    kind="segmented",
                    breakpoints_cm=(0.0, 5.0e-4, LENGTH_CM),
                    per_segment_nodes=(6, 11),
                )
            )
        )
    assert "junction" in str(excinfo.value)


def test_segmented_mesh_puts_the_junction_on_a_node():
    device = build_device(
        _config(
            mesh=MeshSpec(
                kind="segmented",
                breakpoints_cm=(0.0, JUNCTION_CM, LENGTH_CM),
                per_segment_nodes=(11, 31),
            )
        )
    )
    assert np.isclose(device.mesh.x_nodes_cm, JUNCTION_CM, rtol=1e-12, atol=0.0).any()


def test_uniform_mesh_keeps_junction_inside_a_control_volume():
    # 均匀网格的结区可以落在控制体内部，由 cbar_nd 的控制体积分吸收（第 5.1 节），
    # 因此这项约束只对 segmented 生效。
    device = build_device(_config(mesh=MeshSpec(kind="uniform", nodes=7)))
    assert not np.isclose(device.mesh.x_nodes_cm, JUNCTION_CM, rtol=1e-12, atol=0.0).any()
    assert np.all(np.isfinite(device.cbar_nd))


def test_right_terminal_must_be_grounded():
    with pytest.raises(ValueError):
        _config(right_voltage_V=0.1)


def test_mesh_spec_requires_integer_node_counts():
    with pytest.raises(ValueError):
        MeshSpec(kind="uniform", nodes=21.5)
    with pytest.raises(ValueError):
        MeshSpec(kind="segmented", breakpoints_cm=(0.0, LENGTH_CM), per_segment_nodes=(2.5,))
