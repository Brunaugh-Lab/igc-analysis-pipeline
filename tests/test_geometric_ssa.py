"""Tests for cumulative-Q3 to D32/geometric-SSA workflow."""

from pathlib import Path

import numpy as np
import pytest

from igc_sea.analysis.geometric_ssa import calculate_d32
from igc_sea.io.particle_size import VolumeDistribution, read_cumulative_q3


def _distribution(q=(0.0, 50.0, 100.0)):
    return VolumeDistribution(
        source_file=Path("sample.csv"),
        metadata={"formulation_id": "F1", "Dispersing system": "dry"},
        boundaries_um=np.array([1.0, 2.0, 4.0]),
        cumulative_q3_percent=np.array(q),
    )


def test_d32_uses_q3_increments_at_geometric_bin_midpoints():
    result = calculate_d32(_distribution())
    midpoints = np.sqrt([1.0 * 2.0, 2.0 * 4.0])
    expected = 1.0 / np.sum(np.array([0.5, 0.5]) / midpoints)
    assert result.d32_um == pytest.approx(expected)
    assert result.d32_bin_lower_um == pytest.approx(1.0 / (0.5 / 1.0 + 0.5 / 2.0))
    assert result.d32_bin_upper_um == pytest.approx(1.0 / (0.5 / 2.0 + 0.5 / 4.0))


def test_density_adds_geometric_ssa_and_conservative_dose_value():
    result = calculate_d32(_distribution(), density_g_cm3=1.25, density_basis="skeletal")
    assert result.ssa_geo_m2_g == pytest.approx(6 / (1.25 * result.d32_um))
    assert result.ssa_dose_m2_g == pytest.approx(6 / (1.25 * result.d32_bin_upper_um))
    assert result.density_basis == "skeletal"


def test_nonmonotonic_q3_is_rejected():
    with pytest.raises(ValueError, match="monotonic"):
        calculate_d32(_distribution(q=(0.0, 60.0, 50.0)))


def test_first_populated_bin_qc():
    result = calculate_d32(_distribution(), first_bin_warning_percent=40.0)
    assert "LOW_END_SENSITIVE" in result.qc_flags


def test_reader_selects_uppercase_cumulative_q3(tmp_path):
    path = tmp_path / "paq.csv"
    path.write_text(
        "Product,Dispersing system,formulation_id\n"
        "Sample,dry,F1\n"
        "xo / µm,Q₃ / %,xm / µm,q₃ / mm⁻¹,xm / µm,q₃ lg\n"
        "0.5,0,0.7,999,0.67,9\n"
        "0.9,40,1.0,888,0.99,8\n"
        "1.1,100,1.2,777,1.19,7\n",
        encoding="utf-8",
    )
    distribution = read_cumulative_q3(path)
    assert distribution.metadata["formulation_id"] == "F1"
    assert distribution.cumulative_q3_percent.tolist() == [0.0, 40.0, 100.0]
