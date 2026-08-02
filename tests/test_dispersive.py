"""Tests for dispersive surface energy analysis (Dorris-Gray method)."""

import pytest
import numpy as np
import pandas as pd

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from igc_sea.constants import R_GAS, N_AVOGADRO, A_CH2, gamma_ch2
from igc_sea.analysis.dispersive import dorris_gray_gamma_d, calculate_net_retention_volume
from igc_sea.analysis.retention import james_martin_correction, net_retention_volume, specific_retention_volume
from igc_sea.analysis.peak_detection import (
    detect_baseline, subtract_baseline, find_peak_max, find_peak_cofm,
    integrate_peak, process_chromatogram,
    asymmetry_factor, tailing_factor, _half_widths_at_fraction,
)
from igc_sea.analysis.interpolation import interpolate_to_coverage
from igc_sea.analysis.calibration import (
    moles_from_area, monolayer_capacity, actual_coverage,
)


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

class TestConstants:
    def test_gamma_ch2_at_20C(self):
        """γ_CH₂ at 293.15 K should be 35.6 mJ/m²."""
        assert gamma_ch2(293.15) == pytest.approx(35.6e-3, rel=1e-6)

    def test_gamma_ch2_at_30C(self):
        """γ_CH₂ at 303.15 K should be 35.02 mJ/m²."""
        assert gamma_ch2(303.15) == pytest.approx(35.02e-3, rel=1e-3)

    def test_gamma_ch2_decreases_with_temperature(self):
        assert gamma_ch2(310) < gamma_ch2(293.15)


class TestDorrisGrayInputValidation:
    @pytest.mark.parametrize("bad_vn", [0.0, -1.0, np.nan, np.inf])
    def test_nonpositive_or_nonfinite_vn_raises(self, bad_vn):
        df = pd.DataFrame({
            "coverage": [0.01, 0.01],
            "carbon_number": [8, 9],
            "VN": [10.0, bad_vn],
        })
        with pytest.raises(ValueError, match="strictly positive"):
            dorris_gray_gamma_d(df, 303.15)


# ---------------------------------------------------------------------------
# Peak detection
# ---------------------------------------------------------------------------

class TestPeakDetection:
    @pytest.fixture
    def synthetic_peak(self):
        """Generate a synthetic Gaussian peak on a linear baseline."""
        time = np.linspace(0, 5, 3000)
        # Baseline: y = 100 + 2*t
        baseline = 100 + 2 * time
        # Gaussian peak centered at t=1.5, sigma=0.1, height=5000
        peak = 5000 * np.exp(-0.5 * ((time - 1.5) / 0.1) ** 2)
        signal = baseline + peak
        return time, signal, 100.0, 2.0, 1.5

    def test_baseline_detection(self, synthetic_peak):
        time, signal, true_intercept, true_gradient, _ = synthetic_peak
        intercept, gradient = detect_baseline(time, signal)
        # Peak-aware baseline uses narrower regions flanking the peak,
        # so allow wider tolerance than a simple edge-fit would need.
        assert intercept == pytest.approx(true_intercept, rel=0.1)
        assert gradient == pytest.approx(true_gradient, rel=0.5)

    def test_peak_max_time(self, synthetic_peak):
        time, signal, intercept, gradient, true_peak_time = synthetic_peak
        corrected = subtract_baseline(time, signal, intercept, gradient)
        t_max = find_peak_max(time, corrected)
        assert t_max == pytest.approx(true_peak_time, abs=0.01)

    def test_peak_cofm(self, synthetic_peak):
        time, signal, intercept, gradient, true_peak_time = synthetic_peak
        corrected = subtract_baseline(time, signal, intercept, gradient)
        t_cofm = find_peak_cofm(time, corrected)
        # CoM should be very close to peak max for symmetric Gaussian
        assert t_cofm == pytest.approx(true_peak_time, abs=0.01)

    def test_peak_area_positive(self, synthetic_peak):
        time, signal, intercept, gradient, _ = synthetic_peak
        corrected = subtract_baseline(time, signal, intercept, gradient)
        area = integrate_peak(time, corrected)
        assert area > 0

    def test_process_chromatogram(self, synthetic_peak):
        time, signal, _, _, _ = synthetic_peak
        result = process_chromatogram(time, signal)
        assert "peak_max_time" in result
        assert "peak_cofm" in result
        assert "peak_area" in result
        assert result["peak_area"] > 0
        # New asymmetry keys
        assert "asymmetry_factor" in result
        assert "tailing_factor" in result
        assert "com_max_divergence_min" in result
        assert "com_max_divergence_frac" in result


# ---------------------------------------------------------------------------
# Peak asymmetry
# ---------------------------------------------------------------------------

class TestPeakAsymmetry:
    """Test peak asymmetry and tailing factor computations."""

    @pytest.fixture
    def symmetric_peak(self):
        """Symmetric Gaussian peak: A_s should be ~1.0."""
        time = np.linspace(0, 5, 5000)
        signal = 5000 * np.exp(-0.5 * ((time - 2.5) / 0.1) ** 2)
        return time, signal

    @pytest.fixture
    def tailing_peak(self):
        """Asymmetric peak with exponential tail: A_s > 1.0."""
        time = np.linspace(0, 5, 5000)
        # Exponentially modified Gaussian (EMG) — sharp rise, slow decay
        t0 = 2.0
        sigma = 0.05
        tau = 0.15  # tail time constant
        from scipy.special import erfc
        z = (1 / tau) * (sigma**2 / (2 * tau) - (time - t0))
        # Clamp z to prevent overflow
        z = np.clip(z, -50, 50)
        signal = (5000 * sigma / tau * np.sqrt(np.pi / 2)
                  * np.exp(0.5 * (sigma / tau)**2 - (time - t0) / tau)
                  * erfc(z / np.sqrt(2)))
        signal = np.maximum(signal, 0.0)
        return time, signal

    def test_symmetric_asymmetry_near_one(self, symmetric_peak):
        time, signal = symmetric_peak
        af = asymmetry_factor(time, signal)
        assert af == pytest.approx(1.0, abs=0.05)

    def test_symmetric_tailing_near_one(self, symmetric_peak):
        time, signal = symmetric_peak
        tf = tailing_factor(time, signal)
        assert tf == pytest.approx(1.0, abs=0.05)

    def test_tailing_asymmetry_greater_than_one(self, tailing_peak):
        time, signal = tailing_peak
        af = asymmetry_factor(time, signal)
        assert af > 1.3  # EMG with tau=0.15 should tail significantly

    def test_tailing_factor_greater_than_one(self, tailing_peak):
        time, signal = tailing_peak
        tf = tailing_factor(time, signal)
        assert tf > 1.1

    def test_half_widths_symmetric(self, symmetric_peak):
        time, signal = symmetric_peak
        a, b = _half_widths_at_fraction(time, signal, 0.5)
        # For symmetric Gaussian, leading ≈ trailing
        assert a == pytest.approx(b, rel=0.05)

    def test_half_widths_tailing(self, tailing_peak):
        time, signal = tailing_peak
        a, b = _half_widths_at_fraction(time, signal, 0.10)
        # Trailing should be wider than leading for tailing peak
        assert b > a

    def test_zero_signal(self):
        time = np.linspace(0, 5, 100)
        signal = np.zeros(100)
        af = asymmetry_factor(time, signal)
        assert np.isnan(af)


# ---------------------------------------------------------------------------
# Retention volume
# ---------------------------------------------------------------------------

class TestRetention:
    def test_james_martin_equal_pressure(self):
        """When P_in = P_out, J should be 1.0."""
        assert james_martin_correction(1.0, 1.0) == pytest.approx(1.0)

    def test_james_martin_typical(self):
        """Typical column with 2:1 pressure ratio."""
        j = james_martin_correction(2.0, 1.0)
        # J = 3/2 * (4-1)/(8-1) = 3/2 * 3/7 = 9/14 ≈ 0.6429
        assert j == pytest.approx(9 / 14, rel=1e-6)

    def test_net_retention_volume_basic(self):
        """V_N = (t_R - t_0) * F."""
        v = net_retention_volume(t_R=2.0, t_0=0.5, flow_rate=10.0)
        assert v == pytest.approx(15.0)

    def test_specific_retention_volume(self):
        v_spec = specific_retention_volume(v_n=15.0, sample_mass_g=0.1)
        assert v_spec == pytest.approx(150.0)


# ---------------------------------------------------------------------------
# Interpolation
# ---------------------------------------------------------------------------

class TestInterpolation:
    def test_linear_interpolation(self):
        df = pd.DataFrame({
            "solvent_name": ["OCTANE"] * 3,
            "actual_coverage": [0.01, 0.05, 0.10],
            "net_retention_volume": [100.0, 80.0, 60.0],
        })
        result = interpolate_to_coverage(df, [0.03, 0.075])
        assert len(result) == 2
        # At 0.03: linear interp between (0.01, 100) and (0.05, 80)
        # = 100 + (0.03-0.01)/(0.05-0.01) * (80-100) = 100 - 10 = 90
        assert result.iloc[0]["interpolated_VN"] == pytest.approx(90.0)

    def test_extrapolation_returns_nan(self):
        df = pd.DataFrame({
            "solvent_name": ["OCTANE"] * 2,
            "actual_coverage": [0.05, 0.10],
            "net_retention_volume": [80.0, 60.0],
        })
        result = interpolate_to_coverage(df, [0.01, 0.20])
        assert np.isnan(result.iloc[0]["interpolated_VN"])
        assert np.isnan(result.iloc[1]["interpolated_VN"])


# ---------------------------------------------------------------------------
# Dorris-Gray γ_d
# ---------------------------------------------------------------------------

class TestDorrisGray:
    @pytest.fixture
    def alkane_data(self):
        """Synthetic alkane data at a single coverage."""
        # Use realistic V_N values that should give ~35 mJ/m²
        T = 303.15
        return pd.DataFrame({
            "coverage": [0.005, 0.005, 0.005],
            "carbon_number": [8, 9, 10],
            "VN": [17.0, 50.0, 132.0],  # mL/g — realistic
        }), T

    def test_returns_dataframe(self, alkane_data):
        df, T = alkane_data
        result = dorris_gray_gamma_d(df, T)
        assert isinstance(result, pd.DataFrame)
        assert "gamma_d_mJm2" in result.columns
        assert "r_squared" in result.columns

    def test_positive_gamma_d(self, alkane_data):
        df, T = alkane_data
        result = dorris_gray_gamma_d(df, T)
        assert result["gamma_d_mJm2"].iloc[0] > 0

    def test_r_squared_near_unity(self, alkane_data):
        df, T = alkane_data
        result = dorris_gray_gamma_d(df, T)
        # Alkane series should give very high R²
        assert result["r_squared"].iloc[0] > 0.95

    def test_empty_alkanes_raises(self):
        df = pd.DataFrame({
            "coverage": [0.005],
            "carbon_number": [np.nan],
            "VN": [17.0],
        })
        with pytest.raises(ValueError, match="No alkane data"):
            dorris_gray_gamma_d(df, 303.15)


# ---------------------------------------------------------------------------
# Calibration → actual surface coverage
# ---------------------------------------------------------------------------

class TestCalibration:
    """Test source-neutral power-law calibration calculations."""

    def test_moles_from_area_positive(self):
        """Moles should be positive for positive area."""
        C1, C2 = 5.0e-11, 0.98
        moles = moles_from_area(6875.23, C1, C2)
        assert moles > 0

    def test_moles_from_area_zero(self):
        """Zero area should give zero moles."""
        assert moles_from_area(0.0, 1e-10, 1.0) == 0.0

    def test_monolayer_capacity(self):
        """Monolayer capacity with known BET/mass/cross-section."""
        n_mono = monolayer_capacity(10.0, 0.090, 6.3e-19)
        # BET=10 m²/g, mass=0.090 g, octane cross-section=6.3e-19 m²
        expected = 10.0 * 0.090 / (6.02214076e23 * 6.3e-19)
        assert n_mono == pytest.approx(expected, rel=1e-6)

    def test_actual_coverage_matches_declared_equation(self):
        C1 = 5.0e-11
        C2 = 0.98
        peak_area = 7000.0
        bet_ssa = 10.0
        mass_g = 0.090
        a_cross = 6.3e-19

        cov = actual_coverage(peak_area, C1, C2, bet_ssa, mass_g, a_cross)
        expected_moles = C1 * peak_area**C2
        expected_capacity = bet_ssa * mass_g / (N_AVOGADRO * a_cross)
        assert cov == pytest.approx(expected_moles / expected_capacity)
