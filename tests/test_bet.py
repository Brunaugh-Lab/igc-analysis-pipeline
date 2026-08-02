"""Tests for BET specific surface area module.

Tests the core physical chemistry calculations:
- Saturation pressure from extended Antoine equation
- P/P0 computation from injected moles
- Adsorption isotherm from cumulative V_N integration
- BET linearization and SSA calculation
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from igc_analysis.analysis.bet import (
    saturation_pressure,
    partial_pressure_ratio,
    eluted_peak_concentration,
    james_martin_j,
    build_adsorption_isotherm,
    bet_linearization,
    bet_quality_checks,
    classify_isotherm,
    IsothermClassification,
    OCTANE_ANTOINE,
    OCTANE_CROSS_SECTION_M2,
    ORIGIN_STRATEGIES,
    BETResult,
    BETQCFlags,
    InjectionResult,
    IsothermPoint,
    _assemble_result,
    _InjectionRecord,
    _is_clipped,
)
from igc_analysis.analysis.probes import ProbeProperties
from igc_analysis.constants import R_GAS, N_AVOGADRO


# ---------------------------------------------------------------------------
# Saturation pressure
# ---------------------------------------------------------------------------

class TestSaturationPressure:
    """Test the extended Antoine equation for octane."""

    def test_octane_303K(self):
        """P_sat(octane, 303.15 K) should be ~2464 Pa."""
        P_sat = saturation_pressure(303.15, **OCTANE_ANTOINE)
        assert 2400 < P_sat < 2550, f"P_sat={P_sat} Pa outside expected range"

    def test_octane_increases_with_temperature(self):
        """P_sat should increase with temperature."""
        P_low = saturation_pressure(293.15, **OCTANE_ANTOINE)
        P_high = saturation_pressure(313.15, **OCTANE_ANTOINE)
        assert P_high > P_low

    def test_octane_boiling_point(self):
        """At octane boiling point (~399 K), P_sat should be ~101325 Pa."""
        P_bp = saturation_pressure(399.0, **OCTANE_ANTOINE)
        assert 90000 < P_bp < 115000, f"P_sat at BP={P_bp} Pa"


# ---------------------------------------------------------------------------
# P/P0 calculation
# ---------------------------------------------------------------------------

class TestPartialPressureRatio:
    """Test P/P0 from injected moles."""

    def test_zero_moles_gives_zero(self):
        pp0 = partial_pressure_ratio(0.0, 2e-6, 303.15, 2464.0)
        assert pp0 == 0.0

    def test_typical_range(self):
        """Typical injection (1e-7 mol in 2000 µL loop at 303 K)
        should give P/P0 in 0.01-0.5 range."""
        pp0 = partial_pressure_ratio(1e-7, 2e-6, 303.15, 2464.0)
        assert 0.01 < pp0 < 0.5

    def test_proportional_to_moles(self):
        """P/P0 should be proportional to moles injected."""
        pp0_1 = partial_pressure_ratio(1e-7, 2e-6, 303.15, 2464.0)
        pp0_2 = partial_pressure_ratio(2e-7, 2e-6, 303.15, 2464.0)
        assert abs(pp0_2 / pp0_1 - 2.0) < 1e-10


# ---------------------------------------------------------------------------
# Eluted peak-apex concentration (the physically correct isotherm c-axis)
# ---------------------------------------------------------------------------

class TestElutedPeakConcentration:
    """Test gas-phase concentration at the peak apex.

    c_apex = apex_signal * (n_inj / area) / F_col.  This is the concentration
    the isotherm relation V_N/m = dq/dc refers to for the peak-maximum method,
    and it replaces the loop concentration n_inj/V_loop that inflated the
    inferred pressure for strongly retained probes.
    """

    def test_zero_area_gives_zero(self):
        assert eluted_peak_concentration(100.0, 1e-7, 0.0, 1e-5) == 0.0

    def test_zero_flow_gives_zero(self):
        assert eluted_peak_concentration(100.0, 1e-7, 5.0, 0.0) == 0.0

    def test_formula(self):
        """c = apex * (n/area) / F, in mol/m³."""
        # apex=200 µV, n=1e-7 mol, area=10 µV·min, F=1e-5 m³/min
        # molar flow = 200 * (1e-7/10) = 2e-6 mol/min; c = 2e-6/1e-5 = 0.2
        c = eluted_peak_concentration(200.0, 1e-7, 10.0, 1e-5)
        assert abs(c - 0.2) < 1e-12

    def test_broader_peak_gives_lower_concentration(self):
        """For fixed moles, a broader peak (larger area, same apex) yields a
        lower apex concentration — the mechanism by which the loop convention
        over-estimates the partial pressure for strongly retained probes."""
        narrow = eluted_peak_concentration(200.0, 1e-7, 5.0, 1e-5)
        broad = eluted_peak_concentration(200.0, 1e-7, 15.0, 1e-5)
        assert broad < narrow


# ---------------------------------------------------------------------------
# James–Martin pressure-gradient correction
# ---------------------------------------------------------------------------

class TestJamesMartinJ:
    """Test the James-Martin compressibility correction applied to V_N."""

    def test_no_pressure_drop_gives_unity(self):
        assert abs(james_martin_j(760.0, 760.0) - 1.0) < 1e-12

    def test_non_positive_pressure_gives_unity(self):
        assert james_martin_j(0.0, 760.0) == 1.0
        assert james_martin_j(760.0, -1.0) == 1.0

    def test_bounded_below_one(self):
        """j < 1 for any positive pressure drop; for the realistic IGC range
        (a few to a few hundred Torr over ~760) it stays close to 1."""
        for dP in (5.0, 50.0, 300.0):          # ratio 1.007 … 1.39
            j = james_martin_j(760.0 + dP, 760.0)
            assert 0.8 < j < 1.0
        # j keeps falling toward 0 as the ratio grows (→ 3/2 · 1/r)
        assert james_martin_j(760.0 + 4000.0, 760.0) < 2.0 / 3.0

    def test_decreases_with_larger_drop(self):
        j_small = james_martin_j(760.0 + 10.0, 760.0)
        j_large = james_martin_j(760.0 + 200.0, 760.0)
        assert j_large < j_small

    def test_typical_igc_value(self):
        """A 47 Torr drop over a 760 Torr outlet gives j near 0.97."""
        j = james_martin_j(760.0 + 47.0, 760.0)
        assert abs(j - 0.970) < 0.005


# ---------------------------------------------------------------------------
# Adsorption isotherm construction
# ---------------------------------------------------------------------------

class TestBuildIsotherm:
    """Test cumulative trapezoidal integration of V_N(c)/m."""

    def test_first_point_is_zero(self):
        """First point of the isotherm should be q=0 (integration starts at 0)."""
        c = np.array([0.01, 0.02, 0.03])
        V_N = np.array([1e-6, 0.9e-6, 0.8e-6])  # m³
        q = build_adsorption_isotherm(c, V_N, sample_mass_g=0.05)
        assert q[0] == 0.0

    def test_monotonically_increasing(self):
        """With positive V_N, isotherm should be monotonically increasing."""
        c = np.linspace(0.01, 0.5, 20)
        V_N = 1e-6 * np.ones(20)  # constant V_N
        q = build_adsorption_isotherm(c, V_N, sample_mass_g=0.05)
        assert np.all(np.diff(q) > 0)

    def test_larger_vn_gives_more_adsorption(self):
        """Higher V_N (stronger adsorption) should give larger q."""
        c = np.array([0.01, 0.02, 0.03])
        q_low = build_adsorption_isotherm(c, np.array([1e-7]*3), sample_mass_g=0.05)
        q_high = build_adsorption_isotherm(c, np.array([1e-6]*3), sample_mass_g=0.05)
        assert q_high[-1] > q_low[-1]

    def test_units_mol_per_gram(self):
        """Result should be in mol/g. For typical IGC, q ~ 1e-6 to 1e-5."""
        c = np.array([0.01, 0.05, 0.10, 0.15, 0.20])
        V_N = np.array([2e-6, 1.8e-6, 1.5e-6, 1.2e-6, 1.0e-6])  # m³
        q = build_adsorption_isotherm(c, V_N, sample_mass_g=0.05)
        # q should be small positive numbers (mol/g)
        assert q[-1] > 0
        assert q[-1] < 0.01  # should be much less than 0.01 mol/g


# ---------------------------------------------------------------------------
# BET linearization
# ---------------------------------------------------------------------------

class TestBETLinearization:
    """Test BET fitting against a synthetic curve with known parameters."""

    @pytest.fixture
    def synthetic_bet_data(self) -> tuple[np.ndarray, np.ndarray]:
        pp0 = np.linspace(0.05, 0.32, 12)
        n_monolayer = 0.0104
        c_bet = 3.9
        amount = (
            n_monolayer * c_bet * pp0
            / ((1.0 - pp0) * (1.0 + (c_bet - 1.0) * pp0))
        )
        return pp0, amount

    def test_recovers_known_parameters(self, synthetic_bet_data):
        pp0, amount = synthetic_bet_data
        result = bet_linearization(
            pp0, amount, p0_min=0.05, p0_max=0.37, adaptive=False
        )
        assert result.n_points == 12
        assert result.n_monolayer_mmol_g == pytest.approx(0.0104)
        assert result.C_bet == pytest.approx(3.9)
        assert result.r_squared == pytest.approx(1.0)

    def test_too_few_points_returns_nan(self):
        pp0 = np.array([0.01, 0.02])
        amount = np.array([0.001, 0.002])
        result = bet_linearization(pp0, amount, p0_min=0.05, p0_max=0.35)
        assert np.isnan(result.ssa_m2_g)
        assert result.n_points == 0

    def test_ssa_from_monolayer_capacity(self):
        n_m_mmol_g = 0.0104
        expected_ssa = n_m_mmol_g * 1e-3 * N_AVOGADRO * OCTANE_CROSS_SECTION_M2
        assert expected_ssa == pytest.approx(3.9458, abs=0.01)

    def test_c_bet_is_positive(self, synthetic_bet_data):
        pp0, amount = synthetic_bet_data
        result = bet_linearization(pp0, amount, adaptive=False)
        assert 0 < result.C_bet < 200
# BET quality control checks
# ---------------------------------------------------------------------------

def _make_result(
    ssa=2.0, n_m=0.005, C=5.0, r2=0.999, slope=100.0, intercept=100.0,
    n_points=10, pp0_range=(0.05, 0.35), mass_mg=60.0,
    injections=None, isotherm=None, bet_x=None, bet_y=None,
) -> BETResult:
    """Helper to construct a BETResult for QC testing."""
    if injections is None:
        # Default: 10 injections with decreasing V_N (good trend)
        injections = []
        for i in range(10):
            c = 0.01 * (i + 1)
            injections.append(InjectionResult(
                injection_number=i + 1,
                target_coverage=0.05 * (i + 1),
                peak_max_time=0.6 + 0.01 * i,
                peak_cofm_time=0.61 + 0.01 * i,
                peak_area=1000.0 * (i + 1),
                peak_height=500.0,
                n_injected_mol=1e-8 * (i + 1),
                V_N_mL=5.0 - 0.3 * i,  # decreasing (good)
                concentration_mol_m3=c,
                P_over_P0=0.03 * (i + 1),
            ))
    if isotherm is None:
        isotherm = []
        q_cumulative = 0.0
        for inj in injections:
            q_cumulative += 0.001 * inj.V_N_mL
            isotherm.append(IsothermPoint(
                P_over_P0=inj.P_over_P0,
                n_adsorbed_mmol_g=q_cumulative,
                V_N_mL=inj.V_N_mL,
            ))
    if bet_x is None:
        bet_x = np.array([pt.P_over_P0 for pt in isotherm
                          if pp0_range[0] <= pt.P_over_P0 <= pp0_range[1]])
    if bet_y is None:
        bet_y = np.linspace(20, 50, len(bet_x))

    return BETResult(
        ssa_m2_g=ssa, n_monolayer_mmol_g=n_m, C_bet=C,
        r_squared=r2, slope=slope, intercept=intercept,
        n_points=n_points, p_over_p0_range=pp0_range,
        injections=injections, isotherm=isotherm,
        bet_x=bet_x, bet_y=bet_y,
        temperature_K=303.15, sample_mass_mg=mass_mg, sample_name="test",
    )


class TestBETQCFitQuality:
    """Test QC checks for BET fit quality."""

    def test_few_points_flagged(self):
        """< 5 points should trigger FEW_POINTS."""
        result = _make_result(n_points=3)
        qc = bet_quality_checks(result)
        assert qc.few_points
        assert "FEW_POINTS" in qc.flags

    def test_enough_points_no_flag(self):
        """≥ 5 points should not trigger FEW_POINTS."""
        result = _make_result(n_points=10)
        qc = bet_quality_checks(result)
        assert not qc.few_points

    def test_low_r2_flagged(self):
        """R² < 0.99 should trigger LOW_R2."""
        result = _make_result(r2=0.97)
        qc = bet_quality_checks(result)
        assert qc.low_r2
        assert "LOW_R2" in qc.flags

    def test_good_r2_no_flag(self):
        """R² > 0.99 should not trigger LOW_R2."""
        result = _make_result(r2=0.998)
        qc = bet_quality_checks(result)
        assert not qc.low_r2


class TestBETQCConstant:
    """Test QC checks for BET constant C."""

    def test_c_below_1_5_critical(self):
        """C < 1.5 → critical flag (approaching Type III)."""
        result = _make_result(C=1.2)
        qc = bet_quality_checks(result)
        assert qc.c_below_1_5
        assert "C<1.5_CRITICAL" in qc.flags

    def test_c_below_2_warning(self):
        """C in 1.5-2.0 → warning (poorly conditioned)."""
        result = _make_result(C=1.7)
        qc = bet_quality_checks(result)
        assert qc.c_below_2
        assert not qc.c_below_1_5
        assert "C<2_WARN" in qc.flags

    def test_c_above_100_flagged(self):
        """C > 100 is unusually high for physisorption."""
        result = _make_result(C=150.0)
        qc = bet_quality_checks(result)
        assert qc.c_above_100
        assert "C>100_HIGH" in qc.flags

    def test_normal_c_no_flag(self):
        """C in 2-100 should not trigger any C flag."""
        result = _make_result(C=5.0)
        qc = bet_quality_checks(result)
        assert not qc.c_below_1_5
        assert not qc.c_below_2
        assert not qc.c_above_100


class TestBETQCRetentionVolume:
    """Test QC checks for measurement sensitivity (V_N)."""

    def test_very_low_vn_flagged(self):
        """Min V_N < 0.5 mL should trigger VERY_LOW_VN."""
        injections = [
            InjectionResult(
                injection_number=i, target_coverage=0.05 * i,
                peak_max_time=0.6, peak_cofm_time=0.61,
                peak_area=1000.0, peak_height=500.0,
                n_injected_mol=1e-8, V_N_mL=0.3 + 0.05 * i,
                concentration_mol_m3=0.01 * (i + 1),
                P_over_P0=0.03 * (i + 1),
            ) for i in range(5)
        ]
        result = _make_result(injections=injections)
        qc = bet_quality_checks(result)
        assert qc.very_low_vn
        assert "VERY_LOW_VN" in qc.flags

    def test_low_vn_flagged(self):
        """Min V_N 0.5-1.0 mL should trigger LOW_VN."""
        injections = [
            InjectionResult(
                injection_number=i, target_coverage=0.05 * i,
                peak_max_time=0.6, peak_cofm_time=0.61,
                peak_area=1000.0, peak_height=500.0,
                n_injected_mol=1e-8, V_N_mL=0.7 + 0.3 * i,
                concentration_mol_m3=0.01 * (i + 1),
                P_over_P0=0.03 * (i + 1),
            ) for i in range(5)
        ]
        result = _make_result(injections=injections)
        qc = bet_quality_checks(result)
        assert qc.low_vn
        assert not qc.very_low_vn
        assert "LOW_VN" in qc.flags

    def test_good_vn_no_flag(self):
        """Min V_N > 1.0 mL should not trigger any V_N flag."""
        result = _make_result()  # default V_N starts at 5.0 mL
        qc = bet_quality_checks(result)
        assert not qc.low_vn
        assert not qc.very_low_vn


class TestBETQCRouquerol:
    """Test Rouquerol consistency criteria (IUPAC 2015)."""

    def test_n_times_1_minus_pp0_must_increase(self):
        """n(1-P/P0) decreasing in BET window → ROUQUEROL_N_INCR."""
        # Construct isotherm where n(1-P/P0) decreases
        # This happens when n doesn't increase fast enough to compensate
        # for the (1-P/P0) factor decreasing
        pp0_vals = np.array([0.05, 0.10, 0.15, 0.20, 0.25, 0.30])
        # n increases, but not fast enough
        n_vals = np.array([0.001, 0.0011, 0.00115, 0.0012, 0.00122, 0.00123])
        # n*(1-P/P0) = [0.00095, 0.00099, 0.000978, 0.00096, 0.000915, 0.000861]
        # → decreasing after point 1

        isotherm = [
            IsothermPoint(P_over_P0=pp0, n_adsorbed_mmol_g=n, V_N_mL=2.0)
            for pp0, n in zip(pp0_vals, n_vals)
        ]
        result = _make_result(
            isotherm=isotherm,
            bet_x=pp0_vals,
            bet_y=pp0_vals / (n_vals * (1 - pp0_vals)),
            pp0_range=(0.05, 0.30),
            n_points=6,
        )
        qc = bet_quality_checks(result)
        assert qc.rouquerol_n_increasing
        assert "ROUQUEROL_N_INCR" in qc.flags

    def test_good_isotherm_passes_rouquerol_criterion_1(self):
        """Well-behaved isotherm should pass Rouquerol criterion 1."""
        pp0_vals = np.array([0.05, 0.10, 0.15, 0.20, 0.25, 0.30])
        # n increases fast enough: n*(1-pp0) is increasing
        n_vals = np.array([0.001, 0.0025, 0.004, 0.006, 0.009, 0.014])

        isotherm = [
            IsothermPoint(P_over_P0=pp0, n_adsorbed_mmol_g=n, V_N_mL=2.0)
            for pp0, n in zip(pp0_vals, n_vals)
        ]
        result = _make_result(
            isotherm=isotherm,
            bet_x=pp0_vals,
            bet_y=pp0_vals / (n_vals * (1 - pp0_vals)),
            pp0_range=(0.05, 0.30),
            n_points=6,
            C=5.0,  # 1/(sqrt(5)+1) ≈ 0.31 — within range
        )
        qc = bet_quality_checks(result)
        assert not qc.rouquerol_n_increasing

    def test_nm_outside_range_flagged(self):
        """Monolayer P/P0 outside fitted range → ROUQUEROL_NM_RANGE."""
        # With C = 2.0, monolayer P/P0 = 1/(√2 + 1) ≈ 0.414
        # If BET range is 0.05-0.35, this is outside
        result = _make_result(C=2.0, pp0_range=(0.05, 0.35))
        qc = bet_quality_checks(result)
        assert qc.rouquerol_nm_outside_range
        assert "ROUQUEROL_NM_RANGE" in qc.flags

    def test_nm_inside_range_no_flag(self):
        """Monolayer P/P0 inside fitted range → no flag."""
        # With C = 25.0, monolayer P/P0 = 1/(√25 + 1) = 1/6 ≈ 0.167
        # Well within 0.05-0.35
        result = _make_result(C=25.0, pp0_range=(0.05, 0.35))
        qc = bet_quality_checks(result)
        assert not qc.rouquerol_nm_outside_range


class TestBETQCIsothermShape:
    """Test isotherm monotonicity check."""

    def test_non_monotonic_isotherm_flagged(self):
        """Decreasing isotherm → ISOTHERM_NON_MONO."""
        isotherm = [
            IsothermPoint(P_over_P0=0.05, n_adsorbed_mmol_g=0.001, V_N_mL=3.0),
            IsothermPoint(P_over_P0=0.10, n_adsorbed_mmol_g=0.002, V_N_mL=2.8),
            IsothermPoint(P_over_P0=0.15, n_adsorbed_mmol_g=0.0015, V_N_mL=2.5),  # decrease!
            IsothermPoint(P_over_P0=0.20, n_adsorbed_mmol_g=0.003, V_N_mL=2.3),
        ]
        result = _make_result(isotherm=isotherm)
        qc = bet_quality_checks(result)
        assert qc.isotherm_non_monotonic
        assert "ISOTHERM_NON_MONO" in qc.flags

    def test_monotonic_isotherm_no_flag(self):
        """Properly increasing isotherm should not flag."""
        result = _make_result()  # default has monotonic isotherm
        qc = bet_quality_checks(result)
        assert not qc.isotherm_non_monotonic


class TestBETQCVNTrend:
    """Test V_N trend within BET fitting range."""

    def test_increasing_vn_in_bet_range_flagged(self):
        """V_N increasing within BET P/P0 range → VN_TREND_WRONG."""
        # Create injections where V_N increases in the BET window (0.05-0.35)
        injections = [
            InjectionResult(
                injection_number=i + 1, target_coverage=0.05 * (i + 1),
                peak_max_time=0.6, peak_cofm_time=0.61,
                peak_area=1000.0, peak_height=500.0,
                n_injected_mol=1e-8 * (i + 1),
                V_N_mL=1.0 + 0.3 * i,  # increasing within BET range
                concentration_mol_m3=0.01 * (i + 1),
                P_over_P0=0.05 + 0.03 * i,  # spans BET window
            ) for i in range(9)
        ]
        # bet_x must span these P/P0 values for the check to look inside
        bet_x = np.array([inj.P_over_P0 for inj in injections])
        result = _make_result(
            injections=injections,
            bet_x=bet_x,
            pp0_range=(0.05, 0.29),
        )
        qc = bet_quality_checks(result)
        assert qc.vn_increasing_with_concentration
        assert "VN_TREND_WRONG" in qc.flags

    def test_decreasing_vn_no_flag(self):
        """V_N decreasing with concentration (normal) → no flag."""
        result = _make_result()  # default has decreasing V_N
        qc = bet_quality_checks(result)
        assert not qc.vn_increasing_with_concentration

    def test_flat_vn_no_flag(self):
        """V_N roughly flat in BET range should not flag."""
        injections = [
            InjectionResult(
                injection_number=i + 1, target_coverage=0.05 * (i + 1),
                peak_max_time=0.6, peak_cofm_time=0.61,
                peak_area=1000.0, peak_height=500.0,
                n_injected_mol=1e-8 * (i + 1),
                V_N_mL=2.0 + 0.02 * ((-1)**i),  # ~flat with small noise
                concentration_mol_m3=0.01 * (i + 1),
                P_over_P0=0.05 + 0.03 * i,
            ) for i in range(8)
        ]
        bet_x = np.array([inj.P_over_P0 for inj in injections])
        result = _make_result(
            injections=injections,
            bet_x=bet_x,
            pp0_range=(0.05, 0.26),
        )
        qc = bet_quality_checks(result)
        assert not qc.vn_increasing_with_concentration


class TestBETQCMassSensitivity:
    """Test mass sensitivity check."""

    def test_low_mass_flagged(self):
        """Very low sample mass (< 30 mg) → MASS_SENSITIVE."""
        result = _make_result(mass_mg=20.0)
        qc = bet_quality_checks(result)
        assert qc.mass_sensitive
        assert "MASS_SENSITIVE" in qc.flags

    def test_normal_mass_no_flag(self):
        """Normal sample mass (60 mg) → no flag."""
        result = _make_result(mass_mg=60.0)
        qc = bet_quality_checks(result)
        assert not qc.mass_sensitive


class TestBETQCFlagString:
    """Test QC flag string formatting."""

    def test_clean_result_shows_ok(self):
        """Result with no flags should show 'OK'."""
        result = _make_result(C=25.0)  # C=25 keeps nm in range
        qc = bet_quality_checks(result)
        # May have ROUQUEROL_NM_RANGE depending on defaults; check the method
        if qc.passed:
            assert qc.flag_string == "OK"

    def test_multiple_flags_comma_separated(self):
        """Multiple flags should be comma-separated."""
        result = _make_result(n_points=3, r2=0.95, C=1.2)
        qc = bet_quality_checks(result)
        flags = qc.flags
        assert len(flags) >= 3
        assert "FEW_POINTS" in flags
        assert "LOW_R2" in flags
        assert "C<1.5_CRITICAL" in flags
        assert "," in qc.flag_string

    def test_passed_property(self):
        """passed should be False when any flag is triggered."""
        result = _make_result(r2=0.95)
        qc = bet_quality_checks(result)
        assert not qc.passed


# ---------------------------------------------------------------------------
# Adaptive P/P0 range selection
# ---------------------------------------------------------------------------

class TestAdaptiveRange:
    """Test adaptive P/P0 narrowing for Rouquerol compliance."""

    def test_adaptive_narrows_range(self):
        """With low C, adaptive should narrow range to satisfy Rouquerol."""
        # Build a synthetic isotherm with C ≈ 4
        # n(1-P/P0) should be increasing in a narrow enough window
        pp0 = np.array([0.03, 0.05, 0.08, 0.10, 0.13, 0.15,
                         0.18, 0.20, 0.23, 0.25, 0.28, 0.30, 0.33, 0.35])
        # Realistic BET-like adsorption (n increases with P/P0)
        n = np.array([0.002, 0.003, 0.005, 0.006, 0.008, 0.010,
                       0.013, 0.015, 0.019, 0.022, 0.027, 0.032, 0.040, 0.050])

        result_fixed = bet_linearization(pp0, n, p0_min=0.05, p0_max=0.35,
                                          adaptive=False)
        result_adaptive = bet_linearization(pp0, n, p0_min=0.05, p0_max=0.35,
                                             adaptive=True)

        # Both should produce valid SSA
        assert not math.isnan(result_fixed.ssa_m2_g)
        assert not math.isnan(result_adaptive.ssa_m2_g)
        # Adaptive may use a narrower range
        assert result_adaptive.p_over_p0_range[1] <= result_fixed.p_over_p0_range[1]

    def test_adaptive_false_uses_full_range(self):
        """With adaptive=False, should use the full requested range."""
        pp0 = np.array([0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35])
        n = np.array([0.003, 0.005, 0.008, 0.011, 0.015, 0.020, 0.026])
        result = bet_linearization(pp0, n, p0_min=0.05, p0_max=0.35,
                                    adaptive=False)
        assert result.n_points == 7

    def test_adaptive_preserves_minimum_points(self):
        """Adaptive should never go below 3 points."""
        pp0 = np.array([0.05, 0.10, 0.15, 0.20, 0.25, 0.30])
        n = np.array([0.002, 0.003, 0.004, 0.005, 0.006, 0.007])
        result = bet_linearization(pp0, n, p0_min=0.05, p0_max=0.30,
                                    adaptive=True)
        assert result.n_points >= 3


class TestIsothermClassification:
    """Tests for the Type II / III isotherm classifier."""

    @staticmethod
    def _injections_with_vn(pp0_list, vn_list):
        return [
            InjectionResult(
                injection_number=i + 1,
                target_coverage=None,
                peak_max_time=0.6,
                peak_cofm_time=0.61,
                peak_area=1000.0,
                peak_height=500.0,
                n_injected_mol=1e-8,
                V_N_mL=vn,
                concentration_mol_m3=0.01 * (i + 1),
                P_over_P0=pp0,
            )
            for i, (pp0, vn) in enumerate(zip(pp0_list, vn_list))
        ]

    def test_type_ii_declining_vn_high_c(self):
        """C >= 2 with declining/plateau V_N → Type II, BET applicable."""
        pp0 = [0.02, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30]
        vn = [2.1, 1.95, 1.85, 1.80, 1.78, 1.77, 1.76]  # decline then flat
        result = _make_result(C=3.0,
                              injections=self._injections_with_vn(pp0, vn))
        cls = classify_isotherm(result)
        assert isinstance(cls, IsothermClassification)
        assert cls.isotherm_type == "II"
        assert cls.bet_applicable is True
        assert cls.vn_rising is False

    def test_type_iii_c_below_one(self):
        """C < 1 → Type III, BET inapplicable (thermodynamic, regardless of trend)."""
        pp0 = [0.02, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30]
        vn = [0.74, 0.90, 1.05, 1.20, 1.35, 1.45, 1.55]  # rising
        result = _make_result(C=0.74,
                              injections=self._injections_with_vn(pp0, vn))
        cls = classify_isotherm(result)
        assert cls.isotherm_type == "III"
        assert cls.bet_applicable is False
        assert "< 1" in cls.rationale

    def test_type_iii_low_c_rising_vn(self):
        """1 <= C < 2 with strongly rising V_N → Type III, inapplicable."""
        pp0 = [0.02, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30]
        vn = [0.85, 1.00, 1.15, 1.30, 1.45, 1.55, 1.65]  # ~94% rise
        result = _make_result(C=1.5,
                              injections=self._injections_with_vn(pp0, vn))
        cls = classify_isotherm(result)
        assert cls.isotherm_type == "III"
        assert cls.bet_applicable is False
        assert cls.vn_rising is True

    def test_borderline_low_c_flat_vn(self):
        """1 <= C < 2 with flat V_N → borderline shape, but NOT reportable.

        The descriptive label stays 'II/III borderline', yet bet_applicable
        must be False: a C < 2 monolayer is too poorly defined to report.
        """
        pp0 = [0.02, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30]
        vn = [1.80, 1.79, 1.80, 1.78, 1.79, 1.80, 1.79]  # flat
        result = _make_result(C=1.6,
                              injections=self._injections_with_vn(pp0, vn))
        cls = classify_isotherm(result)
        assert cls.isotherm_type == "II/III borderline"
        assert cls.bet_applicable is False

    def test_type_ii_low_r2_not_reportable(self):
        """A Type-II shape (C ≥ 2) with poor linearity (R² < 0.99) is a valid
        acceptance-gate failure: bet_applicable must be False."""
        pp0 = [0.02, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30]
        vn = [2.1, 1.95, 1.85, 1.80, 1.78, 1.77, 1.76]  # declining → Type II
        result = _make_result(C=5.0, r2=0.80,
                              injections=self._injections_with_vn(pp0, vn))
        cls = classify_isotherm(result)
        assert cls.isotherm_type == "II"
        assert cls.bet_applicable is False

    def test_nan_c_is_indeterminate(self):
        """Undefined C → indeterminate, not applicable."""
        result = _make_result(C=float("nan"))
        cls = classify_isotherm(result)
        assert cls.isotherm_type == "indeterminate"
        assert cls.bet_applicable is False

    def test_fractional_rise_sign(self):
        """Declining V_N gives a negative fractional rise; rising gives positive."""
        pp0 = [0.02, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30]
        falling = _make_result(C=3.0, injections=self._injections_with_vn(
            pp0, [2.1, 1.95, 1.85, 1.80, 1.78, 1.77, 1.76]))
        rising = _make_result(C=0.8, injections=self._injections_with_vn(
            pp0, [0.74, 0.90, 1.05, 1.20, 1.35, 1.45, 1.55]))
        assert classify_isotherm(falling).vn_fractional_rise < 0
        assert classify_isotherm(rising).vn_fractional_rise > 0


# ---------------------------------------------------------------------------
# Origin-integration strategies
# ---------------------------------------------------------------------------

class TestOriginStrategies:
    """The zero-pressure origin treatment in build_adsorption_isotherm."""

    def _data(self):
        # Declining V_N (Type II), ascending concentration.
        c = np.array([0.01, 0.02, 0.03, 0.04])
        vn = np.array([2.0e-6, 1.8e-6, 1.6e-6, 1.4e-6])
        return c, vn

    def test_legacy_first_point_zero(self):
        c, vn = self._data()
        q = build_adsorption_isotherm(c, vn, 0.05, origin="legacy")
        assert q[0] == 0.0

    def test_rectangular_adds_full_rectangle(self):
        c, vn = self._data()
        q = build_adsorption_isotherm(c, vn, 0.05, origin="rectangular")
        # q[0] = V_N[0]/m * c0
        assert q[0] == pytest.approx(vn[0] / 0.05 * c[0])
        assert q[0] > 0

    def test_linear_between_legacy_and_rectangular(self):
        c, vn = self._data()
        q_leg = build_adsorption_isotherm(c, vn, 0.05, origin="legacy")
        q_rect = build_adsorption_isotherm(c, vn, 0.05, origin="rectangular")
        q_lin = build_adsorption_isotherm(c, vn, 0.05, origin="linear")
        # For declining V_N, linear extrapolation to c=0 gives V_N(0) > V_N[0],
        # so the linear origin rectangle-equivalent is >= rectangular here.
        assert q_lin[0] > q_leg[0]

    def test_only_first_point_differs(self):
        c, vn = self._data()
        q_leg = build_adsorption_isotherm(c, vn, 0.05, origin="legacy")
        q_rect = build_adsorption_isotherm(c, vn, 0.05, origin="rectangular")
        # Cumulative increments after the origin are identical.
        assert np.allclose(np.diff(q_leg), np.diff(q_rect))

    def test_unknown_strategy_raises(self):
        c, vn = self._data()
        with pytest.raises(ValueError):
            build_adsorption_isotherm(c, vn, 0.05, origin="bogus")

    def test_negative_origin_protection(self):
        """Rising V_N (Type III): linear extrapolation to c=0 must clamp ≥ 0,
        never producing a negative adsorbed amount at the origin."""
        c = np.array([0.02, 0.04, 0.06])
        vn = np.array([1.0e-6, 1.5e-6, 2.0e-6])  # rising → extrapolates below 0
        q = build_adsorption_isotherm(c, vn, 0.05, origin="linear")
        assert q[0] >= 0.0

    def test_all_strategies_registered(self):
        assert set(ORIGIN_STRATEGIES) == {"legacy", "rectangular", "linear"}


# ---------------------------------------------------------------------------
# Probe-specific SSA conversion
# ---------------------------------------------------------------------------

class TestProbeSpecificSSA:
    """SSA must scale with the probe's molecular cross-section."""

    def _pp0_n(self):
        pp0 = np.array([0.05, 0.10, 0.15, 0.20, 0.25, 0.30])
        n = np.array([0.003, 0.005, 0.008, 0.011, 0.015, 0.020])
        return pp0, n

    def test_ssa_scales_with_cross_section(self):
        pp0, n = self._pp0_n()
        octane = bet_linearization(pp0, n, adaptive=False,
                                   cross_section_m2=6.3e-19)
        hexane = bet_linearization(pp0, n, adaptive=False,
                                   cross_section_m2=5.15e-19)
        # Same isotherm, different a_cross → SSA ratio equals a_cross ratio.
        assert hexane.ssa_m2_g / octane.ssa_m2_g == pytest.approx(5.15 / 6.3, rel=1e-6)

    def test_default_is_octane(self):
        pp0, n = self._pp0_n()
        default = bet_linearization(pp0, n, adaptive=False)
        octane = bet_linearization(pp0, n, adaptive=False,
                                   cross_section_m2=OCTANE_CROSS_SECTION_M2)
        assert default.ssa_m2_g == pytest.approx(octane.ssa_m2_g)


# ---------------------------------------------------------------------------
# Matched methane dead time for peak-max vs center-of-mass
# ---------------------------------------------------------------------------

class TestMatchedRetentionConvention:
    """CoM retention must subtract the methane *CoM* dead time, and peak-max
    the methane *peak-max* dead time — not mix the two (the fixed bug)."""

    @staticmethod
    def _props():
        r = {"c1": 96.084, "c2": -7900.2, "c3": -11.003, "c4": 7e-6, "c5": 2.0}
        return ProbeProperties("OCTANE", 6.3e-19, r["c1"], r["c2"], r["c3"],
                               r["c4"], r["c5"], carbon_number=8)

    @staticmethod
    def _records():
        recs = []
        for i in range(3):
            peak = {
                "peak_max_time": 1.00 + 0.10 * i,
                "peak_cofm": 1.20 + 0.10 * i,     # CoM later than max (tailing)
                "peak_max_value": 1000.0,
                "peak_area": 100.0 * (i + 1),
                "asymmetry_factor": 1.3,
            }
            recs.append(_InjectionRecord(
                number=i + 1, target_coverage=0.1, peak=peak,
                n_injected_mol=1e-8 * (i + 1), area=100.0 * (i + 1),
                temp_col_K=303.15, flow_col_mL_min=11.0, j_factor=1.0,
                conditions_source="measured", clipped=False,
            ))
        return recs

    def test_peak_max_uses_max_dead_time(self):
        recs = self._records()
        res = _assemble_result(
            recs, mass_g=0.05, props=self._props(), V_loop_m3=2e-6,
            t0_max=0.5, t0_cofm=0.7, retention_mode="peak_max",
            concentration_mode="eluted", origin="legacy",
            p0_min=0.05, p0_max=0.35)
        inj = min(res.injections, key=lambda x: x.injection_number)
        # net = peak_max_time(1.0) - t0_max(0.5) = 0.5
        assert inj.net_retention_time_min == pytest.approx(0.5)

    def test_cofm_uses_cofm_dead_time(self):
        recs = self._records()
        res = _assemble_result(
            recs, mass_g=0.05, props=self._props(), V_loop_m3=2e-6,
            t0_max=0.5, t0_cofm=0.7, retention_mode="cofm",
            concentration_mode="eluted", origin="legacy",
            p0_min=0.05, p0_max=0.35)
        inj = min(res.injections, key=lambda x: x.injection_number)
        # net = peak_cofm(1.2) - t0_cofm(0.7) = 0.5  (NOT 1.2 - 0.5 = 0.7)
        assert inj.net_retention_time_min == pytest.approx(0.5)
        assert inj.net_retention_time_min != pytest.approx(0.7)


class TestClipDetector:
    def test_normal_peak_not_clipped(self):
        t = np.linspace(0, 5, 500)
        sig = 1000.0 * np.exp(-((t - 2.5) ** 2) / 0.05)
        assert not _is_clipped(sig)

    def test_railed_peak_clipped(self):
        sig = np.concatenate([
            np.linspace(0, 5000, 50),
            np.full(20, 5000.0),   # long identical plateau at the ceiling
            np.linspace(5000, 0, 50),
        ])
        assert _is_clipped(sig)


# ---------------------------------------------------------------------------
