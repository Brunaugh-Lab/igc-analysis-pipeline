"""Tests for acid-base analysis: Schultz reference line, ΔG_sp, Gutmann Ka/Kb."""

import math
import pytest
import numpy as np
import pandas as pd

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from igc_sea.analysis.acid_base import (
    schultz_parameter,
    schultz_reference_line,
    calculate_delta_g_sp,
    gutmann_ka_kb,
    leave_one_out_influence,
    acid_base_quality_checks,
    check_polar_vn_consistency,
    check_dg_sp_variability,
    van_oss_gamma_sp,
    compute_van_oss,
    VAN_OSS_PROBE_PARAMS,
    run_acid_base_analysis,
    _resolve_probe,
)
from igc_sea.constants import R_GAS, N_AVOGADRO
from igc_sea.utils import get_probe


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

T_COL = 303.15  # K — standard measurement temperature


def _make_alkane_vn(vn_octane, vn_nonane, vn_decane):
    """Build a minimal alkane DataFrame at one coverage."""
    return pd.DataFrame([
        {"solvent_name": "OCTANE", "VN": vn_octane},
        {"solvent_name": "NONANE", "VN": vn_nonane},
        {"solvent_name": "DECANE", "VN": vn_decane},
    ])


def _make_polar_vn(**kwargs):
    """Build a minimal polar probe DataFrame.  kwargs: solvent_name → VN."""
    rows = [{"solvent_name": k, "VN": v} for k, v in kwargs.items()]
    return pd.DataFrame(rows)


# Schultz parameter
# ---------------------------------------------------------------------------

class TestSchultzParameter:
    def test_octane(self):
        probe = get_probe("octane")
        x = schultz_parameter(probe["a_cross"], probe["gamma_l_d"])
        # a=63e-20, gld=21.3 mJ/m² → x = 63e-20 * sqrt(0.0213) = 9.19e-20
        assert x == pytest.approx(9.194e-20, rel=1e-3)

    def test_dcm(self):
        probe = get_probe("dichloromethane")
        x = schultz_parameter(probe["a_cross"], probe["gamma_l_d"])
        # a=24.5e-20, gld=27.6 → x = 24.5e-20 * sqrt(0.0276) = 4.07e-20
        assert x == pytest.approx(4.07e-20, rel=1e-2)

    def test_zero_gamma(self):
        assert schultz_parameter(1e-19, 0.0) == 0.0


# ---------------------------------------------------------------------------
# Schultz reference line
# ---------------------------------------------------------------------------

class TestSchultzReferenceLine:
    def test_three_alkanes_perfect_line(self):
        """Synthetic data: alkanes on a perfect line → R² = 1."""
        # Use realistic VN values (mL/g) that give a clean linear relationship
        # VN increases with carbon number (heavier alkanes retain longer)
        # RT·ln(VN) should be linear vs Schultz parameter
        vn_vals = {"OCTANE": 10.0, "NONANE": 25.0, "DECANE": 60.0}

        # Compute the actual Schultz parameters and RT·ln(VN)
        xs, ys = [], []
        for name in ["OCTANE", "NONANE", "DECANE"]:
            probe = get_probe(name.lower())
            xs.append(schultz_parameter(probe["a_cross"], probe["gamma_l_d"]))
            ys.append(R_GAS * T_COL * math.log(vn_vals[name]))

        # Create perfectly linear data by adjusting VN for decane
        # Fit line through octane and nonane, predict decane
        slope_2pt = (ys[1] - ys[0]) / (xs[1] - xs[0])
        intercept_2pt = ys[0] - slope_2pt * xs[0]
        y_decane_pred = slope_2pt * xs[2] + intercept_2pt
        vn_decane_perfect = math.exp(y_decane_pred / (R_GAS * T_COL))

        rows = [
            {"solvent_name": "OCTANE", "VN": 10.0},
            {"solvent_name": "NONANE", "VN": 25.0},
            {"solvent_name": "DECANE", "VN": vn_decane_perfect},
        ]

        df = pd.DataFrame(rows)
        result = schultz_reference_line(df, T_COL)

        assert result["n_alkanes"] == 3
        assert result["r_squared"] == pytest.approx(1.0, abs=1e-8)
        assert result["slope"] == pytest.approx(slope_2pt, rel=1e-6)

    def test_gamma_d_from_slope(self):
        """Verify γ_d extraction from Schultz slope."""
        # slope = 2 · N_A · √(γ_d) → γ_d = (slope / (2·N_A))²
        gamma_d_target = 38.0e-3  # 38 mJ/m² in J/m²
        slope = 2 * N_AVOGADRO * math.sqrt(gamma_d_target)

        rows = []
        for name in ["OCTANE", "NONANE", "DECANE"]:
            probe = get_probe(name.lower())
            x = schultz_parameter(probe["a_cross"], probe["gamma_l_d"])
            y = slope * x + 1000.0
            vn = math.exp(y / (R_GAS * T_COL))
            rows.append({"solvent_name": name, "VN": vn})

        result = schultz_reference_line(pd.DataFrame(rows), T_COL)
        assert result["gamma_d_schultz_mJm2"] == pytest.approx(38.0, rel=1e-4)

    def test_single_alkane_returns_nan(self):
        df = pd.DataFrame([{"solvent_name": "OCTANE", "VN": 10.0}])
        result = schultz_reference_line(df, T_COL)
        assert np.isnan(result["slope"])
        assert result["n_alkanes"] == 1

    def test_negative_vn_skipped(self):
        df = _make_alkane_vn(10.0, -5.0, 20.0)
        result = schultz_reference_line(df, T_COL)
        assert result["n_alkanes"] == 2  # nonane skipped


# ---------------------------------------------------------------------------
# ΔG_sp calculation
# ---------------------------------------------------------------------------

class TestDeltaGSp:
    def test_probe_above_alkane_line_positive(self):
        """A polar probe above the alkane line has positive ΔG_sp."""
        # Build a Schultz line
        alkanes = _make_alkane_vn(10.0, 25.0, 60.0)
        sline = schultz_reference_line(alkanes, T_COL)

        # DCM with VN that puts it above the line
        polar = _make_polar_vn(DICHLOROMETHANE=50.0)
        result = calculate_delta_g_sp(polar, sline, T_COL)

        assert len(result) == 1
        assert result.iloc[0]["probe"] == "dichloromethane"
        assert result.iloc[0]["delta_g_sp_kJmol"] > 0

    def test_alkanes_excluded(self):
        """Alkane probes in the input are not returned as polar results."""
        alkanes = _make_alkane_vn(10.0, 25.0, 60.0)
        sline = schultz_reference_line(alkanes, T_COL)

        # Pass alkanes as "polar" input — should be excluded
        result = calculate_delta_g_sp(alkanes, sline, T_COL)
        assert len(result) == 0

    def test_nan_vn_skipped(self):
        alkanes = _make_alkane_vn(10.0, 25.0, 60.0)
        sline = schultz_reference_line(alkanes, T_COL)

        polar = pd.DataFrame([
            {"solvent_name": "DICHLOROMETHANE", "VN": 50.0},
            {"solvent_name": "ETHYL ACETATE", "VN": np.nan},
        ])
        result = calculate_delta_g_sp(polar, sline, T_COL)
        assert len(result) == 1

    def test_nan_schultz_line_returns_empty(self):
        sline = {"slope": np.nan, "intercept": np.nan}
        polar = _make_polar_vn(DICHLOROMETHANE=50.0)
        result = calculate_delta_g_sp(polar, sline, T_COL)
        assert len(result) == 0

    def test_dn_an_star_populated(self):
        """Verify probe properties are attached to results."""
        alkanes = _make_alkane_vn(10.0, 25.0, 60.0)
        sline = schultz_reference_line(alkanes, T_COL)
        polar = _make_polar_vn(DICHLOROMETHANE=50.0)
        result = calculate_delta_g_sp(polar, sline, T_COL)

        row = result.iloc[0]
        assert row["dn"] == 0.0  # DCM is pure acid
        assert row["an_star"] == 16.3


# ---------------------------------------------------------------------------
# Gutmann Ka/Kb regression
# ---------------------------------------------------------------------------

class TestGutmannKaKb:
    def test_known_regression(self):
        """Verify Ka/Kb extraction from known ΔG_sp values."""
        # Set up: Ka=0.1, Kb=0.5
        # ΔG_sp = Ka·DN + Kb·AN*
        ka_true = 0.1
        kb_true = 0.5

        probes = [
            ("dichloromethane", 0.0, 16.4),
            ("ethyl acetate", 71.5, 6.3),
            ("acetone", 71.1, 10.5),
            ("ethanol", 80.0, 37.1),
        ]

        rows = []
        for name, dn, an_star in probes:
            dg_sp = ka_true * dn + kb_true * an_star
            rows.append({
                "probe": name,
                "delta_g_sp_kJmol": dg_sp,
                "dn": dn,
                "an_star": an_star,
            })

        result = gutmann_ka_kb(pd.DataFrame(rows))
        assert result["Ka"] == pytest.approx(ka_true, rel=1e-6)
        assert result["Kb"] == pytest.approx(kb_true, rel=1e-6)
        assert result["r_squared"] == pytest.approx(1.0, abs=1e-10)
        assert result["n_probes"] == 4
        assert result["fit_method"] == "regression"

    def test_exclude_probes(self):
        """Excluding a probe reduces n_probes."""
        rows = [
            {"probe": "dichloromethane", "delta_g_sp_kJmol": 14.0, "dn": 0.0, "an_star": 16.4},
            {"probe": "ethyl acetate", "delta_g_sp_kJmol": 9.0, "dn": 71.5, "an_star": 6.3},
            {"probe": "acetone", "delta_g_sp_kJmol": 9.5, "dn": 71.1, "an_star": 10.5},
            {"probe": "ethanol", "delta_g_sp_kJmol": 12.0, "dn": 80.0, "an_star": 37.1},
        ]
        result = gutmann_ka_kb(pd.DataFrame(rows), exclude_probes=["ethanol"])
        assert result["n_probes"] == 3
        assert "ethanol" not in result["probes_used"]

    def test_empty_input(self):
        result = gutmann_ka_kb(pd.DataFrame())
        assert np.isnan(result["Ka"])
        assert result["n_probes"] == 0

    def test_single_probe(self):
        df = pd.DataFrame([{
            "probe": "dichloromethane",
            "delta_g_sp_kJmol": 14.0,
            "dn": 0.0,
            "an_star": 16.4,
        }])
        result = gutmann_ka_kb(df)
        assert np.isnan(result["Ka"])
        assert result["fit_method"] == "insufficient"

    def test_two_probe_deterministic(self):
        """Two probes give a deterministic solve, not a regression."""
        rows = [
            {"probe": "dichloromethane", "delta_g_sp_kJmol": 14.0, "dn": 0.0, "an_star": 16.3},
            {"probe": "ethyl acetate", "delta_g_sp_kJmol": 9.0, "dn": 71.5, "an_star": 6.3},
        ]
        result = gutmann_ka_kb(pd.DataFrame(rows))
        assert result["fit_method"] == "deterministic"
        assert np.isnan(result["r_squared"])  # R² meaningless for 2 points
        assert result["n_probes"] == 2
        assert not np.isnan(result["Ka"])
        assert not np.isnan(result["Kb"])

    def test_residuals_returned(self):
        rows = [
            {"probe": "dichloromethane", "delta_g_sp_kJmol": 14.0, "dn": 0.0, "an_star": 16.4},
            {"probe": "ethyl acetate", "delta_g_sp_kJmol": 9.0, "dn": 71.5, "an_star": 6.3},
            {"probe": "acetone", "delta_g_sp_kJmol": 9.5, "dn": 71.1, "an_star": 10.5},
        ]
        result = gutmann_ka_kb(pd.DataFrame(rows))
        assert len(result["residuals"]) == 3
        assert all("residual" in r for r in result["residuals"])


# ---------------------------------------------------------------------------
# Leave-one-out influence
# ---------------------------------------------------------------------------

class TestLeaveOneOut:
    def test_identifies_outlier(self):
        """An outlier probe should be flagged."""
        # 3 probes on a perfect line, one outlier
        ka, kb = 0.07, 0.5
        rows = []
        for name, dn, an_star in [
            ("dichloromethane", 0.0, 16.4),
            ("ethyl acetate", 71.5, 6.3),
            ("acetone", 71.1, 10.5),
        ]:
            rows.append({
                "probe": name,
                "delta_g_sp_kJmol": ka * dn + kb * an_star,
                "dn": dn,
                "an_star": an_star,
            })

        # Add an outlier: ethanol with wrong ΔG_sp
        rows.append({
            "probe": "ethanol",
            "delta_g_sp_kJmol": 2.0,  # way off
            "dn": 80.0,
            "an_star": 37.1,
        })

        df = pd.DataFrame(rows)
        influence = leave_one_out_influence(df)

        ethanol_inf = [i for i in influence if i["probe"] == "ethanol"]
        assert len(ethanol_inf) == 1
        assert ethanol_inf[0]["is_outlier"]
        # Removing ethanol should improve R²
        assert ethanol_inf[0]["delta_r_squared"] > 0

    def test_too_few_probes(self):
        df = pd.DataFrame([
            {"probe": "dcm", "delta_g_sp_kJmol": 14.0, "dn": 0.0, "an_star": 16.4},
            {"probe": "etac", "delta_g_sp_kJmol": 9.0, "dn": 71.5, "an_star": 6.3},
        ])
        result = leave_one_out_influence(df)
        assert result == []


# ---------------------------------------------------------------------------
# Quality checks
# ---------------------------------------------------------------------------

class TestQualityChecks:
    def test_low_r2_flag(self):
        gutmann = {"Ka": 0.07, "Kb": 0.5, "r_squared": 0.4, "n_probes": 4}
        dg_sp = pd.DataFrame([
            {"probe": "dcm", "delta_g_sp_kJmol": 14.0},
        ])
        qc = acid_base_quality_checks(gutmann, dg_sp)
        codes = [f["code"] for f in qc["flags"]]
        assert "GUTMANN_R2_LOW" in codes

    def test_critical_r2_flag(self):
        gutmann = {"Ka": 0.07, "Kb": 0.5, "r_squared": 0.2, "n_probes": 4}
        dg_sp = pd.DataFrame({"probe": [], "delta_g_sp_kJmol": []})
        qc = acid_base_quality_checks(gutmann, dg_sp)
        codes = [f["code"] for f in qc["flags"]]
        assert "GUTMANN_R2_CRITICAL" in codes
        assert "FAIL" in qc["summary"]

    def test_few_probes_flag(self):
        """1 probe triggers FEW_PROBES."""
        gutmann = {"Ka": np.nan, "Kb": np.nan, "r_squared": np.nan,
                    "n_probes": 1, "fit_method": "insufficient"}
        dg_sp = pd.DataFrame({"probe": [], "delta_g_sp_kJmol": []})
        qc = acid_base_quality_checks(gutmann, dg_sp)
        codes = [f["code"] for f in qc["flags"]]
        assert "FEW_PROBES" in codes

    def test_deterministic_fit_flag(self):
        """2-probe fit triggers DETERMINISTIC_FIT, not FEW_PROBES."""
        gutmann = {"Ka": 0.07, "Kb": 0.5, "r_squared": np.nan,
                    "n_probes": 2, "fit_method": "deterministic"}
        dg_sp = pd.DataFrame({"probe": [], "delta_g_sp_kJmol": []})
        qc = acid_base_quality_checks(gutmann, dg_sp)
        codes = [f["code"] for f in qc["flags"]]
        assert "DETERMINISTIC_FIT" in codes
        assert "FEW_PROBES" not in codes

    def test_negative_dg_sp_flag(self):
        gutmann = {"Ka": 0.07, "Kb": 0.5, "r_squared": 0.9, "n_probes": 4}
        dg_sp = pd.DataFrame([
            {"probe": "ethanol", "delta_g_sp_kJmol": -1.5},
        ])
        qc = acid_base_quality_checks(gutmann, dg_sp)
        codes = [f["code"] for f in qc["flags"]]
        assert "NEGATIVE_DG_SP" in codes

    def test_pass_summary(self):
        gutmann = {"Ka": 0.07, "Kb": 0.5, "r_squared": 0.85, "n_probes": 4}
        dg_sp = pd.DataFrame([
            {"probe": "dcm", "delta_g_sp_kJmol": 14.0},
        ])
        qc = acid_base_quality_checks(gutmann, dg_sp)
        assert qc["summary"] == "PASS"

    def test_negative_ka_flag(self):
        gutmann = {"Ka": -0.01, "Kb": 0.5, "r_squared": 0.9, "n_probes": 4}
        dg_sp = pd.DataFrame({"probe": [], "delta_g_sp_kJmol": []})
        qc = acid_base_quality_checks(gutmann, dg_sp)
        codes = [f["code"] for f in qc["flags"]]
        assert "NEGATIVE_KA" in codes


# ---------------------------------------------------------------------------
# Probe resolution
# ---------------------------------------------------------------------------

class TestProbeResolution:
    def test_uppercase_solvent(self):
        probe = _resolve_probe("DICHLOROMETHANE")
        assert probe is not None
        assert probe["name"] == "dichloromethane"

    def test_ethanol(self):
        probe = _resolve_probe("ETHANOL")
        assert probe is not None
        assert probe["dn"] == 80.0
        assert probe["an_star"] == 37.1

    def test_unknown_returns_none(self):
        assert _resolve_probe("WATER") is None

    def test_alkane_resolves(self):
        probe = _resolve_probe("octane")
        assert probe is not None
        assert probe["category"] == "alkane"


# ---------------------------------------------------------------------------
# Full orchestrator
# ---------------------------------------------------------------------------

class TestRunAcidBaseAnalysis:
    def test_multi_coverage(self):
        """Verify analysis runs across multiple coverages."""
        # Build synthetic data at two coverages
        rows = []
        for cov in [0.01, 0.02]:
            for name, vn_base in [("OCTANE", 10.0), ("NONANE", 25.0), ("DECANE", 60.0)]:
                rows.append({"solvent_name": name, "VN": vn_base * (1 - cov), "coverage": cov})
            for name, vn_base in [("DICHLOROMETHANE", 50.0), ("ETHYL ACETATE", 15.0),
                                   ("ACETONE", 12.0)]:
                rows.append({"solvent_name": name, "VN": vn_base * (1 - cov), "coverage": cov})

        df = pd.DataFrame(rows)
        result = run_acid_base_analysis(df, T_COL)

        assert "profile" in result
        assert "delta_g_sp" in result
        assert "qc" in result
        assert len(result["profile"]) == 2  # two coverages

    def test_no_polar_probes(self):
        """Only alkanes → empty acid-base results."""
        df = pd.DataFrame([
            {"solvent_name": "OCTANE", "VN": 10.0, "coverage": 0.01},
            {"solvent_name": "NONANE", "VN": 25.0, "coverage": 0.01},
            {"solvent_name": "DECANE", "VN": 60.0, "coverage": 0.01},
        ])
        result = run_acid_base_analysis(df, T_COL)
        # Profile should have a row but no probes for Gutmann
        assert result["delta_g_sp"].empty

    def test_fit_method_in_profile(self):
        """Profile DataFrame includes fit_method column."""
        rows = []
        for cov in [0.01]:
            for name, vn in [("OCTANE", 10.0), ("NONANE", 25.0), ("DECANE", 60.0)]:
                rows.append({"solvent_name": name, "VN": vn, "coverage": cov})
            for name, vn in [("DICHLOROMETHANE", 50.0), ("ETHYL ACETATE", 15.0),
                              ("ACETONE", 12.0)]:
                rows.append({"solvent_name": name, "VN": vn, "coverage": cov})
        result = run_acid_base_analysis(pd.DataFrame(rows), T_COL)
        assert "fit_method" in result["profile"].columns
        assert result["profile"].iloc[0]["fit_method"] == "regression"

    def test_two_probe_sample_van_oss(self):
        """Sample with only DCM + EtAc gets van Oss results."""
        rows = []
        for cov in [0.01]:
            for name, vn in [("OCTANE", 10.0), ("NONANE", 25.0), ("DECANE", 60.0)]:
                rows.append({"solvent_name": name, "VN": vn, "coverage": cov})
            rows.append({"solvent_name": "DICHLOROMETHANE", "VN": 50.0, "coverage": cov})
            rows.append({"solvent_name": "ETHYL ACETATE", "VN": 15.0, "coverage": cov})
        # Without include_van_oss, van Oss should not be computed
        result_no_vo = run_acid_base_analysis(pd.DataFrame(rows), T_COL)
        assert result_no_vo["profile"].iloc[0]["fit_method"] == "deterministic"
        assert len(result_no_vo["van_oss_results"]) == 0
        assert "gamma_s_minus_mJm2" not in result_no_vo["profile"].columns

        # With include_van_oss=True, van Oss should be computed
        result = run_acid_base_analysis(
            pd.DataFrame(rows), T_COL, include_van_oss=True,
        )
        assert result["profile"].iloc[0]["fit_method"] == "deterministic"
        assert len(result["van_oss_results"]) == 1
        assert not np.isnan(result["profile"].iloc[0]["gamma_s_minus_mJm2"])


# ---------------------------------------------------------------------------
# Probe categories
# ---------------------------------------------------------------------------

class TestProbeCategories:
    def test_ethyl_acetate_is_base(self):
        assert get_probe("ethyl acetate")["category"] == "base"

    def test_acetone_is_base(self):
        assert get_probe("acetone")["category"] == "base"

    def test_ethanol_is_amphoteric(self):
        assert get_probe("ethanol")["category"] == "amphoteric"

    def test_dcm_is_acid(self):
        assert get_probe("dichloromethane")["category"] == "acid"

    def test_acetonitrile_is_amphoteric(self):
        assert get_probe("acetonitrile")["category"] == "amphoteric"


# ---------------------------------------------------------------------------
# Polar V_N consistency
# ---------------------------------------------------------------------------

class TestPolarVNConsistency:
    def test_erratic_probe_flagged(self):
        """V_N varying >10× non-monotonically is flagged."""
        rows = []
        vn_values = [1800, 99, 467, 83, 35, 2072, 57]
        for i, vn in enumerate(vn_values):
            rows.append({
                "solvent_name": "ACETONITRILE",
                "VN": float(vn),
                "coverage": 0.005 + i * 0.02,
            })
        checks = check_polar_vn_consistency(pd.DataFrame(rows))
        assert len(checks) == 1
        assert checks[0]["flagged"]
        assert checks[0]["vn_ratio"] > 10
        assert not checks[0]["is_monotonic"]
        assert "acetonitrile" in checks[0]["message"]

    def test_monotonic_not_flagged(self):
        """V_N decreasing monotonically is not flagged even with large ratio."""
        rows = [
            {"solvent_name": "ACETONITRILE", "VN": 2000.0, "coverage": 0.005},
            {"solvent_name": "ACETONITRILE", "VN": 500.0, "coverage": 0.01},
            {"solvent_name": "ACETONITRILE", "VN": 100.0, "coverage": 0.02},
        ]
        checks = check_polar_vn_consistency(pd.DataFrame(rows))
        assert len(checks) == 1
        assert not checks[0]["flagged"]

    def test_small_ratio_not_flagged(self):
        """V_N varying <10× is not flagged."""
        rows = [
            {"solvent_name": "ETHYL ACETATE", "VN": 20.0, "coverage": 0.005},
            {"solvent_name": "ETHYL ACETATE", "VN": 25.0, "coverage": 0.01},
            {"solvent_name": "ETHYL ACETATE", "VN": 15.0, "coverage": 0.02},
        ]
        checks = check_polar_vn_consistency(pd.DataFrame(rows))
        assert len(checks) == 1
        assert not checks[0]["flagged"]

    def test_alkanes_excluded(self):
        """Alkane probes are not checked."""
        rows = [
            {"solvent_name": "OCTANE", "VN": 2000.0, "coverage": 0.005},
            {"solvent_name": "OCTANE", "VN": 50.0, "coverage": 0.01},
            {"solvent_name": "OCTANE", "VN": 1500.0, "coverage": 0.02},
        ]
        checks = check_polar_vn_consistency(pd.DataFrame(rows))
        assert len(checks) == 0


# ---------------------------------------------------------------------------
# Van Oss acid-base components
# ---------------------------------------------------------------------------

class TestVanOss:
    def test_known_gamma_roundtrip(self):
        """Known γ_S⁺/γ_S⁻ → ΔG_sp → van_oss_gamma_sp → recover originals."""
        gamma_s_plus = 0.5e-3   # J/m² = 0.5 mJ/m²
        gamma_s_minus = 15.0e-3  # J/m² = 15.0 mJ/m²

        dcm = get_probe("dichloromethane")
        etac = get_probe("ethyl acetate")
        dcm_vo = VAN_OSS_PROBE_PARAMS["dichloromethane"]
        etac_vo = VAN_OSS_PROBE_PARAMS["ethyl acetate"]

        # ΔG_sp = 2 · a · N_A · √(γ_L⁺ · γ_S⁻)  for acid probe
        dg_dcm = 2.0 * dcm["a_cross"] * N_AVOGADRO * np.sqrt(
            dcm_vo["gamma_l_plus"] * 1e-3 * gamma_s_minus)
        # ΔG_sp = 2 · a · N_A · √(γ_L⁻ · γ_S⁺)  for base probe
        dg_etac = 2.0 * etac["a_cross"] * N_AVOGADRO * np.sqrt(
            etac_vo["gamma_l_minus"] * 1e-3 * gamma_s_plus)

        result = van_oss_gamma_sp(
            delta_g_sp_acid=dg_dcm,
            delta_g_sp_base=dg_etac,
            a_acid=dcm["a_cross"],
            a_base=etac["a_cross"],
            gamma_l_plus_acid=dcm_vo["gamma_l_plus"],
            gamma_l_minus_base=etac_vo["gamma_l_minus"],
        )

        assert result["gamma_s_minus_mJm2"] == pytest.approx(15.0, rel=1e-6)
        assert result["gamma_s_plus_mJm2"] == pytest.approx(0.5, rel=1e-6)
        assert result["gamma_s_ab_mJm2"] == pytest.approx(
            2.0 * np.sqrt(0.5 * 15.0), rel=1e-6)

    def test_zero_dg_sp_gives_zero_gamma(self):
        """ΔG_sp = 0 → γ component = 0 (no interaction)."""
        result = van_oss_gamma_sp(
            delta_g_sp_acid=0.0,
            delta_g_sp_base=0.0,
            a_acid=24.5e-20,
            a_base=33.0e-20,
            gamma_l_plus_acid=5.2,
            gamma_l_minus_base=19.2,
        )
        assert result["gamma_s_minus_mJm2"] == 0.0
        assert result["gamma_s_plus_mJm2"] == 0.0
        assert result["gamma_s_ab_mJm2"] == 0.0

    def test_compute_van_oss_with_both_probes(self):
        """compute_van_oss returns results when DCM + EtAc are present."""
        alkanes = _make_alkane_vn(10.0, 25.0, 60.0)
        sline = schultz_reference_line(alkanes, T_COL)
        polar = pd.DataFrame([
            {"solvent_name": "DICHLOROMETHANE", "VN": 50.0},
            {"solvent_name": "ETHYL ACETATE", "VN": 15.0},
        ])
        dg_sp = calculate_delta_g_sp(polar, sline, T_COL)
        result = compute_van_oss(dg_sp)
        assert result is not None
        assert "gamma_s_minus_mJm2" in result
        assert "gamma_s_plus_mJm2" in result
        assert result["gamma_s_minus_mJm2"] > 0  # DCM ΔG_sp > 0 → surface is basic

    def test_compute_van_oss_missing_probe(self):
        """compute_van_oss returns None if either probe is missing."""
        alkanes = _make_alkane_vn(10.0, 25.0, 60.0)
        sline = schultz_reference_line(alkanes, T_COL)
        polar = pd.DataFrame([{"solvent_name": "DICHLOROMETHANE", "VN": 50.0}])
        dg_sp = calculate_delta_g_sp(polar, sline, T_COL)
        assert compute_van_oss(dg_sp) is None

    def test_W_a_sp_returned(self):
        """Van Oss result includes work of adhesion per unit area."""
        result = van_oss_gamma_sp(
            delta_g_sp_acid=14000.0,
            delta_g_sp_base=8000.0,
            a_acid=24.5e-20,
            a_base=33.0e-20,
            gamma_l_plus_acid=5.2,
            gamma_l_minus_base=19.2,
        )
        assert "W_a_acid_mJm2" in result
        assert "W_a_base_mJm2" in result
        # W_a = ΔG_sp / (N_A * a), should be ~95 mJ/m² for DCM
        assert result["W_a_acid_mJm2"] == pytest.approx(
            14000.0 / (N_AVOGADRO * 24.5e-20) * 1e3, rel=1e-6)


# ---------------------------------------------------------------------------
# Van Oss QC flags
# ---------------------------------------------------------------------------

class TestVanOssQC:
    def test_high_gamma_flagged(self):
        """γ_S⁻ > 100 mJ/m² triggers VAN_OSS_GAMMA_HIGH."""
        gutmann = {"Ka": 0.07, "Kb": 0.5, "r_squared": 0.9,
                    "n_probes": 4, "fit_method": "regression"}
        dg_sp = pd.DataFrame({"probe": ["dichloromethane"], "delta_g_sp_kJmol": [14.0]})
        vo = {"gamma_s_minus_mJm2": 433.0, "gamma_s_plus_mJm2": 24.0,
              "gamma_s_ab_mJm2": 200.0}
        qc = acid_base_quality_checks(gutmann, dg_sp, van_oss_result=vo)
        codes = [f["code"] for f in qc["flags"]]
        assert "VAN_OSS_GAMMA_HIGH" in codes

    def test_reasonable_gamma_not_flagged(self):
        """γ_S⁻ < 100 mJ/m² does not trigger flag."""
        gutmann = {"Ka": 0.07, "Kb": 0.5, "r_squared": 0.9,
                    "n_probes": 4, "fit_method": "regression"}
        dg_sp = pd.DataFrame({"probe": ["dichloromethane"], "delta_g_sp_kJmol": [5.0]})
        vo = {"gamma_s_minus_mJm2": 30.0, "gamma_s_plus_mJm2": 5.0,
              "gamma_s_ab_mJm2": 24.5}
        qc = acid_base_quality_checks(gutmann, dg_sp, van_oss_result=vo)
        codes = [f["code"] for f in qc["flags"]]
        assert "VAN_OSS_GAMMA_HIGH" not in codes

    def test_incomplete_van_oss_flag(self):
        """Missing DCM triggers VAN_OSS_INCOMPLETE."""
        gutmann = {"Ka": 0.07, "Kb": 0.5, "r_squared": 0.9,
                    "n_probes": 3, "fit_method": "regression"}
        dg_sp = pd.DataFrame([
            {"probe": "ethyl acetate", "delta_g_sp_kJmol": 9.0},
            {"probe": "acetone", "delta_g_sp_kJmol": 10.0},
            {"probe": "ethanol", "delta_g_sp_kJmol": 12.0},
        ])
        qc = acid_base_quality_checks(gutmann, dg_sp, van_oss_result=None)
        codes = [f["code"] for f in qc["flags"]]
        assert "VAN_OSS_INCOMPLETE" in codes
        msg = [f["message"] for f in qc["flags"] if f["code"] == "VAN_OSS_INCOMPLETE"][0]
        assert "DCM" in msg


# ---------------------------------------------------------------------------
# ΔG_sp variability check
# ---------------------------------------------------------------------------

class TestDGSpVariability:
    def test_high_variability_flagged(self):
        """Probe with 3× median std is flagged."""
        rows = []
        # Well-behaved probes: std ~0.5
        for cov in [0.005, 0.01, 0.02, 0.04]:
            rows.append({"probe": "ethyl acetate", "delta_g_sp_kJmol": 9.0 + cov * 10, "coverage": cov})
            rows.append({"probe": "acetone", "delta_g_sp_kJmol": 10.0 + cov * 10, "coverage": cov})
        # Erratic probe: std ~7
        acn_vals = [21.0, 8.0, 18.0, 3.0]
        for i, cov in enumerate([0.005, 0.01, 0.02, 0.04]):
            rows.append({"probe": "acetonitrile", "delta_g_sp_kJmol": acn_vals[i], "coverage": cov})

        checks = check_dg_sp_variability(pd.DataFrame(rows))
        flagged = [c for c in checks if c["flagged"]]
        flagged_probes = [c["probe"] for c in flagged]
        assert "acetonitrile" in flagged_probes

    def test_uniform_variability_not_flagged(self):
        """All probes with similar std are not flagged."""
        rows = []
        for cov in [0.005, 0.01, 0.02]:
            rows.append({"probe": "ethyl acetate", "delta_g_sp_kJmol": 9.0 - cov * 20, "coverage": cov})
            rows.append({"probe": "acetone", "delta_g_sp_kJmol": 10.0 - cov * 20, "coverage": cov})
            rows.append({"probe": "ethanol", "delta_g_sp_kJmol": 12.0 - cov * 20, "coverage": cov})
        checks = check_dg_sp_variability(pd.DataFrame(rows))
        flagged = [c for c in checks if c["flagged"]]
        assert len(flagged) == 0
