"""Tests for quality control checks on dispersive surface energy results."""

import pytest
import numpy as np
import pandas as pd

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from igc_sea.analysis.quality import (
    check_alkane_count,
    check_com_divergence,
    check_gamma_d_bounds,
    classify_profile_shape,
    check_profile_shape,
    check_vn_ordering,
    run_qc_checks,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_gamma_d(coverages, gamma_d_values, n_alkanes=3):
    """Build a minimal gamma_d DataFrame for testing."""
    return pd.DataFrame({
        "coverage": coverages,
        "gamma_d_mJm2": gamma_d_values,
        "r_squared": [0.999] * len(coverages),
        "n_alkanes": [n_alkanes] * len(coverages),
    })


def _make_vn(rows):
    """Build a minimal vn_results DataFrame.

    Each row: (target_coverage, carbon_number, VN_cofm)
    """
    return pd.DataFrame(rows, columns=["target_coverage", "carbon_number", "VN_cofm"])


# ---------------------------------------------------------------------------
# Check 1: Alkane count
# ---------------------------------------------------------------------------

class TestAlkaneCount:
    def test_three_alkanes_clean(self):
        gd = _make_gamma_d([0.05, 0.10, 0.15], [35, 33, 30], n_alkanes=3)
        assert check_alkane_count(gd) == []

    def test_two_alkanes_flagged(self):
        gd = _make_gamma_d([0.05, 0.10], [35, 33], n_alkanes=2)
        flags = check_alkane_count(gd)
        assert len(flags) == 2
        assert all(f["severity"] == "warning" for f in flags)
        assert all(f["check"] == "alkane_count" for f in flags)

    def test_one_alkane_critical(self):
        gd = _make_gamma_d([0.05], [35], n_alkanes=1)
        flags = check_alkane_count(gd)
        assert len(flags) == 1
        assert flags[0]["severity"] == "critical"

    def test_mixed_counts(self):
        gd = pd.DataFrame({
            "coverage": [0.05, 0.10, 0.15],
            "gamma_d_mJm2": [35, 33, 30],
            "r_squared": [0.999, 0.999, 0.999],
            "n_alkanes": [3, 2, 3],
        })
        flags = check_alkane_count(gd)
        assert len(flags) == 1
        assert flags[0]["coverage"] == 0.10


# ---------------------------------------------------------------------------
# Check 2: Physical bounds
# ---------------------------------------------------------------------------

class TestGammaDBounds:
    def test_normal_range(self):
        gd = _make_gamma_d([0.05, 0.10, 0.15], [35, 33, 30])
        assert check_gamma_d_bounds(gd) == []

    def test_low_warning(self):
        gd = _make_gamma_d([0.05], [12.0])
        flags = check_gamma_d_bounds(gd)
        assert len(flags) == 1
        assert flags[0]["severity"] == "warning"

    def test_extreme_low_critical(self):
        gd = _make_gamma_d([0.05], [3.0])
        flags = check_gamma_d_bounds(gd)
        assert len(flags) == 1
        assert flags[0]["severity"] == "critical"

    def test_high_warning(self):
        gd = _make_gamma_d([0.05], [90.0])
        flags = check_gamma_d_bounds(gd)
        assert len(flags) == 1
        assert flags[0]["severity"] == "warning"

    def test_extreme_high_critical(self):
        gd = _make_gamma_d([0.05], [297.0])
        flags = check_gamma_d_bounds(gd)
        assert len(flags) == 1
        assert flags[0]["severity"] == "critical"

    def test_custom_bounds(self):
        gd = _make_gamma_d([0.05], [90.0])
        # With wider bounds, 90 should be clean
        assert check_gamma_d_bounds(gd, low=10.0, high=100.0) == []

    def test_nan_skipped(self):
        gd = _make_gamma_d([0.05], [float("nan")])
        assert check_gamma_d_bounds(gd) == []


# ---------------------------------------------------------------------------
# Check 3: Profile shape
# ---------------------------------------------------------------------------

class TestProfileShape:
    def test_normal_decline(self):
        gd = _make_gamma_d([0.005, 0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.15],
                           [40, 38, 36, 35, 34, 33, 32, 31, 29])
        assert classify_profile_shape(gd) == "normal_decline"

    def test_flat(self):
        gd = _make_gamma_d([0.005, 0.01, 0.02, 0.04, 0.06],
                           [40.5, 40.2, 40.8, 40.3, 40.6])
        assert classify_profile_shape(gd) == "flat"

    def test_u_shaped(self):
        gd = _make_gamma_d([0.005, 0.01, 0.02, 0.04, 0.06, 0.08, 0.10],
                           [45, 42, 38, 36, 37, 40, 43])
        assert classify_profile_shape(gd) == "u_shaped"

    def test_increasing(self):
        gd = _make_gamma_d([0.005, 0.01, 0.02, 0.04, 0.06],
                           [30, 32, 35, 38, 41])
        assert classify_profile_shape(gd) == "increasing"

    def test_collapse(self):
        gd = _make_gamma_d([0.005, 0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.15],
                           [37, 36, 37, 38, 36, 32, 26, 19, 4])
        assert classify_profile_shape(gd) == "collapse"

    def test_insufficient(self):
        gd = _make_gamma_d([0.05, 0.10], [35, 33])
        assert classify_profile_shape(gd) == "insufficient"

    def test_check_collapse_is_critical(self):
        gd = _make_gamma_d([0.02, 0.04, 0.06, 0.10, 0.15],
                           [40, 38, 35, 20, 5])
        shape, flags = check_profile_shape(gd)
        assert shape == "collapse"
        assert len(flags) == 1
        assert flags[0]["severity"] == "critical"

    def test_check_normal_no_flags(self):
        gd = _make_gamma_d([0.02, 0.04, 0.06, 0.10, 0.15],
                           [40, 38, 36, 34, 30])
        shape, flags = check_profile_shape(gd)
        assert shape == "normal_decline"
        assert flags == []


# ---------------------------------------------------------------------------
# Check 4: V_N ordering
# ---------------------------------------------------------------------------

class TestVNOrdering:
    def test_monotonic_pass(self):
        vn = _make_vn([
            (0.05, 8, 15.0),
            (0.05, 9, 40.0),
            (0.05, 10, 110.0),
        ])
        assert check_vn_ordering(vn) == []

    def test_non_monotonic_critical(self):
        vn = _make_vn([
            (0.05, 8, 50.0),   # octane higher than nonane
            (0.05, 9, 40.0),
            (0.05, 10, 110.0),
        ])
        flags = check_vn_ordering(vn)
        crits = [f for f in flags if f["severity"] == "critical"]
        assert len(crits) >= 1

    def test_outlier_ratio_warning(self):
        vn = _make_vn([
            (0.05, 8, 15.0),
            (0.05, 9, 16.0),   # ratio = 1.07, too close
            (0.05, 10, 110.0),
        ])
        flags = check_vn_ordering(vn)
        warns = [f for f in flags if f["severity"] == "warning"]
        assert len(warns) >= 1

    def test_multiple_coverages(self):
        vn = _make_vn([
            (0.05, 8, 15.0), (0.05, 9, 40.0), (0.05, 10, 110.0),  # ok
            (0.10, 8, 12.0), (0.10, 9, 11.0), (0.10, 10, 90.0),   # nonmonotonic
        ])
        flags = check_vn_ordering(vn)
        crits = [f for f in flags if f["severity"] == "critical"]
        assert len(crits) >= 1
        # Only the 0.10 coverage should be flagged
        assert all(f["coverage"] == 0.10 for f in crits)

    def test_empty_alkanes(self):
        vn = pd.DataFrame(columns=["target_coverage", "carbon_number", "VN_cofm"])
        assert check_vn_ordering(vn) == []


# ---------------------------------------------------------------------------
# Check 5: CoM divergence
# ---------------------------------------------------------------------------

class TestComDivergence:
    def test_small_divergence_no_flag(self):
        gd = pd.DataFrame({
            "coverage": [0.05, 0.10, 0.15],
            "gamma_d_mJm2": [35.5, 33.8, 30.3],
            "gamma_d_pm_mJm2": [35, 33, 30],
            "delta_cofm_pm": [0.5, 0.8, 0.3],
        })
        assert check_com_divergence(gd) == []

    def test_large_divergence_flagged(self):
        gd = pd.DataFrame({
            "coverage": [0.05, 0.10, 0.15],
            "gamma_d_mJm2": [38, 37, 34],
            "gamma_d_pm_mJm2": [35, 33, 30],
            "delta_cofm_pm": [3.0, 4.0, 4.0],
        })
        flags = check_com_divergence(gd)
        assert len(flags) == 1
        assert flags[0]["severity"] == "warning"
        assert flags[0]["check"] == "com_divergence"

    def test_no_cofm_column(self):
        gd = _make_gamma_d([0.05, 0.10], [35, 33])
        assert check_com_divergence(gd) == []


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------

class TestRunQCChecks:
    def test_clean_sample(self):
        gd = _make_gamma_d([0.02, 0.04, 0.06, 0.08, 0.10, 0.15],
                           [40, 38, 36, 34, 33, 30])
        vn = _make_vn([
            (0.05, 8, 15.0), (0.05, 9, 40.0), (0.05, 10, 110.0),
        ])
        result = run_qc_checks(gd, vn)
        assert result["pass"] is True
        assert result["profile_shape"] == "normal_decline"
        assert "0 critical" in result["summary"]

    def test_bad_sample(self):
        gd = pd.DataFrame({
            "coverage": [0.05, 0.10, 0.15],
            "gamma_d_mJm2": [297, 253, 65],
            "r_squared": [1.0, 1.0, 1.0],
            "n_alkanes": [2, 2, 2],
        })
        vn = _make_vn([
            (0.05, 8, 50.0), (0.05, 9, 40.0),  # nonmonotonic
        ])
        result = run_qc_checks(gd, vn)
        assert result["pass"] is False
        assert result["profile_shape"] == "collapse"
        # Should have critical flags from bounds and alkane count
        crits = [f for f in result["flags"] if f["severity"] == "critical"]
        assert len(crits) >= 1

    def test_return_structure(self):
        gd = _make_gamma_d([0.05, 0.10, 0.15], [35, 33, 30])
        vn = _make_vn([])
        result = run_qc_checks(gd, vn)
        assert "flags" in result
        assert "profile_shape" in result
        assert "summary" in result
        assert "pass" in result
        assert isinstance(result["flags"], list)
        assert isinstance(result["pass"], bool)
