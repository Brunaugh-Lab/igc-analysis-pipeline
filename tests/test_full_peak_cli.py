"""Regression tests for full-peak figures and generated interpretation."""

from types import SimpleNamespace

import numpy as np
import pandas as pd

from igc_sea.cli.full_peak import (
    _plot_bracket_assignments,
    _scientific_interpretation_lines,
)


def test_fixed_block_mean_bracket_plot_accepts_missing_endpoints(tmp_path):
    assignments = pd.DataFrame({
        "block": ["block-001", "block-001"],
        "assigned_t0_min": [1.0, 1.0],
        "pre_bracket_t0_min": [None, None],
        "post_bracket_t0_min": [None, None],
        "n_injected_mol": [2e-7, 1e-7],
        "net_first_moment_min": [0.2, 0.1],
    })

    _plot_bracket_assignments(assignments, tmp_path)

    assert (tmp_path / "full_peak_methane_bracket_diagnostics.pdf").is_file()
    assert (tmp_path / "full_peak_methane_bracket_diagnostics.png").is_file()


def test_interpretation_uses_actual_ranked_models_and_avoids_overclaim():
    table = pd.DataFrame([
        {"model": "none", "r_squared": 0.60, "lodo_rmse": 0.30,
         "all_params_identifiable": True},
        {"model": "langmuir", "r_squared": 0.90, "lodo_rmse": 0.100,
         "all_params_identifiable": False},
        {"model": "freundlich", "r_squared": 0.91, "lodo_rmse": 0.101,
         "all_params_identifiable": True},
    ])
    fit = SimpleNamespace(
        cooperative=False,
        params=np.array([1.0, 1.0]),
    )

    text = "\n".join(_scientific_interpretation_lines(
        fit, table, "langmuir"))

    assert "`langmuir`" in text and "`freundlich`" in text
    assert "Henry and Freundlich" not in text
    assert "adding adsorption is required" not in text
    assert "collapses onto its linear" not in text
    assert "never approaches saturation" not in text
    assert "mass-transfer" not in text
    assert "not mechanistically diagnostic" in text
