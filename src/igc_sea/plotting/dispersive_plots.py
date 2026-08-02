"""Visualization for dispersive surface energy analysis.

Plots:
  - Alkane reference line (RT·ln(V_N) vs carbon number)
  - γ_d vs surface coverage heterogeneity profile
  - Comparison against an independently supplied reference profile
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from igc_sea.constants import R_GAS
from igc_sea.plotting.theme import get_palette, SINGLE_COL, DOUBLE_COL


def plot_alkane_line(
    df: pd.DataFrame,
    coverage: float,
    temperature: float,
    fit_params: dict | None = None,
    polar_probes: pd.DataFrame | None = None,
    title: str = "",
) -> plt.Figure:
    """Plot the alkane reference line at a single coverage.

    Parameters
    ----------
    df : pd.DataFrame
        Alkane data with columns ``carbon_number``, ``VN``.
    coverage : float
        Surface coverage value (for labeling).
    temperature : float
        Column temperature (K).
    fit_params : dict, optional
        Keys ``slope_Jmol``, ``intercept``, ``r_squared``, ``gamma_d_mJm2``.
    polar_probes : pd.DataFrame, optional
        Polar probe data with columns ``solvent_name``, ``VN``.
    title : str
        Plot title override.

    Returns
    -------
    plt.Figure
    """
    fig, ax = plt.subplots(figsize=(SINGLE_COL, SINGLE_COL * 0.85))
    palette = get_palette()

    # Compute RT·ln(V_N) for alkanes
    alkanes = df.dropna(subset=["VN"]).copy()
    alkanes["RT_ln_VN"] = R_GAS * temperature * np.log(alkanes["VN"])

    ax.scatter(alkanes["carbon_number"], alkanes["RT_ln_VN"] / 1000,
               color=palette[0], s=40, zorder=5, label="Alkanes")

    # Fit line
    if fit_params is not None:
        x_range = np.array([alkanes["carbon_number"].min() - 0.5,
                            alkanes["carbon_number"].max() + 0.5])
        y_fit = (fit_params["slope_Jmol"] * x_range + fit_params["intercept"]) / 1000
        ax.plot(x_range, y_fit, "--", color=palette[0], linewidth=0.8)

        # Annotation
        r2 = fit_params["r_squared"]
        gd = fit_params["gamma_d_mJm2"]
        ax.text(0.05, 0.95,
                f"R² = {r2:.6f}\nγd = {gd:.2f} mJ/m²",
                transform=ax.transAxes, fontsize=7, va="top",
                bbox=dict(facecolor="white", alpha=0.8, edgecolor="none"))

    # Polar probes (if provided)
    if polar_probes is not None and not polar_probes.empty:
        pp = polar_probes.dropna(subset=["VN"]).copy()
        # Use molecular-weight proxy for x-position (for visualization only)
        # Plot as triangles above the alkane line
        for i, (_, row) in enumerate(pp.iterrows()):
            rt_ln = R_GAS * temperature * np.log(row["VN"]) / 1000
            ax.scatter([], [], color=palette[1 + i % 6], marker="^", s=30)
            # We don't have a natural x-axis position for polar probes
            # in Dorris-Gray — skip plotting them on the carbon number axis

    if not title:
        title = f"Alkane reference line (n/nm = {coverage})"
    ax.set_title(title, fontsize=9)
    ax.set_xlabel("Carbon number")
    ax.set_ylabel("RT·ln(V_N) (kJ/mol)")
    ax.legend(fontsize=7, frameon=False)

    fig.tight_layout()
    return fig


def plot_gamma_d_profile(
    profile: pd.DataFrame,
    reference: pd.DataFrame | None = None,
    title: str = "Dispersive surface energy profile",
) -> plt.Figure:
    """Plot γ_d as a function of surface coverage.

    Parameters
    ----------
    profile : pd.DataFrame
        Columns: ``coverage``, ``gamma_d_mJm2``.
    reference : pd.DataFrame, optional
        Independent reference values with the same columns (for overlay).
    title : str
        Plot title.

    Returns
    -------
    plt.Figure
    """
    fig, ax = plt.subplots(figsize=(SINGLE_COL, SINGLE_COL * 0.85))
    palette = get_palette()

    ax.plot(profile["coverage"], profile["gamma_d_mJm2"],
            "o-", color=palette[0], markersize=4, linewidth=1,
            label="This pipeline")

    if reference is not None:
        ax.plot(reference["coverage"], reference["gamma_d_mJm2"],
                "s--", color=palette[1], markersize=4, linewidth=0.8,
                label="Reference", alpha=0.7)

    ax.set_xlabel("Surface coverage (n/nm)")
    ax.set_ylabel("γd (mJ/m²)")
    ax.set_title(title, fontsize=9)
    ax.legend(fontsize=7, frameon=False)

    fig.tight_layout()
    return fig


def plot_validation_comparison(
    our_values: pd.DataFrame,
    reference_values: pd.DataFrame,
    title: str = "Pipeline validation",
) -> plt.Figure:
    """Compare a computed profile with an independent reference profile.

    Parameters
    ----------
    our_values : pd.DataFrame
        Columns: ``coverage``, ``gamma_d_mJm2``.
    reference_values : pd.DataFrame
        Columns: ``coverage``, ``gamma_d_mJm2``.
    title : str
        Plot title.

    Returns
    -------
    plt.Figure
    """
    fig, axes = plt.subplots(1, 2, figsize=(DOUBLE_COL, DOUBLE_COL * 0.4))
    palette = get_palette()

    # Merge on coverage
    merged = pd.merge(
        our_values[["coverage", "gamma_d_mJm2"]],
        reference_values[["coverage", "gamma_d_mJm2"]].rename(
            columns={"gamma_d_mJm2": "gamma_d_reference_mJm2"}
        ),
        on="coverage", how="inner",
    )

    # Panel 1: overlay
    ax = axes[0]
    ax.plot(merged["coverage"], merged["gamma_d_mJm2"],
            "o-", color=palette[0], markersize=4, label="Our pipeline")
    ax.plot(merged["coverage"], merged["gamma_d_reference_mJm2"],
            "s--", color=palette[1], markersize=4, label="Reference")
    ax.set_xlabel("Coverage (n/nm)")
    ax.set_ylabel("γd (mJ/m²)")
    ax.set_title("Overlay", fontsize=8)
    ax.legend(fontsize=7, frameon=False)

    # Panel 2: residuals
    ax = axes[1]
    residuals = merged["gamma_d_mJm2"] - merged["gamma_d_reference_mJm2"]
    ax.bar(range(len(merged)), residuals, color=palette[2], alpha=0.7)
    ax.axhline(0, color="k", linewidth=0.5)
    ax.set_xticks(range(len(merged)))
    ax.set_xticklabels([f"{c:.3f}" for c in merged["coverage"]], rotation=45, fontsize=6)
    ax.set_xlabel("Coverage (n/nm)")
    ax.set_ylabel("Δγd (mJ/m²)")
    ax.set_title("Residuals", fontsize=8)

    fig.suptitle(title, fontsize=9, y=1.02)
    fig.tight_layout()
    return fig
