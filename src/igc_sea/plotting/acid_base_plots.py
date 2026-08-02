"""Visualization for acid-base surface chemistry analysis.

Planned plots:
  - Gutmann plot (ΔG_sp/AN* vs DN/AN*)
  - ΔG_sp bar chart by probe
  - K_a / K_b comparison across samples
"""

import matplotlib.pyplot as plt
import pandas as pd


def plot_gutmann(df: pd.DataFrame) -> plt.Figure:
    """Plot the Gutmann regression (ΔG_sp/AN* vs DN/AN*).

    Parameters
    ----------
    df : pd.DataFrame
        Columns: sample_name, probe, dn_over_an_star, delta_g_over_an_star

    Returns
    -------
    matplotlib.figure.Figure

    Raises
    ------
    NotImplementedError
    """
    raise NotImplementedError("Requires acid-base analysis data.")


def plot_delta_g_bars(df: pd.DataFrame) -> plt.Figure:
    """Bar chart of ΔG_sp per polar probe.

    Parameters
    ----------
    df : pd.DataFrame
        Columns: sample_name, probe, delta_G_sp_kJmol

    Returns
    -------
    matplotlib.figure.Figure

    Raises
    ------
    NotImplementedError
    """
    raise NotImplementedError("Requires acid-base analysis data.")
