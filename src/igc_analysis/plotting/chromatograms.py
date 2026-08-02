"""Raw detector-trace visualization for inverse gas chromatography data."""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from igc_analysis.plotting.theme import get_palette, SINGLE_COL, DOUBLE_COL


def plot_chromatogram(
    time: np.ndarray,
    signal: np.ndarray,
    peak_params: dict | None = None,
    title: str = "",
    ax: plt.Axes | None = None,
) -> plt.Figure:
    """Plot a single FID chromatogram with optional peak annotations.

    Parameters
    ----------
    time : np.ndarray
        Time axis (minutes).
    signal : np.ndarray
        FID signal (µV).
    peak_params : dict, optional
        Output of ``peak_detection.process_chromatogram()``.
        If provided, marks the baseline, peak max, and center of mass.
    title : str
        Plot title.
    ax : plt.Axes, optional
        Existing axes to plot on.

    Returns
    -------
    plt.Figure
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(DOUBLE_COL, DOUBLE_COL * 0.45))
    else:
        fig = ax.figure

    palette = get_palette()

    ax.plot(time, signal, color=palette[0], linewidth=0.6, label="FID signal")

    if peak_params is not None:
        # Baseline
        bl_intercept = peak_params["baseline_intercept"]
        bl_gradient = peak_params["baseline_gradient"]
        bl_line = bl_intercept + bl_gradient * time
        ax.plot(time, bl_line, "--", color=palette[7], linewidth=0.5,
                label="Baseline")

        # Peak max
        t_max = peak_params["peak_max_time"]
        ax.axvline(t_max, color=palette[5], linewidth=0.5, linestyle=":",
                    label=f"Peak max ({t_max:.4f} min)")

        # Center of mass
        t_cofm = peak_params["peak_cofm"]
        ax.axvline(t_cofm, color=palette[4], linewidth=0.5, linestyle="-.",
                    label=f"CoM ({t_cofm:.4f} min)")

    ax.set_xlabel("Time (min)")
    ax.set_ylabel("FID signal (µV)")
    if title:
        ax.set_title(title, fontsize=9)
    ax.legend(fontsize=7, frameon=False)

    fig.tight_layout()
    return fig


def plot_chromatogram_overlay(
    chromatograms: dict[str, tuple[np.ndarray, np.ndarray]],
    labels: dict[str, str] | None = None,
    normalize: bool = False,
    title: str = "",
) -> plt.Figure:
    """Overlay multiple chromatograms for comparison.

    Parameters
    ----------
    chromatograms : dict
        Keys are identifiers (e.g. ``"injection1"``), values are
        ``(time, signal)`` tuples.
    labels : dict, optional
        Human-readable labels for each key. If None, uses the keys.
    normalize : bool
        If True, normalize each trace to unit peak height.
    title : str
        Plot title.

    Returns
    -------
    plt.Figure
    """
    fig, ax = plt.subplots(figsize=(DOUBLE_COL, DOUBLE_COL * 0.5))
    palette = get_palette(len(chromatograms))

    for i, (key, (time, signal)) in enumerate(chromatograms.items()):
        label = (labels or {}).get(key, key)
        y = signal
        if normalize and np.max(signal) > 0:
            y = signal / np.max(signal)

        ax.plot(time, y, color=palette[i % len(palette)],
                linewidth=0.6, label=label)

    ax.set_xlabel("Time (min)")
    ax.set_ylabel("Normalized signal" if normalize else "FID signal (µV)")
    if title:
        ax.set_title(title, fontsize=9)
    ax.legend(fontsize=7, frameon=False, loc="upper right")

    fig.tight_layout()
    return fig
