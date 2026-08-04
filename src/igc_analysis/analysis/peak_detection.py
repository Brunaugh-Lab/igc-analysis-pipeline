"""Peak detection and integration for inverse chromatography traces.

Processes raw FID signal traces to extract:
- Baseline (linear fit to signal endpoints)
- Peak maximum time
- Center-of-mass (first moment) retention time
- Integrated peak area
- Peak width
"""

from __future__ import annotations

import numpy as np

# np.trapz was removed in NumPy 2.0 in favour of np.trapezoid. Bind whichever
# exists so the pipeline runs on both NumPy 1.x and 2.x.
_trapezoid = getattr(np, "trapezoid", None) or np.trapz


def detect_baseline(
    time: np.ndarray,
    signal: np.ndarray,
    start_frac: float = 0.05,
    end_frac: float = 0.05,
    peak_aware: bool = True,
    dead_time_min: float = 0.0,
) -> tuple[float, float]:
    """Estimate a linear baseline from regions flanking the peak.

    When *peak_aware* is True (default), the baseline is fit through
    the quiet regions immediately before and after the chromatographic
    peak. This avoids
    contamination from the injection transient at t≈0 and from
    far-tail drift.

    Parameters
    ----------
    time : np.ndarray
        Time axis (minutes).
    signal : np.ndarray
        FID signal (µV).
    start_frac : float
        Legacy: fraction of total points to use from the start.
        Only used when *peak_aware* is False.
    end_frac : float
        Legacy: fraction of total points from the end.
    peak_aware : bool
        If True, find the peak and fit baseline from the pre-peak
        and post-peak flat regions.
    dead_time_min : float
        Approximate dead time (min).  Points before this are skipped
        when selecting the pre-peak baseline region.

    Returns
    -------
    tuple[float, float]
        ``(intercept, gradient)`` of the baseline line:
        ``baseline(t) = intercept + gradient * t``.
    """
    if not peak_aware:
        return _baseline_from_edges(time, signal, start_frac, end_frac)

    return _baseline_peak_aware(time, signal, dead_time_min)


def _baseline_from_edges(
    time: np.ndarray,
    signal: np.ndarray,
    start_frac: float,
    end_frac: float,
) -> tuple[float, float]:
    """Legacy baseline: fit from first/last fraction of trace."""
    n = len(time)
    n_start = max(3, int(n * start_frac))
    n_end = max(3, int(n * end_frac))

    t_bl = np.concatenate([time[:n_start], time[-n_end:]])
    s_bl = np.concatenate([signal[:n_start], signal[-n_end:]])

    coeffs = np.polyfit(t_bl, s_bl, 1)
    return float(coeffs[1]), float(coeffs[0])


def _baseline_peak_aware(
    time: np.ndarray,
    signal: np.ndarray,
    dead_time_min: float,
) -> tuple[float, float]:
    """Two-point linear baseline anchored in the pre-peak and post-peak
    regions of the chromatogram.

    A detector trace has three distinct regions:

    1. **Injection transient** (t ≈ 0 to ~0.2 min) — pressure/flow
       disturbance from the solvent pulse injection.
    2. **Analyte peak** — the probe molecule eluting through the column,
       detected by the FID.
    3. **Post-peak tail & drift** — the signal returns toward baseline
       but can drift over long traces due to FID temperature changes.

    The baseline is the FID response in the absence of analyte.  We
    estimate it as a linear function connecting two anchor points:

    - **Pre-peak anchor** at ~48% of the dead time — after the injection
      transient has settled but before the analyte peak begins.  The
      carrier gas (no analyte) is passing through the FID here.
    - **Post-peak anchor** at ~60% of the total trace time — after the
      bulk of the analyte has eluted.  Some residual tail may remain,
      but this point represents the best compromise between capturing
      the full peak and avoiding long-term FID drift.

    Both anchors use the median signal in a ±0.02 min window for noise
    robustness.
    """
    n = len(time)

    # Find peak index
    peak_idx = np.argmax(signal)
    peak_time = time[peak_idx]

    # --- Pre-peak anchor ---
    # At ~48% of dead time: the injection transient (valve actuation,
    # pressure pulse) has dissipated and only carrier gas flows through
    # the FID.  The analyte peak has not yet begun to elute.
    if dead_time_min > 0:
        pre_anchor_time = dead_time_min * 0.48
    else:
        pre_anchor_time = peak_time * 0.25

    pre_mask = (time >= pre_anchor_time - 0.02) & (time <= pre_anchor_time + 0.02)
    if pre_mask.sum() < 3:
        pre_mask = (time >= pre_anchor_time - 0.05) & (time <= pre_anchor_time + 0.05)
    if pre_mask.sum() < 2:
        return _baseline_from_edges(time, signal, 0.05, 0.05)

    pre_signal_val = np.median(signal[pre_mask])
    pre_time_val = np.median(time[pre_mask])

    # --- Post-peak anchor ---
    # At ~60% of total trace time: the analyte peak has largely eluted.
    # For short-retention probes (octane, stage_time=4 min) this lands
    # at ~2.4 min; for long-retention probes (decane, stage_time=20 min)
    # at ~12 min.  Placing it here rather than at the trace end avoids
    # long-term FID baseline drift that accumulates over the full
    # recording window.
    trace_duration = time[-1] - time[0]
    post_anchor_time = time[0] + trace_duration * 0.60

    post_mask = (time >= post_anchor_time - 0.02) & (time <= post_anchor_time + 0.02)
    if post_mask.sum() < 3:
        post_mask = (time >= post_anchor_time - 0.05) & (time <= post_anchor_time + 0.05)
    if post_mask.sum() < 2:
        return _baseline_from_edges(time, signal, 0.05, 0.05)

    post_signal_val = np.median(signal[post_mask])
    post_time_val = np.median(time[post_mask])

    # --- Two-point baseline ---
    if abs(post_time_val - pre_time_val) < 1e-10:
        return float(pre_signal_val), 0.0

    gradient = (post_signal_val - pre_signal_val) / (post_time_val - pre_time_val)
    intercept = pre_signal_val - gradient * pre_time_val
    return float(intercept), float(gradient)


def subtract_baseline(
    time: np.ndarray,
    signal: np.ndarray,
    intercept: float,
    gradient: float,
) -> np.ndarray:
    """Subtract a linear baseline from the signal.

    Parameters
    ----------
    time, signal : np.ndarray
        Chromatogram data.
    intercept, gradient : float
        Baseline parameters.

    Returns
    -------
    np.ndarray
        Baseline-corrected signal.
    """
    baseline = intercept + gradient * time
    return signal - baseline


def find_peak_max(
    time: np.ndarray,
    corrected_signal: np.ndarray,
) -> float:
    """Return the time at which the baseline-corrected signal is maximum.

    Parameters
    ----------
    time : np.ndarray
        Time axis (minutes).
    corrected_signal : np.ndarray
        Baseline-subtracted FID signal.

    Returns
    -------
    float
        Time of peak maximum (minutes).
    """
    idx = np.argmax(corrected_signal)
    return float(time[idx])


def find_peak_cofm(
    time: np.ndarray,
    corrected_signal: np.ndarray,
) -> float:
    """Compute the center-of-mass (first moment) retention time.

    .. math::

        t_{CoM} = \\frac{\\sum_i s_i \\cdot t_i}{\\sum_i s_i}

    Only positive signal values contribute (noise below baseline excluded).

    Parameters
    ----------
    time : np.ndarray
        Time axis (minutes).
    corrected_signal : np.ndarray
        Baseline-subtracted FID signal.

    Returns
    -------
    float
        Center-of-mass retention time (minutes).
    """
    # Only use positive signal (above baseline). Integrate with the time axis so
    # nonuniform but contract-valid sampling does not bias the first moment.
    mask = corrected_signal > 0
    if not mask.any():
        return find_peak_max(time, corrected_signal)

    signal_positive = np.maximum(corrected_signal, 0.0)
    denominator = float(_trapezoid(signal_positive, time))
    if denominator <= 0:
        return find_peak_max(time, corrected_signal)
    numerator = float(_trapezoid(signal_positive * time, time))
    return numerator / denominator


def integrate_peak(
    time: np.ndarray,
    corrected_signal: np.ndarray,
) -> float:
    """Integrate the baseline-corrected peak using the trapezoidal rule.

    Parameters
    ----------
    time : np.ndarray
        Time axis (minutes).
    corrected_signal : np.ndarray
        Baseline-subtracted FID signal.

    Returns
    -------
    float
        Peak area (µV·min).
    """
    # Only integrate positive signal
    s_clipped = np.maximum(corrected_signal, 0.0)
    return float(_trapezoid(s_clipped, time))


def peak_width(
    time: np.ndarray,
    corrected_signal: np.ndarray,
    fraction: float = 0.5,
) -> float:
    """Estimate peak width at a given fraction of peak height.

    Parameters
    ----------
    time : np.ndarray
        Time axis (minutes).
    corrected_signal : np.ndarray
        Baseline-subtracted signal.
    fraction : float
        Fraction of peak height (default 0.5 for FWHM).

    Returns
    -------
    float
        Width in minutes.
    """
    peak_height = np.max(corrected_signal)
    threshold = fraction * peak_height

    above = corrected_signal >= threshold
    if not above.any():
        return 0.0

    indices = np.where(above)[0]
    return float(time[indices[-1]] - time[indices[0]])


def _half_widths_at_fraction(
    time: np.ndarray,
    corrected_signal: np.ndarray,
    fraction: float,
) -> tuple[float, float]:
    """Return (leading_half_width, trailing_half_width) at a given
    fraction of peak height.

    Uses linear interpolation to find the exact crossing times on
    each side of the peak maximum.
    """
    peak_idx = np.argmax(corrected_signal)
    peak_height = corrected_signal[peak_idx]
    threshold = fraction * peak_height

    if peak_height <= 0 or threshold <= 0:
        return 0.0, 0.0

    t_peak = time[peak_idx]

    # Leading edge: walk backward from peak to find threshold crossing
    t_lead = t_peak  # default if not found
    for i in range(peak_idx, 0, -1):
        if corrected_signal[i - 1] < threshold <= corrected_signal[i]:
            # Linear interpolation between points i-1 and i
            frac_interp = ((threshold - corrected_signal[i - 1])
                           / (corrected_signal[i] - corrected_signal[i - 1]))
            t_lead = time[i - 1] + frac_interp * (time[i] - time[i - 1])
            break

    # Trailing edge: walk forward from peak to find threshold crossing
    t_trail = t_peak  # default if not found
    for i in range(peak_idx, len(corrected_signal) - 1):
        if corrected_signal[i] >= threshold > corrected_signal[i + 1]:
            frac_interp = ((corrected_signal[i] - threshold)
                           / (corrected_signal[i] - corrected_signal[i + 1]))
            t_trail = time[i] + frac_interp * (time[i + 1] - time[i])
            break

    leading = t_peak - t_lead
    trailing = t_trail - t_peak
    return max(leading, 0.0), max(trailing, 0.0)


def asymmetry_factor(
    time: np.ndarray,
    corrected_signal: np.ndarray,
) -> float:
    """Compute the peak asymmetry factor at 10% of peak height.

    .. math::

        A_s = \\frac{b}{a}

    where *a* is the leading half-width and *b* is the trailing
    half-width, both measured at 10% of the peak height from the
    baseline.

    Returns 1.0 for a perfectly symmetric peak.  Values > 1.0
    indicate tailing; values > 2.0 indicate severe tailing.
    """
    a, b = _half_widths_at_fraction(time, corrected_signal, 0.10)
    if a <= 0:
        return float("nan")
    return b / a


def tailing_factor(
    time: np.ndarray,
    corrected_signal: np.ndarray,
) -> float:
    """Compute the USP tailing factor at 5% of peak height.

    .. math::

        T_f = \\frac{a + b}{2a}

    where *a* is the leading half-width and *b* is the trailing
    half-width, both measured at 5% of peak height.

    Returns 1.0 for a perfectly symmetric peak.
    """
    a, b = _half_widths_at_fraction(time, corrected_signal, 0.05)
    if a <= 0:
        return float("nan")
    return (a + b) / (2 * a)


def process_chromatogram(
    time: np.ndarray,
    signal: np.ndarray,
    baseline_start_frac: float = 0.05,
    baseline_end_frac: float = 0.05,
    dead_time_min: float = 0.0,
) -> dict:
    """Full peak detection pipeline for a single chromatogram.

    Parameters
    ----------
    time : np.ndarray
        Time axis (minutes).
    signal : np.ndarray
        Raw FID signal (µV).
    baseline_start_frac, baseline_end_frac : float
        Fractions of the trace used for baseline estimation (legacy mode).
    dead_time_min : float
        Approximate dead time (minutes).  Used by the peak-aware baseline
        to skip the injection transient.

    Returns
    -------
    dict
        Keys: ``baseline_intercept``, ``baseline_gradient``,
        ``peak_max_time``, ``peak_cofm``, ``peak_area``,
        ``peak_max_value``, ``peak_width_half``.
    """
    intercept, gradient = detect_baseline(
        time, signal, baseline_start_frac, baseline_end_frac,
        dead_time_min=dead_time_min,
    )
    corrected = subtract_baseline(time, signal, intercept, gradient)

    t_max = find_peak_max(time, corrected)
    t_cofm = find_peak_cofm(time, corrected)
    w_half = peak_width(time, corrected, 0.5)

    return {
        "baseline_intercept": intercept,
        "baseline_gradient": gradient,
        "peak_max_time": t_max,
        "peak_cofm": t_cofm,
        "peak_area": integrate_peak(time, corrected),
        "peak_max_value": float(np.max(corrected)),
        "peak_width_half": w_half,
        "asymmetry_factor": asymmetry_factor(time, corrected),
        "tailing_factor": tailing_factor(time, corrected),
        "com_max_divergence_min": t_cofm - t_max,
        "com_max_divergence_frac": (t_cofm - t_max) / w_half if w_half > 0 else 0.0,
    }
