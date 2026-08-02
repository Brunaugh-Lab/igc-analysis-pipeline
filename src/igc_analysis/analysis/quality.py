"""Quality control checks for dispersive surface energy results.

All functions are pure — they accept data and return structured flags
without printing or modifying anything.  Each flag is a dict::

    {"check": str, "severity": str, "coverage": float|None,
     "message": str, "value": Any}

Severity levels: ``"info"``, ``"warning"``, ``"critical"``.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Check 1: Alkane count per coverage
# ---------------------------------------------------------------------------

def check_alkane_count(gamma_d: pd.DataFrame) -> list[dict]:
    """Flag coverages where fewer than 3 alkanes were used.

    With exactly 2 alkanes the linear regression fits a line through
    2 points, giving R² = 1.000 trivially — the result is physically
    meaningless because there is no residual to assess goodness-of-fit.
    """
    flags: list[dict] = []
    if "n_alkanes" not in gamma_d.columns:
        return flags

    for _, row in gamma_d.iterrows():
        n = row["n_alkanes"]
        if pd.isna(n):
            continue
        n = int(n)
        cov = row.get("coverage", None)
        if n < 2:
            flags.append({
                "check": "alkane_count",
                "severity": "critical",
                "coverage": cov,
                "message": f"Only {n} alkane(s) at coverage {cov} — cannot compute gamma_d",
                "value": n,
            })
        elif n == 2:
            flags.append({
                "check": "alkane_count",
                "severity": "warning",
                "coverage": cov,
                "message": (
                    f"Only 2 alkanes at coverage {cov} — R² is trivially 1.000, "
                    f"gamma_d is unreliable"
                ),
                "value": n,
            })
    return flags


# ---------------------------------------------------------------------------
# Check 2: Physical bounds on gamma_d
# ---------------------------------------------------------------------------

def check_gamma_d_bounds(
    gamma_d: pd.DataFrame,
    low: float = 15.0,
    high: float = 80.0,
) -> list[dict]:
    """Flag gamma_d values outside a physically reasonable range.

    Default bounds (15–80 mJ/m²) are appropriate for organic
    pharmaceutical powders.  Adjust for other material classes.

    Values outside [5, 150] are flagged as critical (almost certainly
    an artifact regardless of material).
    """
    flags: list[dict] = []
    for _, row in gamma_d.iterrows():
        gd = row.get("gamma_d_mJm2", np.nan)
        if pd.isna(gd):
            continue
        cov = row.get("coverage", None)

        if gd < 5 or gd > 150:
            flags.append({
                "check": "gamma_d_bounds",
                "severity": "critical",
                "coverage": cov,
                "message": f"gamma_d = {gd:.1f} mJ/m² at coverage {cov} — far outside physical range",
                "value": gd,
            })
        elif gd < low or gd > high:
            flags.append({
                "check": "gamma_d_bounds",
                "severity": "warning",
                "coverage": cov,
                "message": f"gamma_d = {gd:.1f} mJ/m² at coverage {cov} — outside expected [{low}, {high}] range",
                "value": gd,
            })
    return flags


# ---------------------------------------------------------------------------
# Check 3: Profile shape classification
# ---------------------------------------------------------------------------

def classify_profile_shape(gamma_d: pd.DataFrame) -> str:
    """Classify the gamma_d-vs-coverage profile shape.

    Returns one of:

    - ``"normal_decline"`` — gamma_d decreases with coverage (common
      for heterogeneous surfaces where high-energy sites fill first)
    - ``"flat"`` — gamma_d roughly constant (homogeneous surface)
    - ``"u_shaped"`` — drops then rises at high coverage (possible
      multilayer effects or adsorbate-adsorbate interactions)
    - ``"increasing"`` — gamma_d increases with coverage (unusual,
      may indicate column conditioning or measurement artifacts)
    - ``"collapse"`` — gamma_d drops sharply at high coverage
      (instrument failure, gas supply issue, or column overloading)
    - ``"insufficient"`` — fewer than 3 valid data points
    """
    valid = gamma_d.dropna(subset=["gamma_d_mJm2"]).sort_values("coverage")
    vals = valid["gamma_d_mJm2"].values

    if len(vals) < 3:
        return "insufficient"

    max_val = vals.max()
    min_val = vals.min()
    span = max_val - min_val

    # Collapse: last value drops below half of the maximum
    if vals[-1] < 0.5 * max_val:
        return "collapse"

    # Flat: range < 3 mJ/m²
    if span < 3.0:
        return "flat"

    # Check monotonicity with noise tolerance
    diffs = np.diff(vals)
    tol = 0.5  # mJ/m² noise tolerance

    # Monotonically decreasing (all diffs ≤ tol)
    if np.all(diffs <= tol):
        return "normal_decline"

    # Monotonically increasing (all diffs ≥ -tol)
    if np.all(diffs >= -tol):
        return "increasing"

    # U-shaped: minimum is in the interior
    min_idx = np.argmin(vals)
    if 0 < min_idx < len(vals) - 1:
        return "u_shaped"

    # Default: mostly declining with some noise
    return "normal_decline"


def check_profile_shape(gamma_d: pd.DataFrame) -> tuple[str, list[dict]]:
    """Classify profile shape and flag unusual shapes."""
    shape = classify_profile_shape(gamma_d)
    flags: list[dict] = []

    if shape == "collapse":
        flags.append({
            "check": "profile_shape",
            "severity": "critical",
            "coverage": None,
            "message": "Profile collapse — gamma_d drops sharply at high coverage (instrument failure?)",
            "value": shape,
        })
    elif shape == "increasing":
        flags.append({
            "check": "profile_shape",
            "severity": "warning",
            "coverage": None,
            "message": "Increasing profile — gamma_d rises with coverage (unusual, check for artifacts)",
            "value": shape,
        })
    elif shape == "u_shaped":
        flags.append({
            "check": "profile_shape",
            "severity": "warning",
            "coverage": None,
            "message": "U-shaped profile — gamma_d dips then rises (possible multilayer or interaction effects)",
            "value": shape,
        })
    elif shape == "insufficient":
        flags.append({
            "check": "profile_shape",
            "severity": "warning",
            "coverage": None,
            "message": "Fewer than 3 valid gamma_d points — cannot classify profile",
            "value": shape,
        })

    return shape, flags


# ---------------------------------------------------------------------------
# Check 4: V_N ordering by carbon number
# ---------------------------------------------------------------------------

def check_vn_ordering(vn_results: pd.DataFrame) -> list[dict]:
    """Check that V_N increases with carbon number at each coverage.

    In a well-behaved IGC experiment, heavier alkanes have longer
    retention times and therefore larger V_N.  The expected ordering
    is V_N(octane) < V_N(nonane) < V_N(decane), with each roughly
    2–4× the previous.

    Flags:
    - Non-monotonic ordering → critical
    - Consecutive ratio outside [1.5, 6.0] → warning
    """
    flags: list[dict] = []
    alkanes = vn_results[vn_results["carbon_number"].notna()].copy()
    if alkanes.empty:
        return flags

    # Group by target_coverage (which is the grid we interpolate to)
    for cov, group in alkanes.groupby("target_coverage"):
        g = group.sort_values("carbon_number")
        cns = g["carbon_number"].values
        vns = g["VN_cofm"].values

        if len(cns) < 2:
            continue

        # Check monotonicity
        for i in range(1, len(vns)):
            if vns[i] <= vns[i - 1]:
                flags.append({
                    "check": "vn_ordering",
                    "severity": "critical",
                    "coverage": cov,
                    "message": (
                        f"V_N not monotonic at coverage {cov}: "
                        f"C{int(cns[i])}={vns[i]:.1f} <= C{int(cns[i-1])}={vns[i-1]:.1f} mL/g"
                    ),
                    "value": {"cn": int(cns[i]), "vn": vns[i], "cn_prev": int(cns[i-1]), "vn_prev": vns[i-1]},
                })

            # Check ratio
            if vns[i - 1] > 0:
                ratio = vns[i] / vns[i - 1]
                if ratio < 1.5 or ratio > 6.0:
                    flags.append({
                        "check": "vn_ordering",
                        "severity": "warning",
                        "coverage": cov,
                        "message": (
                            f"Unusual V_N ratio at coverage {cov}: "
                            f"C{int(cns[i])}/C{int(cns[i-1])} = {ratio:.2f} "
                            f"(expected 1.5–6.0)"
                        ),
                        "value": {"ratio": ratio, "cn": int(cns[i]), "cn_prev": int(cns[i-1])},
                    })
    return flags


# ---------------------------------------------------------------------------
# Check 5: CoM vs Peak Max divergence
# ---------------------------------------------------------------------------

def check_com_divergence(
    gamma_d: pd.DataFrame,
    threshold: float = 3.0,
) -> list[dict]:
    """Flag large divergence between CoM (primary) and Peak Max gamma_d.

    Large divergence indicates severe peak tailing.  When the two
    methods disagree by more than *threshold* mJ/m², the gamma_d
    value is sensitive to the retention time definition and cross-
    formulation comparisons should be interpreted with caution.
    """
    flags: list[dict] = []
    if "delta_cofm_pm" not in gamma_d.columns:
        return flags

    deltas = gamma_d["delta_cofm_pm"].dropna()
    if deltas.empty:
        return flags

    max_abs = deltas.abs().max()
    if max_abs > threshold:
        flags.append({
            "check": "com_divergence",
            "severity": "warning",
            "coverage": None,
            "message": (
                f"CoM vs Peak Max gamma_d divergence up to {max_abs:.1f} mJ/m\u00b2 "
                f"(threshold {threshold:.0f}) \u2014 severe peak tailing; "
                f"gamma_d ranking may depend on method choice"
            ),
            "value": max_abs,
        })
    return flags


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------

def run_qc_checks(
    gamma_d: pd.DataFrame,
    vn_results: pd.DataFrame,
    bounds: tuple[float, float] = (15.0, 80.0),
) -> dict:
    """Run all quality checks and return a structured result.

    Parameters
    ----------
    gamma_d : pd.DataFrame
        Output of :func:`dorris_gray_gamma_d`.
    vn_results : pd.DataFrame
        Per-injection V_N results from the pipeline.
    bounds : tuple
        ``(low, high)`` gamma_d bounds in mJ/m² for the domain.

    Returns
    -------
    dict
        Keys:
        - ``flags`` : list[dict] — all individual flags
        - ``profile_shape`` : str — classification label
        - ``summary`` : str — one-line summary
        - ``pass`` : bool — True if zero critical flags
    """
    all_flags: list[dict] = []

    # Check 1: alkane count
    all_flags.extend(check_alkane_count(gamma_d))

    # Check 2: physical bounds
    all_flags.extend(check_gamma_d_bounds(gamma_d, low=bounds[0], high=bounds[1]))

    # Check 3: profile shape
    shape, shape_flags = check_profile_shape(gamma_d)
    all_flags.extend(shape_flags)

    # Check 4: V_N ordering
    all_flags.extend(check_vn_ordering(vn_results))

    # Check 5: Peak Max vs CoM divergence
    all_flags.extend(check_com_divergence(gamma_d))

    n_warn = sum(1 for f in all_flags if f["severity"] == "warning")
    n_crit = sum(1 for f in all_flags if f["severity"] == "critical")
    passed = n_crit == 0

    status = "PASS" if passed else "FAIL"
    summary = f"QC {status}: {n_warn} warning(s), {n_crit} critical; profile: {shape}"

    return {
        "flags": all_flags,
        "profile_shape": shape,
        "summary": summary,
        "pass": passed,
    }
