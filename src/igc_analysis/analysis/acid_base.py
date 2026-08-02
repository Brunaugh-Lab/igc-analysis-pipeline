"""Acid-base surface chemistry analysis: Gutmann and van Oss approaches.

Physics chain (Gutmann):
    alkane V_N → Schultz reference line (RT·ln(V_N) vs a·√γ_L^d)
    polar probe V_N → RT·ln(V_N) → subtract alkane line → ΔG_sp
    ΔG_sp/AN* vs DN/AN* → linear regression → K_a (slope), K_b (intercept)

Physics chain (van Oss):
    ΔG_sp from monopolar acid probe (DCM) → γ_S⁻ (surface base component)
    ΔG_sp from monopolar base probe (EtAc) → γ_S⁺ (surface acid component)
    γ_S^AB = 2·√(γ_S⁺ · γ_S⁻)

The Schultz approach places all probes (alkane and polar) on a common axis
using the probe's cross-sectional area and dispersive surface tension.  The
alkane reference line captures pure dispersive interactions; a polar probe's
deviation above this line is the specific (acid-base) free energy ΔG_sp.

The Gutmann equation decomposes ΔG_sp into acid and base contributions:
    ΔG_sp = K_a · DN + K_b · AN*
where DN (donor number) and AN* (modified acceptor number) are tabulated
probe properties.  Rearranging for regression:
    ΔG_sp / AN* = K_a · (DN / AN*) + K_b

The van Oss approach uses monopolar probes to directly solve for the acid
(γ_S⁺) and base (γ_S⁻) components of surface energy.  Unlike Gutmann,
this gives results in mJ/m² units and requires only two probes (one acid,
one base).

References:
    - Schultz, Lavielle & Martin (1987) J. Adhesion 23, 45-60
    - Gutmann (1978) The Donor-Acceptor Approach to Molecular Interactions
    - Riddle & Fowkes (1990) — AN* correction removing dispersive contribution
    - van Oss, Good & Chaudhury (1988) Langmuir 4, 884-891
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

from igc_analysis.constants import R_GAS, N_AVOGADRO
from igc_analysis.utils import PROBE_MOLECULES, get_probe


# ---------------------------------------------------------------------------
# Solvent name → probe database lookup
# ---------------------------------------------------------------------------

# Map instrument solvent names (uppercase) to probe database names (lowercase).
# Only probes that appear in our experiments and have DN/AN* values.
_SOLVENT_TO_PROBE = {
    "DICHLOROMETHANE": "dichloromethane",
    "ETHYL ACETATE": "ethyl acetate",
    "ETHANOL": "ethanol",
    "ACETONE": "acetone",
    "ACETONITRILE": "acetonitrile",
    "CHLOROFORM": "chloroform",
    "TETRAHYDROFURAN": "tetrahydrofuran",
    "DIETHYL ETHER": "diethyl ether",
}


def _resolve_probe(solvent_name: str) -> pd.Series | None:
    """Look up probe properties for an instrument solvent name.

    Returns None if the solvent is not a recognized polar probe.
    """
    probe_name = _SOLVENT_TO_PROBE.get(solvent_name.upper())
    if probe_name is None:
        # Try direct case-insensitive lookup
        try:
            return get_probe(solvent_name)
        except KeyError:
            return None
    return get_probe(probe_name)


# ---------------------------------------------------------------------------
# Schultz reference line (alkane)
# ---------------------------------------------------------------------------

def schultz_parameter(a_cross: float, gamma_l_d_mJm2: float) -> float:
    """Compute the Schultz x-axis parameter: a · √(γ_L^d).

    Parameters
    ----------
    a_cross : float
        Molecular cross-sectional area (m²).
    gamma_l_d_mJm2 : float
        Dispersive surface tension of the liquid probe (mJ/m²).

    Returns
    -------
    float
        Schultz parameter in m² · (J/m²)^0.5.
    """
    return a_cross * np.sqrt(gamma_l_d_mJm2 * 1e-3)


def schultz_reference_line(
    alkane_vn: pd.DataFrame,
    temperature: float,
    vn_column: str = "VN",
    solvent_column: str = "solvent_name",
) -> dict:
    """Fit the Schultz alkane reference line at a single coverage.

    Plots RT·ln(V_N) vs a·√(γ_L^d) for the alkane series and fits a
    linear regression.  The slope encodes √(γ_S^d) and the line serves
    as the reference for computing polar probe ΔG_sp.

    Parameters
    ----------
    alkane_vn : pd.DataFrame
        Alkane data at one coverage.  Must contain ``solvent_column`` and
        ``vn_column``.  Solvent names should be uppercase instrument names
        (e.g. "OCTANE", "NONANE", "DECANE").
    temperature : float
        Column temperature (K).
    vn_column : str
        Column containing specific net retention volume (mL/g).
    solvent_column : str
        Column containing solvent name.

    Returns
    -------
    dict
        Keys: ``slope``, ``intercept``, ``r_squared``, ``n_alkanes``,
        ``gamma_d_schultz_mJm2`` (dispersive γ from Schultz slope),
        ``alkane_points`` (list of dicts with probe details).
    """
    alkane_names = {"HEXANE", "HEPTANE", "OCTANE", "NONANE", "DECANE"}

    points = []
    for _, row in alkane_vn.iterrows():
        name = row[solvent_column].upper()
        if name not in alkane_names:
            continue
        vn = row[vn_column]
        if pd.isna(vn) or vn <= 0:
            continue
        probe = get_probe(name.lower())
        x = schultz_parameter(probe["a_cross"], probe["gamma_l_d"])
        y = R_GAS * temperature * np.log(vn)
        points.append({
            "solvent": name,
            "x_schultz": x,
            "RT_ln_VN": y,
            "VN": vn,
        })

    if len(points) < 2:
        return {
            "slope": np.nan,
            "intercept": np.nan,
            "r_squared": np.nan,
            "n_alkanes": len(points),
            "gamma_d_schultz_mJm2": np.nan,
            "alkane_points": points,
        }

    x_arr = np.array([p["x_schultz"] for p in points])
    y_arr = np.array([p["RT_ln_VN"] for p in points])

    slope, intercept, r_value, p_value, std_err = stats.linregress(x_arr, y_arr)

    # γ_d from Schultz slope: slope = 2 · N_A · √(γ_S^d)
    # → γ_S^d = (slope / (2 · N_A))²
    gamma_d = (slope / (2.0 * N_AVOGADRO)) ** 2
    gamma_d_mJm2 = gamma_d * 1e3  # J/m² → mJ/m²

    return {
        "slope": slope,
        "intercept": intercept,
        "r_squared": r_value ** 2,
        "n_alkanes": len(points),
        "gamma_d_schultz_mJm2": gamma_d_mJm2,
        "alkane_points": points,
    }


# ---------------------------------------------------------------------------
# ΔG_sp calculation
# ---------------------------------------------------------------------------

def calculate_delta_g_sp(
    polar_vn: pd.DataFrame,
    schultz_line: dict,
    temperature: float,
    vn_column: str = "VN",
    solvent_column: str = "solvent_name",
) -> pd.DataFrame:
    """Calculate specific free energy of adsorption ΔG_sp for polar probes.

    ΔG_sp is the vertical distance between a polar probe's RT·ln(V_N) and
    the alkane Schultz reference line at the same x = a·√(γ_L^d).

    Parameters
    ----------
    polar_vn : pd.DataFrame
        Polar probe data at one coverage.  Must contain ``solvent_column``
        and ``vn_column``.
    schultz_line : dict
        Output of :func:`schultz_reference_line`.
    temperature : float
        Column temperature (K).

    Returns
    -------
    pd.DataFrame
        Columns: ``probe``, ``solvent_name``, ``VN``, ``RT_ln_VN``,
        ``x_schultz``, ``alkane_predicted``, ``delta_g_sp_Jmol``,
        ``delta_g_sp_kJmol``, ``dn``, ``an_star``, ``category``.
    """
    slope = schultz_line["slope"]
    intercept = schultz_line["intercept"]

    if np.isnan(slope):
        return pd.DataFrame()

    rows = []
    for _, row in polar_vn.iterrows():
        solvent = row[solvent_column]
        probe = _resolve_probe(solvent)
        if probe is None:
            continue

        # Skip alkanes
        if probe["category"] == "alkane":
            continue

        vn = row[vn_column]
        if pd.isna(vn) or vn <= 0:
            continue

        x = schultz_parameter(probe["a_cross"], probe["gamma_l_d"])
        rt_ln_vn = R_GAS * temperature * np.log(vn)
        alkane_predicted = slope * x + intercept
        delta_g_sp = rt_ln_vn - alkane_predicted  # J/mol

        rows.append({
            "probe": probe["name"],
            "solvent_name": solvent,
            "VN": vn,
            "RT_ln_VN": rt_ln_vn,
            "x_schultz": x,
            "alkane_predicted": alkane_predicted,
            "delta_g_sp_Jmol": delta_g_sp,
            "delta_g_sp_kJmol": delta_g_sp / 1000.0,
            "dn": probe["dn"],
            "an_star": probe["an_star"],
            "category": probe["category"],
        })

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Gutmann K_a / K_b regression
# ---------------------------------------------------------------------------

def gutmann_ka_kb(
    delta_g_sp_df: pd.DataFrame,
    exclude_probes: list[str] | None = None,
) -> dict:
    """Calculate K_a and K_b from Gutmann plot regression.

    Regresses ΔG_sp/AN* vs DN/AN* for polar probes at one coverage.

    Parameters
    ----------
    delta_g_sp_df : pd.DataFrame
        Output of :func:`calculate_delta_g_sp`.  Must contain columns
        ``delta_g_sp_kJmol``, ``dn``, ``an_star``, ``probe``.
    exclude_probes : list of str, optional
        Probe names to exclude from regression (e.g. ``["acetonitrile"]``).

    Returns
    -------
    dict
        Keys: ``Ka``, ``Kb``, ``Kb_Ka_ratio``, ``r_squared``, ``n_probes``,
        ``probes_used``, ``residuals`` (list of per-probe residual dicts).
    """
    if delta_g_sp_df.empty:
        return _empty_gutmann_result()

    df = delta_g_sp_df.copy()

    # Exclude specified probes
    if exclude_probes:
        df = df[~df["probe"].str.lower().isin([p.lower() for p in exclude_probes])]

    # Exclude probes with AN* == 0 (can't divide) or negative ΔG_sp
    df = df[df["an_star"] > 0].copy()

    if len(df) < 2:
        return _empty_gutmann_result(n_probes=len(df))

    # Gutmann regression: ΔG_sp/AN* = Ka · (DN/AN*) + Kb
    df["y"] = df["delta_g_sp_kJmol"] / df["an_star"]
    df["x"] = df["dn"] / df["an_star"]

    slope, intercept, r_value, p_value, std_err = stats.linregress(
        df["x"].values, df["y"].values,
    )

    ka = slope
    kb = intercept

    # 2 probes: deterministic solve, R² is meaningless (always 1.0)
    if len(df) == 2:
        fit_method = "deterministic"
        r_squared = np.nan
    else:
        fit_method = "regression"
        r_squared = r_value ** 2

    # Per-probe residuals
    residuals = []
    for _, row in df.iterrows():
        predicted = ka * row["x"] + kb
        residual = row["y"] - predicted
        residuals.append({
            "probe": row["probe"],
            "delta_g_sp_kJmol": row["delta_g_sp_kJmol"],
            "dn": row["dn"],
            "an_star": row["an_star"],
            "x_gutmann": row["x"],
            "y_gutmann": row["y"],
            "predicted": predicted,
            "residual": residual,
        })

    return {
        "Ka": ka,
        "Kb": kb,
        "Kb_Ka_ratio": kb / ka if ka != 0 else np.inf,
        "r_squared": r_squared,
        "fit_method": fit_method,
        "n_probes": len(df),
        "probes_used": list(df["probe"]),
        "residuals": residuals,
    }


def _empty_gutmann_result(n_probes: int = 0) -> dict:
    return {
        "Ka": np.nan,
        "Kb": np.nan,
        "Kb_Ka_ratio": np.nan,
        "r_squared": np.nan,
        "fit_method": "insufficient",
        "n_probes": n_probes,
        "probes_used": [],
        "residuals": [],
    }


# ---------------------------------------------------------------------------
# Leave-one-out outlier detection
# ---------------------------------------------------------------------------

def leave_one_out_influence(delta_g_sp_df: pd.DataFrame) -> list[dict]:
    """Assess each probe's influence on Ka/Kb via leave-one-out.

    For each probe, recompute Ka/Kb without it and report the change.
    Large changes indicate an influential (potentially outlier) probe.

    Parameters
    ----------
    delta_g_sp_df : pd.DataFrame
        Output of :func:`calculate_delta_g_sp`.

    Returns
    -------
    list of dict
        Per-probe influence: ``probe``, ``delta_Ka``, ``delta_Kb``,
        ``delta_r_squared``, ``is_outlier``.
    """
    if len(delta_g_sp_df) < 3:
        return []

    full = gutmann_ka_kb(delta_g_sp_df)
    if np.isnan(full["Ka"]):
        return []

    results = []
    for probe_name in delta_g_sp_df["probe"].unique():
        reduced = gutmann_ka_kb(delta_g_sp_df, exclude_probes=[probe_name])
        if np.isnan(reduced["Ka"]):
            continue

        delta_ka = reduced["Ka"] - full["Ka"]
        delta_kb = reduced["Kb"] - full["Kb"]
        delta_r2 = reduced["r_squared"] - full["r_squared"]

        # Outlier heuristic: removing it improves R² by > 0.15 or
        # changes Ka or Kb by > 50%
        ka_pct = abs(delta_ka / full["Ka"]) if full["Ka"] != 0 else 0
        kb_pct = abs(delta_kb / full["Kb"]) if full["Kb"] != 0 else 0
        is_outlier = delta_r2 > 0.15 or ka_pct > 0.5 or kb_pct > 0.5

        results.append({
            "probe": probe_name,
            "delta_Ka": delta_ka,
            "delta_Kb": delta_kb,
            "delta_r_squared": delta_r2,
            "is_outlier": is_outlier,
        })

    return results


# ---------------------------------------------------------------------------
# Quality checks
# ---------------------------------------------------------------------------

def acid_base_quality_checks(
    gutmann_result: dict,
    delta_g_sp_df: pd.DataFrame,
    influence: list[dict] | None = None,
    van_oss_result: dict | None = None,
) -> dict:
    """Run QC checks on acid-base analysis results.

    Parameters
    ----------
    gutmann_result : dict
        Output of :func:`gutmann_ka_kb`.
    delta_g_sp_df : pd.DataFrame
        Output of :func:`calculate_delta_g_sp`.
    influence : list of dict, optional
        Output of :func:`leave_one_out_influence`.
    van_oss_result : dict, optional
        Output of :func:`compute_van_oss`.

    Returns
    -------
    dict
        Keys: ``flags`` (list of flag dicts), ``summary`` (str).
    """
    flags = []

    fit_method = gutmann_result.get("fit_method", "regression")

    # R² checks — only meaningful for regression (3+ probes)
    if fit_method == "regression":
        r2 = gutmann_result.get("r_squared", np.nan)
        if not np.isnan(r2):
            if r2 < 0.3:
                flags.append({
                    "code": "GUTMANN_R2_CRITICAL",
                    "severity": "critical",
                    "message": f"Gutmann R² = {r2:.3f} < 0.3 — model is unreliable",
                })
            elif r2 < 0.5:
                flags.append({
                    "code": "GUTMANN_R2_LOW",
                    "severity": "warning",
                    "message": f"Gutmann R² = {r2:.3f} < 0.5 — interpret Ka/Kb with caution",
                })

    # Deterministic fit: 2 probes, no degrees of freedom
    if fit_method == "deterministic":
        flags.append({
            "code": "DETERMINISTIC_FIT",
            "severity": "warning",
            "message": "Only 2 probes — Ka/Kb from deterministic solve, not regression",
        })

    # Probe count — only flag if truly insufficient (0-1)
    n = gutmann_result.get("n_probes", 0)
    if n < 2:
        flags.append({
            "code": "FEW_PROBES",
            "severity": "warning",
            "message": f"Only {n} probes — insufficient for Gutmann analysis",
        })

    # Negative ΔG_sp (probe sits below alkane line — physically suspect)
    if not delta_g_sp_df.empty:
        neg = delta_g_sp_df[delta_g_sp_df["delta_g_sp_kJmol"] < 0]
        for _, row in neg.iterrows():
            flags.append({
                "code": "NEGATIVE_DG_SP",
                "severity": "warning",
                "message": (
                    f"{row['probe']}: ΔG_sp = {row['delta_g_sp_kJmol']:.2f} kJ/mol < 0 "
                    f"— probe below alkane line"
                ),
            })

    # Outlier probes from leave-one-out
    if influence:
        for inf in influence:
            if inf["is_outlier"]:
                flags.append({
                    "code": "PROBE_OUTLIER",
                    "severity": "info",
                    "message": (
                        f"{inf['probe']} is an influential outlier "
                        f"(ΔR² = {inf['delta_r_squared']:+.3f} when removed)"
                    ),
                })

    # Ka or Kb negative (physically unexpected for these surfaces)
    ka = gutmann_result.get("Ka", np.nan)
    kb = gutmann_result.get("Kb", np.nan)
    if not np.isnan(ka) and ka < 0:
        flags.append({
            "code": "NEGATIVE_KA",
            "severity": "warning",
            "message": f"Ka = {ka:.4f} < 0 — negative acidity is physically suspect",
        })
    if not np.isnan(kb) and kb < 0:
        flags.append({
            "code": "NEGATIVE_KB",
            "severity": "warning",
            "message": f"Kb = {kb:.4f} < 0 — negative basicity is physically suspect",
        })

    # Van Oss checks
    if van_oss_result is not None:
        gs_minus = van_oss_result.get("gamma_s_minus_mJm2", np.nan)
        gs_plus = van_oss_result.get("gamma_s_plus_mJm2", np.nan)
        if not np.isnan(gs_minus) and gs_minus > 100:
            flags.append({
                "code": "VAN_OSS_GAMMA_HIGH",
                "severity": "warning",
                "message": (
                    f"γ_S⁻ = {gs_minus:.1f} mJ/m² exceeds physical bounds — "
                    f"known limitation of van Oss scale applied to IGC"
                ),
            })
        if not np.isnan(gs_plus) and gs_plus > 100:
            flags.append({
                "code": "VAN_OSS_GAMMA_HIGH",
                "severity": "warning",
                "message": (
                    f"γ_S⁺ = {gs_plus:.1f} mJ/m² exceeds physical bounds — "
                    f"known limitation of van Oss scale applied to IGC"
                ),
            })
    else:
        # Check if van Oss is missing due to incomplete probes
        if not delta_g_sp_df.empty:
            probes = set(delta_g_sp_df["probe"].str.lower())
            has_acid = "dichloromethane" in probes
            has_base = "ethyl acetate" in probes
            if not has_acid or not has_base:
                missing = []
                if not has_acid:
                    missing.append("DCM (monopolar acid)")
                if not has_base:
                    missing.append("ethyl acetate (monopolar base)")
                flags.append({
                    "code": "VAN_OSS_INCOMPLETE",
                    "severity": "info",
                    "message": f"Van Oss not computed — missing {', '.join(missing)}",
                })

    n_crit = sum(1 for f in flags if f["severity"] == "critical")
    n_warn = sum(1 for f in flags if f["severity"] == "warning")
    if n_crit > 0:
        summary = f"FAIL — {n_crit} critical, {n_warn} warnings"
    elif n_warn > 0:
        summary = f"REVIEW — {n_warn} warnings"
    else:
        summary = "PASS"

    return {"flags": flags, "summary": summary}


# ---------------------------------------------------------------------------
# Polar probe V_N consistency check
# ---------------------------------------------------------------------------

def check_polar_vn_consistency(
    vn_results: pd.DataFrame,
    vn_column: str = "VN",
    coverage_column: str = "coverage",
    solvent_column: str = "solvent_name",
    max_ratio: float = 10.0,
) -> list[dict]:
    """Check polar probes for non-monotonic V_N variation across coverages.

    Flags probes where max(V_N)/min(V_N) exceeds ``max_ratio`` and the
    V_N trend is not monotonically decreasing with coverage (the physically
    expected behavior for Type I/II isotherms).

    Parameters
    ----------
    vn_results : pd.DataFrame
        V_N data across all coverages.
    max_ratio : float
        Flag if max/min V_N ratio exceeds this (default 10×).

    Returns
    -------
    list of dict
        Per-probe check results: ``probe``, ``vn_min``, ``vn_max``,
        ``vn_ratio``, ``is_monotonic``, ``flagged``, ``message``.
    """
    alkane_names = {"HEXANE", "HEPTANE", "OCTANE", "NONANE", "DECANE"}
    results = []

    for solvent, grp in vn_results.groupby(solvent_column):
        if solvent.upper() in alkane_names:
            continue
        probe = _resolve_probe(solvent)
        if probe is None:
            continue

        vn_vals = grp[[coverage_column, vn_column]].dropna(subset=[vn_column])
        vn_vals = vn_vals[vn_vals[vn_column] > 0].sort_values(coverage_column)

        if len(vn_vals) < 2:
            continue

        vn_min = vn_vals[vn_column].min()
        vn_max = vn_vals[vn_column].max()
        ratio = vn_max / vn_min if vn_min > 0 else np.inf

        # Check monotonicity (V_N should decrease with increasing coverage)
        diffs = np.diff(vn_vals[vn_column].values)
        is_monotonic = bool(np.all(diffs <= 0) or np.all(diffs >= 0))

        flagged = ratio > max_ratio and not is_monotonic

        entry = {
            "probe": probe["name"],
            "solvent_name": solvent,
            "vn_min": vn_min,
            "vn_max": vn_max,
            "vn_ratio": ratio,
            "is_monotonic": is_monotonic,
            "flagged": flagged,
        }

        if flagged:
            note = ""
            if probe["name"] == "acetonitrile":
                note = " (acetonitrile is known to give erratic V_N)"
            entry["message"] = (
                f"{probe['name']}: V_N varies {ratio:.1f}× across coverages "
                f"with non-monotonic trend{note}"
            )
        else:
            entry["message"] = ""

        results.append(entry)

    return results


def check_dg_sp_variability(
    delta_g_sp_all: pd.DataFrame,
    coverage_column: str = "coverage",
    threshold_multiplier: float = 3.0,
) -> list[dict]:
    """Flag probes with unusually high ΔG_sp variability across coverages.

    Computes the standard deviation of each probe's ΔG_sp across coverages.
    Probes whose std dev exceeds ``threshold_multiplier`` × the median std
    across all probes are flagged as unreliable.

    Parameters
    ----------
    delta_g_sp_all : pd.DataFrame
        ΔG_sp for all probes across all coverages (from run_acid_base_analysis).
    threshold_multiplier : float
        Flag if std > multiplier × median_std (default 3×).

    Returns
    -------
    list of dict
        Per-probe variability: ``probe``, ``std_kJmol``, ``n_coverages``,
        ``median_std``, ``flagged``, ``message``.
    """
    if delta_g_sp_all.empty:
        return []

    # Compute per-probe std across coverages (only positive ΔG_sp)
    positive = delta_g_sp_all[delta_g_sp_all["delta_g_sp_kJmol"] > 0]
    if positive.empty:
        return []

    probe_stats = []
    for probe_name, grp in positive.groupby("probe"):
        if len(grp) < 2:
            continue
        probe_stats.append({
            "probe": probe_name,
            "std_kJmol": grp["delta_g_sp_kJmol"].std(),
            "n_coverages": len(grp),
        })

    if len(probe_stats) < 2:
        return probe_stats

    stds = [p["std_kJmol"] for p in probe_stats]
    median_std = float(np.median(stds))

    results = []
    for ps in probe_stats:
        ps["median_std"] = median_std
        flagged = ps["std_kJmol"] > threshold_multiplier * median_std and median_std > 0
        ps["flagged"] = flagged
        if flagged:
            ps["message"] = (
                f"{ps['probe']}: ΔG_sp std = {ps['std_kJmol']:.2f} kJ/mol "
                f"({ps['std_kJmol']/median_std:.1f}× median) — unreliable"
            )
        else:
            ps["message"] = ""
        results.append(ps)

    return results


# ---------------------------------------------------------------------------
# Van Oss acid-base components
# ---------------------------------------------------------------------------

# Van Oss probe parameters: acid (γ_L⁺) and base (γ_L⁻) components
# of liquid surface energy in mJ/m².
# Source: van Oss, Good & Chaudhury (1988) Langmuir 4, 884-891
#
# KNOWN LIMITATION: The van Oss scale is referenced to water (γ⁺ = γ⁻ =
# 25.5 mJ/m²) and was developed for contact angle measurements.  When
# applied to IGC ΔG_sp data, the small γ_L⁺ for DCM (5.2 mJ/m²) acts as
# a divisor that amplifies γ_S⁻ to physically unreasonable values (>100
# mJ/m²) for strongly basic surfaces.  This is a known limitation of the
# van Oss framework applied to IGC, not a calculation error.  Relative
# trends between formulations are more reliable than absolute values.
# Alternative γ_L⁺/γ_L⁻ scales (e.g. Della Volpe & Siboni 2004) may
# give more physically reasonable magnitudes.
VAN_OSS_PROBE_PARAMS = {
    "dichloromethane": {"gamma_l_plus": 5.2, "gamma_l_minus": 0.0},
    "ethyl acetate": {"gamma_l_plus": 0.0, "gamma_l_minus": 19.2},
}


def van_oss_gamma_sp(
    delta_g_sp_acid: float,
    delta_g_sp_base: float,
    a_acid: float,
    a_base: float,
    gamma_l_plus_acid: float,
    gamma_l_minus_base: float,
) -> dict:
    """Compute γ_S⁺ and γ_S⁻ from monopolar acid and base probe ΔG_sp.

    For a monopolar acid probe (e.g. DCM, γ_L⁻ ≈ 0):
        ΔG_sp = 2 · a · N_A · √(γ_L⁺ · γ_S⁻)
        γ_S⁻ = (ΔG_sp / (2 · a · N_A))² / γ_L⁺

    For a monopolar base probe (e.g. ethyl acetate, γ_L⁺ ≈ 0):
        ΔG_sp = 2 · a · N_A · √(γ_L⁻ · γ_S⁺)
        γ_S⁺ = (ΔG_sp / (2 · a · N_A))² / γ_L⁻

    Parameters
    ----------
    delta_g_sp_acid : float
        ΔG_sp for the acid probe (J/mol). Positive = surface is basic.
    delta_g_sp_base : float
        ΔG_sp for the base probe (J/mol). Positive = surface is acidic.
    a_acid, a_base : float
        Cross-sectional areas (m²).
    gamma_l_plus_acid : float
        Acid component of liquid surface energy for acid probe (mJ/m²).
    gamma_l_minus_base : float
        Base component of liquid surface energy for base probe (mJ/m²).

    Returns
    -------
    dict
        ``gamma_s_minus`` (mJ/m²) — surface base component (from acid probe),
        ``gamma_s_plus`` (mJ/m²) — surface acid component (from base probe),
        ``gamma_s_ab`` (mJ/m²) — acid-base component = 2·√(γ_S⁺·γ_S⁻).
    """
    # γ_S⁻ from acid probe (DCM probes surface basicity)
    if delta_g_sp_acid > 0 and gamma_l_plus_acid > 0:
        term = delta_g_sp_acid / (2.0 * a_acid * N_AVOGADRO)
        gamma_s_minus = (term ** 2) / (gamma_l_plus_acid * 1e-3)  # → J/m²
        gamma_s_minus_mJm2 = gamma_s_minus * 1e3
    else:
        gamma_s_minus_mJm2 = 0.0 if delta_g_sp_acid <= 0 else np.nan

    # γ_S⁺ from base probe (EtAc probes surface acidity)
    if delta_g_sp_base > 0 and gamma_l_minus_base > 0:
        term = delta_g_sp_base / (2.0 * a_base * N_AVOGADRO)
        gamma_s_plus = (term ** 2) / (gamma_l_minus_base * 1e-3)  # → J/m²
        gamma_s_plus_mJm2 = gamma_s_plus * 1e3
    else:
        gamma_s_plus_mJm2 = 0.0 if delta_g_sp_base <= 0 else np.nan

    # γ_S^AB = 2·√(γ_S⁺ · γ_S⁻)
    if (not np.isnan(gamma_s_plus_mJm2) and not np.isnan(gamma_s_minus_mJm2)
            and gamma_s_plus_mJm2 >= 0 and gamma_s_minus_mJm2 >= 0):
        gamma_s_ab = 2.0 * np.sqrt(gamma_s_plus_mJm2 * gamma_s_minus_mJm2)
    else:
        gamma_s_ab = np.nan

    # Work of adhesion per unit area (more physically interpretable)
    W_a_acid = delta_g_sp_acid / (N_AVOGADRO * a_acid) * 1e3 if a_acid > 0 else np.nan
    W_a_base = delta_g_sp_base / (N_AVOGADRO * a_base) * 1e3 if a_base > 0 else np.nan

    return {
        "gamma_s_plus_mJm2": gamma_s_plus_mJm2,
        "gamma_s_minus_mJm2": gamma_s_minus_mJm2,
        "gamma_s_ab_mJm2": gamma_s_ab,
        "W_a_acid_mJm2": W_a_acid,   # specific work of adhesion per area
        "W_a_base_mJm2": W_a_base,
    }


def compute_van_oss(delta_g_sp_df: pd.DataFrame) -> dict | None:
    """Compute van Oss γ_S⁺/γ_S⁻ if DCM and ethyl acetate are available.

    Parameters
    ----------
    delta_g_sp_df : pd.DataFrame
        Output of :func:`calculate_delta_g_sp`.

    Returns
    -------
    dict or None
        Van Oss results, or None if the required probes are missing.
    """
    if delta_g_sp_df.empty:
        return None

    probes_available = set(delta_g_sp_df["probe"].str.lower())
    if not {"dichloromethane", "ethyl acetate"}.issubset(probes_available):
        return None

    dcm_row = delta_g_sp_df[delta_g_sp_df["probe"] == "dichloromethane"].iloc[0]
    etac_row = delta_g_sp_df[delta_g_sp_df["probe"] == "ethyl acetate"].iloc[0]

    dcm_params = VAN_OSS_PROBE_PARAMS["dichloromethane"]
    etac_params = VAN_OSS_PROBE_PARAMS["ethyl acetate"]

    return van_oss_gamma_sp(
        delta_g_sp_acid=dcm_row["delta_g_sp_Jmol"],
        delta_g_sp_base=etac_row["delta_g_sp_Jmol"],
        a_acid=dcm_row.get("a_cross", get_probe("dichloromethane")["a_cross"]),
        a_base=etac_row.get("a_cross", get_probe("ethyl acetate")["a_cross"]),
        gamma_l_plus_acid=dcm_params["gamma_l_plus"],
        gamma_l_minus_base=etac_params["gamma_l_minus"],
    )


# ---------------------------------------------------------------------------
# Top-level orchestrator
# ---------------------------------------------------------------------------

def run_acid_base_analysis(
    vn_results: pd.DataFrame,
    temperature: float,
    coverages: list[float] | None = None,
    vn_column: str = "VN",
    coverage_column: str = "coverage",
    solvent_column: str = "solvent_name",
    exclude_probes: list[str] | None = None,
    include_van_oss: bool = False,
) -> dict:
    """Run complete acid-base analysis across all coverages.

    Parameters
    ----------
    vn_results : pd.DataFrame
        V_N data for all probes (alkane + polar) at target coverages.
        Must contain ``solvent_column``, ``vn_column``, ``coverage_column``.
    temperature : float
        Column temperature (K).
    coverages : list of float, optional
        Coverages to analyze.  Defaults to all unique coverages in data.
    vn_column, coverage_column, solvent_column : str
        Column names.
    exclude_probes : list of str, optional
        Probes to exclude from Gutmann regression.
    include_van_oss : bool
        If True, compute van Oss γ_S⁺/γ_S⁻ from monopolar probes and
        include in profile and QC.  Default False — van Oss values can
        be physically unreasonable for IGC data (see module docstring).

    Returns
    -------
    dict
        Keys:
        - ``profile`` : pd.DataFrame — Ka, Kb, R² at each coverage
        - ``delta_g_sp`` : pd.DataFrame — per-probe ΔG_sp at each coverage
        - ``schultz_lines`` : dict — Schultz line params per coverage
        - ``gutmann_results`` : dict — full Gutmann result per coverage
        - ``van_oss_results`` : dict — van Oss γ_S⁺/γ_S⁻ per coverage
          (only populated when ``include_van_oss=True``)
        - ``polar_vn_checks`` : list — polar probe V_N consistency flags
        - ``qc`` : dict — quality checks (flags, summary)
    """
    if coverages is None:
        coverages = sorted(vn_results[coverage_column].unique())

    # Polar probe V_N consistency check (across all coverages)
    polar_vn_checks = check_polar_vn_consistency(
        vn_results, vn_column=vn_column,
        coverage_column=coverage_column, solvent_column=solvent_column,
    )

    profile_rows = []
    all_delta_g_sp = []
    schultz_lines = {}
    gutmann_results = {}
    van_oss_results = {}
    all_flags = []

    # Add V_N consistency flags
    for check in polar_vn_checks:
        if check["flagged"]:
            all_flags.append({
                "code": "POLAR_VN_INCONSISTENT",
                "severity": "warning",
                "message": check["message"],
                "coverage": "all",
            })

    for cov in coverages:
        cov_data = vn_results[
            abs(vn_results[coverage_column] - cov) < 1e-6
        ].copy()

        if cov_data.empty:
            continue

        # Schultz reference line from alkanes
        sline = schultz_reference_line(
            cov_data, temperature,
            vn_column=vn_column, solvent_column=solvent_column,
        )
        schultz_lines[cov] = sline

        # ΔG_sp for polar probes
        dg_sp = calculate_delta_g_sp(
            cov_data, sline, temperature,
            vn_column=vn_column, solvent_column=solvent_column,
        )

        if not dg_sp.empty:
            dg_sp[coverage_column] = cov
            all_delta_g_sp.append(dg_sp)

        # Gutmann regression
        gutmann = gutmann_ka_kb(dg_sp, exclude_probes=exclude_probes)
        gutmann_results[cov] = gutmann

        # Van Oss (opt-in: only when explicitly requested)
        vo = None
        if include_van_oss:
            vo = compute_van_oss(dg_sp)
            if vo is not None:
                van_oss_results[cov] = vo

        # Leave-one-out (only meaningful for 3+ probes)
        influence = leave_one_out_influence(dg_sp)

        # QC (van Oss checks only when opted in)
        qc_cov = acid_base_quality_checks(
            gutmann, dg_sp, influence,
            van_oss_result=vo if include_van_oss else None,
        )
        for flag in qc_cov["flags"]:
            flag["coverage"] = cov
            all_flags.append(flag)

        # Profile row
        row = {
            coverage_column: cov,
            "Ka": gutmann["Ka"],
            "Kb": gutmann["Kb"],
            "Kb_Ka_ratio": gutmann["Kb_Ka_ratio"],
            "r_squared": gutmann["r_squared"],
            "fit_method": gutmann["fit_method"],
            "n_probes": gutmann["n_probes"],
            "probes_used": ", ".join(gutmann["probes_used"]),
            "gamma_d_schultz_mJm2": sline["gamma_d_schultz_mJm2"],
            "schultz_r_squared": sline["r_squared"],
        }

        # Add van Oss columns only when opted in
        if include_van_oss:
            if vo is not None:
                row["gamma_s_plus_mJm2"] = vo["gamma_s_plus_mJm2"]
                row["gamma_s_minus_mJm2"] = vo["gamma_s_minus_mJm2"]
                row["gamma_s_ab_mJm2"] = vo["gamma_s_ab_mJm2"]
            else:
                row["gamma_s_plus_mJm2"] = np.nan
                row["gamma_s_minus_mJm2"] = np.nan
                row["gamma_s_ab_mJm2"] = np.nan

        profile_rows.append(row)

    profile = pd.DataFrame(profile_rows)
    delta_g_sp_all = pd.concat(all_delta_g_sp, ignore_index=True) if all_delta_g_sp else pd.DataFrame()

    # ΔG_sp variability check across coverages
    variability_checks = check_dg_sp_variability(
        delta_g_sp_all, coverage_column=coverage_column,
    )
    for vc in variability_checks:
        if vc["flagged"]:
            all_flags.append({
                "code": "PROBE_DG_SP_VARIABLE",
                "severity": "warning",
                "message": vc["message"],
                "coverage": "all",
            })

    # Overall QC summary
    n_crit = sum(1 for f in all_flags if f["severity"] == "critical")
    n_warn = sum(1 for f in all_flags if f["severity"] == "warning")
    if n_crit > 0:
        overall = f"FAIL — {n_crit} critical, {n_warn} warnings across {len(coverages)} coverages"
    elif n_warn > 0:
        overall = f"REVIEW — {n_warn} warnings across {len(coverages)} coverages"
    else:
        overall = f"PASS — {len(coverages)} coverages analyzed"

    return {
        "profile": profile,
        "delta_g_sp": delta_g_sp_all,
        "schultz_lines": schultz_lines,
        "gutmann_results": gutmann_results,
        "van_oss_results": van_oss_results,
        "polar_vn_checks": polar_vn_checks,
        "variability_checks": variability_checks,
        "qc": {"flags": all_flags, "summary": overall},
    }
