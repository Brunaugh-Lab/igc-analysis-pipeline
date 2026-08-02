"""Dispersive surface energy (γ_d) analysis using the Dorris-Gray method.

Method:
  For a homologous alkane series (e.g. C8–C10), the net retention volume
  V_N is measured at one or more surface coverages. Plotting
  RT·ln(V_N) vs carbon number yields a linear relationship whose slope
  is ΔG_CH₂ — the free energy of adsorption per methylene group.

  The dispersive surface energy is then:

      γ_d = (ΔG_CH₂)² / (4 · N_A² · a_CH₂² · γ_CH₂)

Input:
  Retention volumes for alkane probes at specified surface coverages.

Output:
  γ_d (mJ/m²) at each surface coverage; regression diagnostics.

References:
  - Dorris & Gray (1980) J. Colloid Interface Sci. 77(2), 353-362
  - Schultz, Lavielle & Martin (1987) J. Adhesion 23, 45-60
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

from igc_sea.constants import R_GAS, N_AVOGADRO, A_CH2, gamma_ch2
from igc_sea.analysis.retention import (
    net_retention_volume,
    specific_retention_volume,
)


# ---------------------------------------------------------------------------
# Core Dorris-Gray calculation
# ---------------------------------------------------------------------------

def dorris_gray_gamma_d(
    df: pd.DataFrame,
    temperature: float,
) -> pd.DataFrame:
    """Calculate dispersive surface energy γ_d via the Dorris-Gray method.

    Performs linear regression of RT·ln(V_N) vs carbon number for the
    alkane series at each coverage.  The slope yields γ_d through:

        γ_d = (slope)² / (4 · N_A² · a_CH₂² · γ_CH₂)

    Parameters
    ----------
    df : pd.DataFrame
        Required columns:

        - ``coverage`` : float — surface coverage n/nm
        - ``carbon_number`` : int — alkane carbon count (8, 9, 10)
        - ``VN`` : float — specific net retention volume (mL/g)

        May also contain ``sample_name`` for labeling.
    temperature : float
        Column temperature (K).

    Returns
    -------
    pd.DataFrame
        Columns: ``coverage``, ``gamma_d_mJm2``, ``r_squared``,
        ``slope_Jmol``, ``intercept``, ``n_alkanes``.
    """
    # Filter to alkanes only
    alkanes = df[df["carbon_number"].notna()].copy()

    if alkanes.empty:
        raise ValueError("No alkane data found (carbon_number column is all NaN)")

    vn = pd.to_numeric(alkanes["VN"], errors="coerce")
    invalid_vn = ~np.isfinite(vn) | (vn <= 0)
    if invalid_vn.any():
        raise ValueError(
            "VN must be finite and strictly positive for every alkane row; "
            f"found {int(invalid_vn.sum())} invalid value(s)"
        )
    alkanes["VN"] = vn

    # RT·ln(V_N) in J/mol
    # V_N is in mL/g — ln(V_N) is dimensionless (the reference state cancels)
    alkanes = alkanes.copy()
    alkanes["RT_ln_VN"] = R_GAS * temperature * np.log(alkanes["VN"])

    results: list[dict] = []

    for coverage, group in alkanes.groupby("coverage"):
        g = group.dropna(subset=["VN", "RT_ln_VN"]).sort_values("carbon_number")

        if len(g) < 2:
            results.append({
                "coverage": coverage,
                "gamma_d_mJm2": np.nan,
                "r_squared": np.nan,
                "slope_Jmol": np.nan,
                "intercept": np.nan,
                "n_alkanes": len(g),
            })
            continue

        # Linear regression: RT·ln(V_N) = slope·n + intercept
        slope, intercept, r_value, p_value, std_err = stats.linregress(
            g["carbon_number"].values, g["RT_ln_VN"].values
        )

        # Dorris-Gray equation
        # ΔG_CH₂ = slope (J/mol per CH₂ group)
        delta_g_ch2 = slope  # J/mol

        # γ_d = (ΔG_CH₂)² / (4 · N_A² · a_CH₂² · γ_CH₂(T))
        gamma_ch2_T = gamma_ch2(temperature)  # J/m²
        gamma_d = delta_g_ch2**2 / (4.0 * N_AVOGADRO**2 * A_CH2**2 * gamma_ch2_T)

        # Convert to mJ/m²
        gamma_d_mJm2 = gamma_d * 1e3

        results.append({
            "coverage": coverage,
            "gamma_d_mJm2": gamma_d_mJm2,
            "r_squared": r_value**2,
            "slope_Jmol": slope,
            "intercept": intercept,
            "n_alkanes": len(g),
        })

    return pd.DataFrame(results).sort_values("coverage").reset_index(drop=True)


# ---------------------------------------------------------------------------
# Convenience wrappers
# ---------------------------------------------------------------------------

def calculate_net_retention_volume(
    retention_time: float,
    dead_time: float,
    flow_rate: float,
    temperature: float,
    james_martin_correction: float = 1.0,
) -> float:
    """Calculate net retention volume V_N from retention time.

    Convenience wrapper around :func:`igc_sea.analysis.retention.net_retention_volume`.

    Parameters
    ----------
    retention_time : float
        Total retention time (min).
    dead_time : float
        Dead time / holdup time (min).
    flow_rate : float
        Carrier gas flow rate (mL/min), corrected to column temperature.
    temperature : float
        Column temperature (K).
    james_martin_correction : float
        James-Martin compressibility correction factor j (default 1.0).

    Returns
    -------
    float
        Net retention volume V_N (mL).
    """
    return net_retention_volume(
        t_R=retention_time,
        t_0=dead_time,
        flow_rate=flow_rate,
        j_correction=james_martin_correction,
    )


def gamma_d_profile(
    df: pd.DataFrame,
    temperature: float,
) -> pd.DataFrame:
    """Calculate γ_d as a function of surface coverage.

    Thin wrapper around :func:`dorris_gray_gamma_d` that returns only
    the coverage–γ_d pairs.

    Parameters
    ----------
    df : pd.DataFrame
        Alkane retention data with columns ``coverage``, ``carbon_number``,
        ``VN``.
    temperature : float
        Column temperature (K).

    Returns
    -------
    pd.DataFrame
        Columns: ``coverage``, ``gamma_d_mJm2``.
    """
    result = dorris_gray_gamma_d(df, temperature)
    return result[["coverage", "gamma_d_mJm2"]].copy()
