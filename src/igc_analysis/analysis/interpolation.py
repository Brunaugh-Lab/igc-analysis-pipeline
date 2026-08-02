"""Coverage interpolation for inverse chromatography retention-volume data.

Raw IGC injections target specific surface coverages (e.g. 0.15, 0.12,
0.10 n/nm), but the *actual* coverage achieved differs slightly from the
target.  This module interpolates the measured retention volumes to exact
target coverage values using piecewise linear interpolation
software.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d


def interpolate_to_coverage(
    injections: pd.DataFrame,
    target_coverages: list[float] | np.ndarray,
    vn_column: str = "net_retention_volume",
    coverage_column: str = "actual_coverage",
    solvent_column: str = "solvent_name",
    method: str = "linear",
    extrapolate: bool = False,
) -> pd.DataFrame:
    """Interpolate retention volumes to exact target coverage values.

    For each solvent, fits an interpolation through V_N vs actual_coverage
    and evaluates at the specified target coverages.

    Parameters
    ----------
    injections : pd.DataFrame
        Must contain columns for solvent name, actual coverage, and V_N.
    target_coverages : array-like
        Desired coverage values to interpolate to.
    vn_column : str
        Column name for retention volume values.
    coverage_column : str
        Column name for actual (measured) coverage values.
    solvent_column : str
        Column name for solvent/probe identity.
    method : str
        Interpolation method: ``"linear"`` (default) or ``"cubic"``.

    Returns
    -------
    pd.DataFrame
        Columns: ``solvent_name``, ``target_coverage``,
        ``interpolated_VN``.  One row per (solvent × target_coverage)
        combination.  Rows where the target coverage is outside the
        measured range are marked NaN.
    """
    target_coverages = np.asarray(target_coverages)
    rows: list[dict] = []

    for solvent, group in injections.groupby(solvent_column):
        # Sort by coverage (ascending)
        g = group.sort_values(coverage_column)
        x = g[coverage_column].values
        y = g[vn_column].values

        # Remove NaN/inf
        mask = np.isfinite(x) & np.isfinite(y)
        x, y = x[mask], y[mask]

        if len(x) < 2:
            # Can't interpolate with fewer than 2 points
            for tc in target_coverages:
                rows.append({
                    "solvent_name": solvent,
                    "target_coverage": tc,
                    "interpolated_VN": np.nan,
                })
            continue

        # Build interpolator
        kind = method if len(x) >= 4 or method == "linear" else "linear"
        if extrapolate:
            f = interp1d(x, y, kind=kind, fill_value="extrapolate", bounds_error=False)
        else:
            f = interp1d(x, y, kind=kind, fill_value=np.nan, bounds_error=False)

        for tc in target_coverages:
            rows.append({
                "solvent_name": solvent,
                "target_coverage": tc,
                "interpolated_VN": float(f(tc)),
            })

    return pd.DataFrame(rows)
