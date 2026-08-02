"""Net retention-volume calculations for inverse gas chromatography.

Converts raw retention times into thermodynamically meaningful retention
volumes, accounting for dead time, gas compressibility, temperature, and
sample mass.
"""

from __future__ import annotations

import numpy as np


def james_martin_correction(p_inlet: float, p_outlet: float) -> float:
    """James-Martin compressibility correction factor *j*.

    Accounts for the pressure drop across the packed column, which causes
    gas velocity (and thus probe residence time) to vary along the column
    length.

    .. math::

        j = \\frac{3}{2} \\cdot
            \\frac{(P_i/P_o)^2 - 1}{(P_i/P_o)^3 - 1}

    Parameters
    ----------
    p_inlet : float
        Inlet pressure (any consistent unit).
    p_outlet : float
        Outlet pressure (same unit as *p_inlet*).

    Returns
    -------
    float
        Dimensionless correction factor (0 < j ≤ 1).

    Notes
    -----
    When *p_inlet* ≈ *p_outlet*, j → 1.0.
    """
    if p_outlet <= 0:
        raise ValueError(f"Outlet pressure must be positive, got {p_outlet}")
    if p_inlet < p_outlet:
        raise ValueError(f"Inlet pressure ({p_inlet}) must be ≥ outlet ({p_outlet})")

    r = p_inlet / p_outlet
    if abs(r - 1.0) < 1e-10:
        return 1.0
    return 1.5 * (r**2 - 1) / (r**3 - 1)


def net_retention_time(
    t_R: float,
    t_0: float,
) -> float:
    """Compute net retention time.

    Parameters
    ----------
    t_R : float
        Total retention time of the probe (minutes).
    t_0 : float
        Dead time / holdup time from methane marker (minutes).

    Returns
    -------
    float
        Net retention time *t_R − t_0* (minutes).
    """
    return t_R - t_0


def net_retention_volume(
    t_R: float,
    t_0: float,
    flow_rate: float,
    j_correction: float = 1.0,
    temperature_col_K: float | None = None,
    temperature_flow_K: float | None = None,
) -> float:
    """Compute net retention volume V_N.

    .. math::

        V_N = (t_R - t_0) \\cdot F \\cdot j \\cdot \\frac{T_{col}}{T_{flow}}

    Parameters
    ----------
    t_R : float
        Total retention time (minutes).
    t_0 : float
        Dead time (minutes).
    flow_rate : float
        Carrier gas exit flow rate (mL/min).
    j_correction : float
        James-Martin compressibility factor (default 1.0).
    temperature_col_K : float, optional
        Column temperature (K). If both temperatures given, applies
        temperature correction.
    temperature_flow_K : float, optional
        Flowmeter temperature (K).

    Returns
    -------
    float
        Net retention volume V_N (mL).
    """
    t_net = t_R - t_0
    v_n = t_net * flow_rate * j_correction

    # Temperature correction: correct flow to column temperature
    if temperature_col_K is not None and temperature_flow_K is not None:
        v_n *= temperature_col_K / temperature_flow_K

    return v_n


def specific_retention_volume(
    v_n: float,
    sample_mass_g: float,
) -> float:
    """Normalize retention volume by sample mass.

    Parameters
    ----------
    v_n : float
        Net retention volume (mL).
    sample_mass_g : float
        Sample mass (g).

    Returns
    -------
    float
        Specific net retention volume V_N (mL/g).
    """
    if sample_mass_g <= 0:
        raise ValueError(f"Sample mass must be positive, got {sample_mass_g}")
    return v_n / sample_mass_g
