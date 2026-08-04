"""Physical constants and reference values for inverse chromatography.

Sources:
    - Dorris & Gray (1980) J. Colloid Interface Sci. 77(2), 353-362
    - Schultz, Lavielle & Martin (1987) J. Adhesion 23, 45-60
    - Fowkes (1964) Ind. Eng. Chem. 56(12), 40-52
"""

# ---------------------------------------------------------------------------
# Fundamental constants
# ---------------------------------------------------------------------------

R_GAS: float = 8.314462618  # J/(mol·K) — NIST 2018 CODATA
"""Universal gas constant."""

N_AVOGADRO: float = 6.02214076e23  # mol⁻¹ — exact (2019 SI)
"""Avogadro's number."""

T_STANDARD_K: float = 273.15
"""Reference temperature for neutral ``flow_standard`` values (K)."""

# ---------------------------------------------------------------------------
# Dorris-Gray reference surface (polyethylene / CH₂)
# ---------------------------------------------------------------------------

A_CH2: float = 6.0e-20  # m²
"""Cross-sectional area of a –CH₂– group on a close-packed polyethylene
surface. Value from Dorris & Gray (1980)."""

GAMMA_CH2_REF: float = 35.6e-3  # J/m² at 293.15 K (20 °C)
"""Dispersive surface energy of the polyethylene reference surface at 20 °C.
Value from Dorris & Gray (1980)."""

GAMMA_CH2_TEMP_COEFF: float = -0.058e-3  # J/(m²·K)
"""Temperature coefficient for γ_CH₂.  dγ/dT ≈ −0.058 mJ/(m²·K)."""


def gamma_ch2(temperature_K: float) -> float:
    """Dispersive surface energy of the CH₂ reference surface at *temperature_K*.

    Uses a linear temperature correction:
        γ_CH₂(T) = 35.6 − 0.058·(T − 293.15) mJ/m²

    Parameters
    ----------
    temperature_K : float
        Column temperature in Kelvin.

    Returns
    -------
    float
        γ_CH₂ in J/m² (SI units — divide by 1e-3 for mJ/m²).
    """
    return GAMMA_CH2_REF + GAMMA_CH2_TEMP_COEFF * (temperature_K - 293.15)
