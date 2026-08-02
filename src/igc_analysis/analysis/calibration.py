"""Calibration calculations from declared neutral inputs.

The legacy helper retained here evaluates a two-parameter power-law model:

    moles_injected = C1 × area^C2           (SI, mol)

Actual fractional surface coverage is then:

    θ = moles_injected / n_monolayer

where:

    n_monolayer = BET_SSA × sample_mass / (N_A × a_cross)

- ``BET_SSA`` is the BET specific surface area (m²/g)
- ``sample_mass`` in grams
- ``N_A`` is Avogadro's number
- ``a_cross`` is the solvent molecular cross-sectional area (m²)
"""

from __future__ import annotations

from igc_analysis.constants import N_AVOGADRO


def moles_from_area(peak_area: float, C1: float, C2: float) -> float:
    """Compute moles from integrated area using a declared power law.

    Parameters
    ----------
    peak_area : float
        Integrated detector area in the calibration's declared area unit.
    C1, C2 : float
        Declared scale and exponent.

    Returns
    -------
    float
        Moles injected (SI, mol).
    """
    if peak_area <= 0:
        return 0.0
    return C1 * (peak_area ** C2)


def monolayer_capacity(bet_ssa_m2g: float, mass_g: float,
                       cross_section_m2: float) -> float:
    """Compute monolayer capacity in moles.

    Parameters
    ----------
    bet_ssa_m2g : float
        BET specific surface area (m²/g).
    mass_g : float
        Sample mass (g).
    cross_section_m2 : float
        Molecular cross-sectional area of the probe solvent (m²).

    Returns
    -------
    float
        Monolayer capacity (mol).
    """
    return bet_ssa_m2g * mass_g / (N_AVOGADRO * cross_section_m2)


def actual_coverage(peak_area: float, C1: float, C2: float,
                    bet_ssa_m2g: float, mass_g: float,
                    cross_section_m2: float) -> float:
    """Compute actual fractional surface coverage from peak area.

    Parameters
    ----------
    peak_area : float
        Integrated FID peak area.
    C1, C2 : float
        C-polynomial coefficients.
    bet_ssa_m2g : float
        BET specific surface area (m²/g).
    mass_g : float
        Sample mass (g).
    cross_section_m2 : float
        Solvent molecular cross-sectional area (m²).

    Returns
    -------
    float
        Fractional surface coverage θ (dimensionless).
    """
    moles = moles_from_area(peak_area, C1, C2)
    n_mono = monolayer_capacity(bet_ssa_m2g, mass_g, cross_section_m2)
    if n_mono <= 0:
        return 0.0
    return moles / n_mono
