"""Shared utilities: probe molecule database, unit conversions, helpers."""

import pandas as pd

# =============================================================================
# Probe Molecule Reference Database
# =============================================================================
#
# Standard IGC probe molecules with physicochemical properties needed for
# surface energy calculations.
#
# Properties:
#   name            — Common name
#   formula         — Molecular formula
#   mw              — Molecular weight (g/mol)
#   category        — alkane / acid / base / amphoteric
#   carbon_number   — Number of carbons (alkanes only, else None)
#   a_cross         — Cross-sectional area (m², for BET/coverage calculations)
#   gamma_l_d       — Dispersive surface energy of liquid probe (mJ/m²)
#   dn              — Gutmann donor number (kJ/mol)
#   an_star         — Modified acceptor number AN* (kJ/mol)
#
# Sources:
#   - Schultz, Lavielle, Martin (1987) J. Adhesion
#   - Gutmann (1978) The Donor-Acceptor Approach to Molecular Interactions
#   - Della Volpe & Siboni (1997) J. Colloid Interface Sci.

PROBE_MOLECULES = pd.DataFrame([
    # Alkane series (Dorris-Gray reference line)
    {
        "name": "hexane",
        "formula": "C6H14",
        "mw": 86.18,
        "category": "alkane",
        "carbon_number": 6,
        "a_cross": 51.5e-20,
        "gamma_l_d": 18.4,
        "dn": 0.0,
        "an_star": 0.0,
    },
    {
        "name": "heptane",
        "formula": "C7H16",
        "mw": 100.21,
        "category": "alkane",
        "carbon_number": 7,
        "a_cross": 57.0e-20,
        "gamma_l_d": 20.3,
        "dn": 0.0,
        "an_star": 0.0,
    },
    {
        "name": "octane",
        "formula": "C8H18",
        "mw": 114.23,
        "category": "alkane",
        "carbon_number": 8,
        "a_cross": 63.0e-20,
        "gamma_l_d": 21.3,
        "dn": 0.0,
        "an_star": 0.0,
    },
    {
        "name": "nonane",
        "formula": "C9H20",
        "mw": 128.26,
        "category": "alkane",
        "carbon_number": 9,
        "a_cross": 69.0e-20,
        "gamma_l_d": 22.7,
        "dn": 0.0,
        "an_star": 0.0,
    },
    {
        "name": "decane",
        "formula": "C10H22",
        "mw": 142.29,
        "category": "alkane",
        "carbon_number": 10,
        "a_cross": 75.0e-20,
        "gamma_l_d": 23.8,
        "dn": 0.0,
        "an_star": 0.0,
    },
    # Polar probes (Gutmann acid-base analysis). These legacy defaults are
    # retained for calculation-level compatibility; release-facing workflows
    # should use source-attributed values from the neutral bundle.
    {
        "name": "dichloromethane",
        "formula": "CH2Cl2",
        "mw": 84.93,
        "category": "acid",
        "carbon_number": None,
        "a_cross": 24.5e-20,
        "gamma_l_d": 24.5,
        "dn": 0.0,
        "an_star": 16.3,
    },
    {
        "name": "chloroform",
        "formula": "CHCl3",
        "mw": 119.38,
        "category": "acid",
        "carbon_number": None,
        "a_cross": 44.0e-20,
        "gamma_l_d": 27.1,
        "dn": 0.0,
        "an_star": 22.6,
    },
    {
        "name": "ethyl acetate",
        "formula": "C4H8O2",
        "mw": 88.11,
        "category": "base",
        "carbon_number": None,
        "a_cross": 33.0e-20,
        "gamma_l_d": 19.6,
        "dn": 71.5,
        "an_star": 6.3,
    },
    {
        "name": "tetrahydrofuran",
        "formula": "C4H8O",
        "mw": 72.11,
        "category": "base",
        "carbon_number": None,
        "a_cross": 45.0e-20,
        "gamma_l_d": 22.5,
        "dn": 83.7,
        "an_star": 2.1,
    },
    {
        "name": "acetone",
        "formula": "C3H6O",
        "mw": 58.08,
        "category": "base",
        "carbon_number": None,
        "a_cross": 34.0e-20,
        "gamma_l_d": 16.5,
        "dn": 71.1,
        "an_star": 10.5,
    },
    {
        "name": "acetonitrile",
        "formula": "C2H3N",
        "mw": 41.05,
        "category": "amphoteric",
        "carbon_number": None,
        "a_cross": 21.4e-20,
        "gamma_l_d": 27.5,
        "dn": 59.0,
        "an_star": 19.7,
    },
    {
        "name": "ethanol",
        "formula": "C2H5OH",
        "mw": 46.07,
        "category": "amphoteric",
        "carbon_number": None,
        "a_cross": 35.3e-20,
        "gamma_l_d": 21.1,
        "dn": 80.0,
        "an_star": 37.1,
    },
    {
        "name": "diethyl ether",
        "formula": "C4H10O",
        "mw": 74.12,
        "category": "base",
        "carbon_number": None,
        "a_cross": 47.0e-20,
        "gamma_l_d": 15.0,
        "dn": 80.3,
        "an_star": 3.9,
    },
])


def get_alkanes() -> pd.DataFrame:
    """Return alkane probe molecules sorted by carbon number."""
    return (
        PROBE_MOLECULES.query("category == 'alkane'")
        .sort_values("carbon_number")
        .reset_index(drop=True)
    )


def get_polar_probes() -> pd.DataFrame:
    """Return polar (non-alkane) probe molecules."""
    return (
        PROBE_MOLECULES.query("category != 'alkane'")
        .sort_values("name")
        .reset_index(drop=True)
    )


def get_probe(name: str) -> pd.Series:
    """Look up a single probe molecule by name (case-insensitive).

    Parameters
    ----------
    name : str
        Probe molecule name (e.g. "octane", "dichloromethane").

    Returns
    -------
    pd.Series
        Row from PROBE_MOLECULES.

    Raises
    ------
    KeyError
        If probe name not found in the database.
    """
    matches = PROBE_MOLECULES[PROBE_MOLECULES["name"].str.lower() == name.lower()]
    if matches.empty:
        available = ", ".join(PROBE_MOLECULES["name"].tolist())
        raise KeyError(f"Probe '{name}' not found. Available: {available}")
    return matches.iloc[0]
