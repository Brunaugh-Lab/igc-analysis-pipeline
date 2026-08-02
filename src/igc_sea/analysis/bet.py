"""BET specific surface area calculations for alkane isotherm data.

Public workflows supply calibrated injections and source-attributed probe
properties through the neutral contract. This module owns only the scientific
transformations and reportability checks.

  1. Peak detection → peak area and retention time for each injection
  2. Declared calibration: area → moles injected
  3. Extended Antoine equation (probe-specific) → P_sat → P/P0
  4. Net retention volume V_N at each concentration, using per-injection
     measured column temperature, flow and pressure drop
  5. Adsorption isotherm by cumulative integration of V_N(c)/m, with a
     selectable zero-pressure origin treatment
  6. BET linearization in the 0.05 ≤ P/P0 ≤ 0.35 window
  7. SSA = n_monolayer × N_A × a_cross (probe-specific a_cross)

References:
    - Brunauer, Emmett & Teller (1938) JACS 60(2), 309-319
    - Thielmann (2004) J. Chromatogr. A 1037, 115-131 (IGC-BET)
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from igc_sea.constants import R_GAS, N_AVOGADRO
from igc_sea.analysis.probes import ProbeProperties, ProbeSelection


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class InjectionResult:
    """Results from processing a single isotherm injection."""

    injection_number: int
    target_coverage: float | None
    peak_max_time: float        # min
    peak_cofm_time: float       # min
    peak_area: float            # µV·min
    peak_height: float          # µV
    n_injected_mol: float       # mol (from declared calibration)
    V_N_mL: float               # net retention volume (mL)
    concentration_mol_m3: float # gas-phase concentration
    P_over_P0: float            # relative partial pressure

    # --- Extended diagnostics (default so legacy constructors still work) ---
    net_retention_time_min: float = float("nan")  # t_R - t0 (matched)
    temp_col_K: float = float("nan")     # per-injection column temperature
    flow_col_mL_min: float = float("nan")  # per-injection column-T flow
    j_factor: float = 1.0                # per-injection James–Martin factor
    conditions_source: str = ""          # "measured" / "method_target"
    asymmetry_factor: float = float("nan")
    peak_clipped: bool = False           # apex at/above the digitiser ceiling


@dataclass
class IsothermPoint:
    """A single point on the adsorption isotherm."""

    P_over_P0: float
    n_adsorbed_mmol_g: float    # amount adsorbed (mMol/g)
    V_N_mL: float               # net retention volume at this point
    concentration_mol_m3: float = float("nan")  # gas-phase concentration


@dataclass
class BETQCFlags:
    """Quality control flags for a BET SSA result.

    Each flag is a (bool, str) tuple: (triggered, description).
    The ``flags`` property returns the list of triggered flag names.
    """

    # --- Fit quality ---
    few_points: bool = False        # < 5 points in BET window
    low_r2: bool = False            # R² < 0.99

    # --- BET constant ---
    c_below_1_5: bool = False       # C < 1.5 → approaching Type III, BET inapplicable
    c_below_2: bool = False         # C < 2 → BET model poorly conditioned
    c_above_100: bool = False       # C > 100 → unusually strong adsorption

    # --- Measurement sensitivity ---
    low_vn: bool = False            # min V_N < 1.0 mL
    very_low_vn: bool = False       # min V_N < 0.5 mL

    # --- Rouquerol criteria (IUPAC 2015) ---
    rouquerol_n_increasing: bool = False   # n(1 - P/P0) not monotonically increasing
    rouquerol_nm_outside_range: bool = False  # n_m falls outside measured P/P0 range

    # --- Isotherm shape ---
    isotherm_non_monotonic: bool = False  # adsorption isotherm decreases (unphysical)

    # --- V_N trend ---
    vn_increasing_with_concentration: bool = False  # V_N increases with c (should decrease in Henry regime)

    # --- Mass sensitivity ---
    mass_sensitive: bool = False    # ±5% mass change shifts SSA by > 10%

    # --- Net retention / dead-time (probe-independent) ---
    negative_net_retention: bool = False  # one or more injections with t_R ≤ t0
    com_before_methane: bool = False      # a probe CoM earlier than methane CoM

    # --- Methane dead-time marker stability ---
    methane_unstable: bool = False        # high methane dead-time SD or drift
    low_signal_to_methane: bool = False   # median probe net-ret ≈ methane noise

    # --- Detector saturation ---
    peak_saturation: bool = False         # one or more peaks clipped at ceiling

    # --- Sensitivity analyses ---
    retention_convention_sensitive: bool = False  # peak-max vs CoM SSA differ
    origin_sensitive: bool = False        # legacy vs corrected-origin SSA differ

    # Descriptive messages for triggered flags
    messages: list[str] = field(default_factory=list)

    @property
    def flags(self) -> list[str]:
        """Return list of triggered flag names."""
        triggered = []
        if self.few_points:
            triggered.append("FEW_POINTS")
        if self.low_r2:
            triggered.append("LOW_R2")
        if self.c_below_1_5:
            triggered.append("C<1.5_CRITICAL")
        elif self.c_below_2:
            triggered.append("C<2_WARN")
        if self.c_above_100:
            triggered.append("C>100_HIGH")
        if self.very_low_vn:
            triggered.append("VERY_LOW_VN")
        elif self.low_vn:
            triggered.append("LOW_VN")
        if self.rouquerol_n_increasing:
            triggered.append("ROUQUEROL_N_INCR")
        if self.rouquerol_nm_outside_range:
            triggered.append("ROUQUEROL_NM_RANGE")
        if self.isotherm_non_monotonic:
            triggered.append("ISOTHERM_NON_MONO")
        if self.vn_increasing_with_concentration:
            triggered.append("VN_TREND_WRONG")
        if self.mass_sensitive:
            triggered.append("MASS_SENSITIVE")
        if self.negative_net_retention:
            triggered.append("NEG_NET_RETENTION")
        if self.com_before_methane:
            triggered.append("COM_BEFORE_METHANE")
        if self.methane_unstable:
            triggered.append("METHANE_UNSTABLE")
        if self.low_signal_to_methane:
            triggered.append("LOW_SIGNAL_TO_METHANE")
        if self.peak_saturation:
            triggered.append("PEAK_SATURATION")
        if self.retention_convention_sensitive:
            triggered.append("RETENTION_SENSITIVE")
        if self.origin_sensitive:
            triggered.append("ORIGIN_SENSITIVE")
        return triggered

    @property
    def passed(self) -> bool:
        """True if no flags were triggered."""
        return len(self.flags) == 0

    @property
    def flag_string(self) -> str:
        """Comma-separated flag names, or 'OK'."""
        return ",".join(self.flags) if self.flags else "OK"


@dataclass
class MethaneStats:
    """Statistics of the methane dead-time markers across an experiment."""

    n: int = 0
    mean_max_min: float = float("nan")   # mean peak-max dead time (min)
    mean_cofm_min: float = float("nan")  # mean center-of-mass dead time (min)
    sd_max_min: float = float("nan")     # SD of peak-max dead time (min)
    range_max_min: float = float("nan")  # max−min of peak-max dead time (min)
    drift_max_min: float = float("nan")  # first→last change in peak-max (min)


@dataclass
class BETDiagnostics:
    """Probe-independent diagnostic scalars for a BET result.

    These are always computed (unlike the boolean QC flags, which only trigger
    on problems) so a report can show the measured envelope even for a clean
    dataset.
    """

    pp0_min: float = float("nan")
    pp0_max: float = float("nan")
    n_points_in_window: int = 0          # injections inside [p0_min, p0_max]
    n_negative_net_retention: int = 0    # injections with t_R ≤ t0
    n_clipped_peaks: int = 0
    median_net_retention_min: float = float("nan")
    signal_to_methane_ratio: float = float("nan")  # median net-ret / methane SD
    methane: MethaneStats = field(default_factory=MethaneStats)

    # Sensitivity analyses (populated when alternative fits are computed)
    ssa_alt_retention: float = float("nan")   # SSA under the other retention mode
    alt_retention_mode: str = ""
    ssa_alt_origin: float = float("nan")      # SSA under the corrected origin
    alt_origin_strategy: str = ""


@dataclass
class BETResult:
    """Results of BET analysis."""

    ssa_m2_g: float             # BET specific surface area (m²/g)
    n_monolayer_mmol_g: float   # monolayer capacity (mMol/g)
    C_bet: float                # BET constant (dimensionless)
    r_squared: float            # R² of the BET linear fit
    slope: float                # slope of BET linear plot
    intercept: float            # intercept of BET linear plot
    n_points: int               # number of points in the BET fit
    p_over_p0_range: tuple[float, float]  # actual P/P0 range used

    # Full data for diagnostics / plotting
    injections: list[InjectionResult] = field(default_factory=list)
    isotherm: list[IsothermPoint] = field(default_factory=list)
    bet_x: np.ndarray = field(default_factory=lambda: np.array([]))
    bet_y: np.ndarray = field(default_factory=lambda: np.array([]))

    # Experimental conditions
    temperature_K: float = 0.0
    sample_mass_mg: float = 0.0
    sample_name: str = ""
    james_martin_j: float = 1.0   # representative (mean) J applied to V_N

    # --- Probe identity and physical parameters used ---
    probe: str = "OCTANE"
    cross_section_m2: float = 6.3e-19   # octane default (OCTANE_CROSS_SECTION_M2)
    antoine: dict = field(default_factory=dict)
    p_sat_Pa: float = float("nan")
    probe_selection: ProbeSelection | None = None

    # --- Method / convention provenance ---
    retention_mode: str = "peak_max"
    concentration_mode: str = "eluted"
    origin_strategy: str = "legacy"
    conditions_source: str = ""   # "measured" / "method_target" / "mixed"
    flow_sccm: float = float("nan")

    # QC (populated by bet_quality_checks)
    qc: BETQCFlags | None = None

    # Diagnostics populated by the calling workflow
    diagnostics: BETDiagnostics | None = None

    # Isotherm classification (populated by classify_isotherm)
    classification: "IsothermClassification | None" = None


def _pct_delta(a: float, b: float | None) -> float:
    """Percent difference |a−b|/|b|·100, or NaN if not computable."""
    if b is None or b == 0 or math.isnan(a) or math.isnan(b):
        return float("nan")
    return abs(a - b) / abs(b) * 100.0


@dataclass
class IsothermClassification:
    """Classification of the adsorption isotherm shape and a BET-applicability
    verdict.

    BET surface area is only defined for a Type II isotherm — one with a
    monolayer "knee".  The two discriminants used here are both already
    computed by the standard pipeline:

    1. **BET constant C** ≈ exp[(E1 − E_L)/RT], the Boltzmann contrast between
       first-layer adsorption (probe-on-surface) and liquefaction
       (probe-on-probe).  C ≥ 2 → surface wins, Type II.  C < 1 → the probe
       binds itself better than the surface, Type III, and the monolayer is
       mathematically undefined.

    2. **V_N trend** — net retention volume is proportional to the local
       isotherm slope dq/dc.  Declining-then-flat V_N is the Type II monolayer
       signature; monotonically *rising* V_N is the Type III (cooperative)
       signature.

    See the theory note "Isotherm Shape as a Probe-Surface Diagnostic".
    """

    isotherm_type: str          # "II", "III", "II/III borderline", or "indeterminate"
    bet_applicable: bool        # whether a BET SSA should be trusted/reported
    vn_trend_slope: float       # slope of V_N vs P/P0 in the low-coverage region (mL per unit P/P0)
    vn_fractional_rise: float   # slope × range / mean(V_N): fractional change of V_N across the region
    vn_trend_r2: float          # R² of that linear trend
    vn_rising: bool             # True if V_N rises meaningfully (Type III tell)
    rationale: str              # one-line physical reason for the verdict
    recommendation: str         # what to do given the verdict


def classify_isotherm(
    result: BETResult,
    pp0_cap: float = 0.35,
    rise_frac_threshold: float = 0.20,
    rise_r2_threshold: float = 0.5,
) -> IsothermClassification:
    """Classify an isotherm as Type II / III / borderline and judge whether
    BET is applicable.

    Parameters
    ----------
    result : BETResult
        A completed BET analysis (needs ``C_bet`` and ``injections``).
    pp0_cap : float
        Upper P/P0 bound for the V_N-trend regression (default 0.35, the
        standard BET ceiling).  The Type II vs III distinction is cleanest in
        this low-coverage region, before high-P/P0 multilayer/turnover.
    rise_frac_threshold : float
        Minimum fractional rise of V_N across the region to call it "rising".
    rise_r2_threshold : float
        Minimum R² of the V_N-vs-P/P0 fit for the trend to count as real
        (vs noise).

    Returns
    -------
    IsothermClassification
    """
    C = result.C_bet

    # --- V_N trend in the low-coverage region (the dq/dc proxy) ---
    inj = [i for i in result.injections
           if 0 < i.P_over_P0 <= pp0_cap and not math.isnan(i.V_N_mL)]
    slope = frac = r2 = float("nan")
    rising = False
    if len(inj) >= 4:
        x = np.array([i.P_over_P0 for i in inj])
        y = np.array([i.V_N_mL for i in inj])
        coeffs = np.polyfit(x, y, 1)
        slope = float(coeffs[0])
        y_pred = np.polyval(coeffs, x)
        ss_res = float(np.sum((y - y_pred) ** 2))
        ss_tot = float(np.sum((y - y.mean()) ** 2))
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 1e-15 else 0.0
        mean_vn = float(y.mean())
        frac = slope * (x.max() - x.min()) / mean_vn if mean_vn > 0 else 0.0
        rising = bool((slope > 0) and (r2 >= rise_r2_threshold)
                      and (frac >= rise_frac_threshold))

    # --- Decision ---
    if math.isnan(C):
        iso_type, applicable = "indeterminate", False
        rationale = "BET constant C is undefined — fit did not converge"
        rec = "Inspect the isotherm and BET linearization plots manually."
    elif C < 1.0 or (C < 2.0 and rising):
        iso_type, applicable = "III", False
        if C < 1.0:
            rationale = (f"C = {C:.2f} < 1 — the probe binds itself better than "
                         f"the surface; no monolayer exists")
        else:
            rationale = (f"1 ≤ C = {C:.2f} < 2 with V_N rising "
                         f"({frac*100:.0f}% across the BET range) — cooperative "
                         f"adsorption, no monolayer knee")
        rec = ("BET inapplicable (Type III). Do not report this SSA. Use the "
               "Henry constant / dispersive γ_d from this run, and obtain SSA "
               "from cryogenic Kr/N₂ physisorption (then pass via bet_ssa_override).")
    elif C >= 2.0 and not rising:
        iso_type, applicable = "II", True
        rationale = f"C = {C:.2f} ≥ 2 and V_N declines/plateaus — clear monolayer"
        rec = "BET applicable."
    else:
        # 1 ≤ C < 2 without a rising V_N trend, or C ≥ 2 with conflicting rise
        iso_type, applicable = "II/III borderline", False
        rationale = (f"C = {C:.2f} (weak monolayer contrast); "
                     f"V_N {'rising' if rising else 'roughly flat'}")
        rec = ("Monolayer poorly defined — not reportable as a BET SSA. Treat "
               "the value as a rough estimate only and confirm with an "
               "orthogonal SSA method (e.g. cryogenic Kr/N₂ physisorption).")

    # --- Reportability gate ---------------------------------------------------
    # ``isotherm_type`` above is a descriptive shape label.  ``bet_applicable``
    # is the *acceptance* verdict and must reflect the complete gate set, not
    # shape alone: a reportable BET SSA needs a genuine Type II monolayer AND a
    # good linear fit AND resolvable retention AND a physical isotherm.  Any
    # failure makes the result non-reportable even if a number came out.
    gate_failures: list[str] = []
    if iso_type != "II":
        gate_failures.append(f"isotherm {iso_type}")
    if not math.isnan(result.r_squared) and result.r_squared < 0.99:
        gate_failures.append(f"R²={result.r_squared:.3f} < 0.99")
    qc = result.qc
    if qc is not None:
        if qc.very_low_vn:
            gate_failures.append("very low / nonpositive V_N")
        if qc.isotherm_non_monotonic:
            gate_failures.append("non-monotonic isotherm")
        if qc.vn_increasing_with_concentration:
            gate_failures.append("V_N rising in BET window")

    if applicable and gate_failures:
        applicable = False
        rationale = (f"{rationale}; acceptance gates failed: "
                     f"{', '.join(gate_failures)}")
        rec = ("Not reportable — BET acceptance gates failed "
               f"({', '.join(gate_failures)}). " + rec)

    return IsothermClassification(
        isotherm_type=iso_type,
        bet_applicable=applicable,
        vn_trend_slope=slope,
        vn_fractional_rise=frac,
        vn_trend_r2=r2,
        vn_rising=rising,
        rationale=rationale,
        recommendation=rec,
    )


# ---------------------------------------------------------------------------
# Physical chemistry calculations
# ---------------------------------------------------------------------------

def saturation_pressure(temperature_K: float,
                        C1: float, C2: float, C3: float,
                        C4: float, C5: float) -> float:
    """Compute saturation vapor pressure from extended Antoine equation.

    .. math::

        \\ln(P_{sat}) = C_1 + C_2/T + C_3 \\ln(T) + C_4 T^{C_5}

    Parameters
    ----------
    temperature_K : float
        Temperature in Kelvin.
    C1, C2, C3, C4, C5 : float
        Extended vapour-pressure coefficients in the stated equation.

    Returns
    -------
    float
        Saturation pressure in Pascals.
    """
    ln_psat = C1 + C2 / temperature_K + C3 * math.log(temperature_K) \
              + C4 * temperature_K ** C5
    return math.exp(ln_psat)


# Legacy octane convenience constants. Release-facing workflows should use
# source-attributed properties carried in the neutral input bundle.
OCTANE_ANTOINE = {
    "C1": 96.084,
    "C2": -7900.2,
    "C3": -11.003,
    "C4": 7.18e-6,
    "C5": 2.0,
}

OCTANE_CROSS_SECTION_M2 = 6.3e-19  # m²


def partial_pressure_ratio(n_injected_mol: float,
                           V_loop_m3: float,
                           temperature_K: float,
                           P_sat_Pa: float) -> float:
    """Compute P/P0 from the amount of probe injected.

    The gas-phase concentration of the probe in the injection plug:

    .. math::

        c = n_{inj} / V_{loop}

    The partial pressure via the ideal gas law:

    .. math::

        P = c \\cdot R \\cdot T

    Parameters
    ----------
    n_injected_mol : float
        Moles of probe injected from the declared calibration.
    V_loop_m3 : float
        Injection loop volume in m³.
    temperature_K : float
        Column temperature (K).
    P_sat_Pa : float
        Saturation vapor pressure of the probe at column temperature.

    Returns
    -------
    float
        P/P0 (dimensionless relative partial pressure).
    """
    c = n_injected_mol / V_loop_m3  # mol/m³
    P = c * R_GAS * temperature_K   # Pa
    return P / P_sat_Pa


def eluted_peak_concentration(
    apex_signal_uV: float,
    n_injected_mol: float,
    peak_area_uV_min: float,
    F_col_m3_min: float,
) -> float:
    """Gas-phase probe concentration at the peak apex in the carrier stream.

    The FID signal is proportional to the molar flow of probe eluting, so the
    calibration constant is ``κ = area / n_inj`` (µV per mol/min).  The molar
    flow at the apex is ``apex / κ = apex · n_inj / area`` and the gas-phase
    concentration is that molar flow divided by the volumetric carrier flow:

    .. math::

        c_{apex} = \\frac{S_{apex} \\cdot n_{inj} / A}{F_{col}}

    This is the concentration the isotherm relation ``V_N/m = dq/dc`` actually
    refers to for the peak-maximum method.  It is typically several-fold lower
    than the pre-injection loop concentration ``n_inj / V_loop``, and the gap
    widens for broader / more strongly retained peaks.

    Parameters
    ----------
    apex_signal_uV : float
        Baseline-subtracted signal at the peak maximum (µV).
    n_injected_mol : float
        Moles of probe injected from the declared calibration.
    peak_area_uV_min : float
        Integrated peak area (µV·min).
    F_col_m3_min : float
        Volumetric carrier flow at column temperature (m³/min).

    Returns
    -------
    float
        Gas-phase concentration at the peak apex (mol/m³), or 0.0 if the
        peak area or flow is non-positive.
    """
    if peak_area_uV_min <= 0 or F_col_m3_min <= 0:
        return 0.0
    molar_flow_mol_min = apex_signal_uV * (n_injected_mol / peak_area_uV_min)
    return molar_flow_mol_min / F_col_m3_min


def james_martin_j(p_inlet: float, p_outlet: float) -> float:
    """James–Martin gas compressibility (pressure-gradient) correction.

    The carrier gas expands along the column as the pressure falls from the
    inlet to the outlet, so the local flow — and hence the retention volume —
    must be corrected for the pressure gradient:

    .. math::

        j = \\frac{3}{2}
            \\frac{(p_i/p_o)^2 - 1}{(p_i/p_o)^3 - 1}

    (James & Martin, 1952). ``j`` is 1 at zero
    pressure drop and decreases toward ~2/3 as the drop grows; the net
    retention volume used for the isotherm is ``V_N = j · F · (t_R − t_0)``.

    Parameters
    ----------
    p_inlet, p_outlet : float
        Column inlet and outlet **absolute** pressures, in any consistent
        unit (only their ratio matters).

    Returns
    -------
    float
        The dimensionless correction factor ``j`` (1.0 if either pressure is
        non-positive or the ratio is ~1).
    """
    if p_inlet <= 0 or p_outlet <= 0:
        return 1.0
    r = p_inlet / p_outlet
    denom = r ** 3 - 1.0
    if abs(denom) < 1e-12:          # zero pressure drop → no correction
        return 1.0
    return 1.5 * (r ** 2 - 1.0) / denom


ORIGIN_STRATEGIES = ("legacy", "rectangular", "linear")


def _origin_vn_at_zero(
    concentrations_mol_m3: np.ndarray,
    V_N_m3: np.ndarray,
    strategy: str,
) -> float:
    """Estimate the net retention volume extrapolated to zero concentration.

    Returns the effective ``V_N(c→0)`` used to close the integral from the
    origin to the first measured point.  Clamped at ``≥ 0`` so a downward
    extrapolation can never inject a negative (unphysical) adsorbed amount.

    - ``legacy``       → 0 (the origin interval is omitted entirely).
    - ``rectangular``  → ``V_N`` at the first measured point (constant/flat
      extrapolation, i.e. a rectangle of height ``V_N[0]``).
    - ``linear``       → linear extrapolation of ``V_N`` vs ``c`` through the
      first two measured points back to ``c = 0``; if that would be negative
      (V_N rising with c, as in Type III), it is clamped to 0.
    """
    n = len(concentrations_mol_m3)
    if strategy == "legacy" or n == 0:
        return 0.0
    if strategy == "rectangular" or n == 1:
        return max(float(V_N_m3[0]), 0.0)
    if strategy == "linear":
        c0, c1 = concentrations_mol_m3[0], concentrations_mol_m3[1]
        v0, v1 = V_N_m3[0], V_N_m3[1]
        if c1 == c0:
            return max(float(v0), 0.0)
        slope = (v1 - v0) / (c1 - c0)
        vn_zero = v0 - slope * c0          # intercept at c = 0
        return max(float(vn_zero), 0.0)    # never allow negative V_N
    raise ValueError(f"Unknown origin strategy {strategy!r}; "
                     f"expected one of {ORIGIN_STRATEGIES}")


def build_adsorption_isotherm(
    concentrations_mol_m3: np.ndarray,
    V_N_m3: np.ndarray,
    sample_mass_g: float,
    origin: str = "legacy",
) -> np.ndarray:
    """Construct adsorption isotherm from net retention volumes.

    In IGC at finite concentration, each injection probes the isotherm
    at a different gas-phase concentration.  The net retention volume
    gives the isotherm derivative:

    .. math::

        dq/dc = V_N / m

    The amount adsorbed at concentration *c_i* is obtained by cumulative
    trapezoidal integration:

    .. math::

        q(c_i) = \\int_0^{c_i} \\frac{V_N(c)}{m} \\, dc

    **Origin treatment.** The historical rectangular convention integrates from *zero*
    pressure, but the measured points start at the first nonzero
    concentration ``c_0``.  The ``origin`` argument selects how the
    ``[0, c_0]`` interval is handled:

    - ``"legacy"`` (default): assign ``q(c_0) = 0`` — the origin interval is
      omitted.  This is the historically validated behaviour.
    - ``"rectangular"``: add a rectangle ``q(c_0) = V_N(c_0)/m · c_0``
      (constant-``V_N`` extrapolation to the origin).
    - ``"linear"``: linearly extrapolate ``V_N`` to ``c = 0`` (clamped at
      ``≥ 0``) and add the trapezoid ``q(c_0) = ½[V_N(0)+V_N(c_0)]/m · c_0``.

    Only the first point differs between strategies; the cumulative integral
    over the measured points is identical.

    Parameters
    ----------
    concentrations_mol_m3 : np.ndarray
        Gas-phase concentrations, sorted ascending (mol/m³).
    V_N_m3 : np.ndarray
        Net retention volumes at each concentration (m³).
    sample_mass_g : float
        Sample mass (g).
    origin : str
        Origin strategy — one of :data:`ORIGIN_STRATEGIES`.

    Returns
    -------
    np.ndarray
        Amount adsorbed q at each concentration (mol/g).
    """
    if origin not in ORIGIN_STRATEGIES:
        raise ValueError(f"Unknown origin strategy {origin!r}; "
                         f"expected one of {ORIGIN_STRATEGIES}")
    n = len(concentrations_mol_m3)
    q = np.zeros(n)
    if n == 0:
        return q

    # Origin contribution for the first measured point.
    #   legacy      → 0 (origin interval omitted)
    #   rectangular → V_N(c0)/m · c0            (rectangle of height V_N[0])
    #   linear      → ½[V_N(0)+V_N(c0)]/m · c0  (trapezoid, V_N(0) clamped ≥ 0)
    c0 = concentrations_mol_m3[0]
    if origin == "legacy":
        q[0] = 0.0
    elif origin == "rectangular":
        q[0] = float(V_N_m3[0]) / sample_mass_g * c0
    else:  # linear
        vn_zero = _origin_vn_at_zero(concentrations_mol_m3, V_N_m3, origin)
        q[0] = 0.5 * (vn_zero + float(V_N_m3[0])) / sample_mass_g * c0

    for i in range(1, n):
        dc = concentrations_mol_m3[i] - concentrations_mol_m3[i - 1]
        dq_dc_avg = 0.5 * (V_N_m3[i - 1] + V_N_m3[i]) / sample_mass_g
        q[i] = q[i - 1] + dq_dc_avg * dc
    return q


def _fit_bet_window(
    P_over_P0: np.ndarray,
    n_adsorbed_mmol_g: np.ndarray,
    p0_min: float,
    p0_max: float,
    cross_section_m2: float = OCTANE_CROSS_SECTION_M2,
) -> BETResult:
    """Fit the BET equation in a single P/P0 window.

    This is the inner workhorse — no adaptive range logic.  ``cross_section_m2``
    is the probe's molecular cross-sectional area used for the SSA conversion.
    """
    mask = (P_over_P0 >= p0_min) & (P_over_P0 <= p0_max) \
           & (n_adsorbed_mmol_g > 0)
    x = P_over_P0[mask]
    n = n_adsorbed_mmol_g[mask]

    if len(x) < 2:
        return BETResult(
            ssa_m2_g=float("nan"),
            n_monolayer_mmol_g=float("nan"),
            C_bet=float("nan"),
            r_squared=float("nan"),
            slope=float("nan"),
            intercept=float("nan"),
            n_points=len(x),
            p_over_p0_range=(p0_min, p0_max),
            bet_x=x,
            bet_y=np.array([]),
            cross_section_m2=cross_section_m2,
        )

    # BET y-axis: P/P0 / [n * (1 - P/P0)]   units: g/mMol
    y = x / (n * (1 - x))

    # Linear regression
    coeffs = np.polyfit(x, y, 1)
    slope, intercept = float(coeffs[0]), float(coeffs[1])

    # Extract BET parameters
    n_m = 1.0 / (slope + intercept)        # mMol/g
    C = 1.0 + slope / intercept if abs(intercept) > 1e-15 else float("nan")

    # SSA (probe-specific molecular cross-section)
    ssa = n_m * 1e-3 * N_AVOGADRO * cross_section_m2  # m²/g

    # R²
    y_pred = np.polyval(coeffs, x)
    ss_res = np.sum((y - y_pred) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    r_sq = 1 - ss_res / ss_tot if ss_tot > 1e-15 else float("nan")

    return BETResult(
        ssa_m2_g=ssa,
        n_monolayer_mmol_g=n_m,
        C_bet=C,
        r_squared=r_sq,
        slope=slope,
        intercept=intercept,
        n_points=len(x),
        p_over_p0_range=(float(x[0]), float(x[-1])),
        bet_x=x,
        bet_y=y,
        cross_section_m2=cross_section_m2,
    )


def _rouquerol_criterion_1_satisfied(
    bet_x: np.ndarray,
    iso_pp0: np.ndarray,
    iso_n: np.ndarray,
) -> bool:
    """Check Rouquerol criterion 1: n(1-P/P0) monotonically increasing.

    This is the criterion used for adaptive range narrowing, because
    removing high-P/P0 points can fix it.  Criterion 2 (monolayer P/P0
    in range) depends on C, which is a material property — narrowing
    the range cannot fix it for low-C materials.
    """
    if len(bet_x) < 3:
        return False

    bet_n = np.interp(bet_x, iso_pp0, iso_n)
    consistency = bet_n * (1 - bet_x)
    return bool(np.all(np.diff(consistency) >= -1e-12))


def bet_linearization(
    P_over_P0: np.ndarray,
    n_adsorbed_mmol_g: np.ndarray,
    p0_min: float = 0.05,
    p0_max: float = 0.35,
    adaptive: bool = True,
    cross_section_m2: float = OCTANE_CROSS_SECTION_M2,
) -> BETResult:
    """Perform BET linear regression on isotherm data.

    The linearized BET equation:

    .. math::

        \\frac{P/P_0}{n (1 - P/P_0)} = \\frac{1}{n_m C}
        + \\frac{C - 1}{n_m C} \\cdot \\frac{P}{P_0}

    where *n* is the amount adsorbed (mMol/g), *n_m* is the monolayer
    capacity, and *C* is the BET constant.

    When ``adaptive=True`` (default), the function first tries the full
    P/P0 window.  If the Rouquerol consistency criteria are not satisfied
    (n(1-P/P0) not increasing, or monolayer P/P0 outside range), it
    progressively narrows the window by removing the highest-P/P0 point
    until either the criteria pass or fewer than 3 points remain.  This
    follows the IUPAC 2015 recommendation for materials with low C.

    Parameters
    ----------
    P_over_P0 : np.ndarray
        Relative partial pressures.
    n_adsorbed_mmol_g : np.ndarray
        Amount adsorbed (mMol/g).
    p0_min, p0_max : float
        Initial P/P0 window for the BET fit.
    adaptive : bool
        If True, narrow the window to satisfy Rouquerol criteria.

    Returns
    -------
    BETResult
        Contains SSA, n_m, C, R², and diagnostic data.
    """
    # First fit with the full requested window
    result = _fit_bet_window(P_over_P0, n_adsorbed_mmol_g, p0_min, p0_max,
                             cross_section_m2)

    if not adaptive or result.n_points < 3:
        return result

    # Check Rouquerol criteria on the initial fit
    iso_pp0 = P_over_P0[P_over_P0 > 0]
    iso_n = n_adsorbed_mmol_g[P_over_P0 > 0]
    sort_idx = np.argsort(iso_pp0)
    iso_pp0 = iso_pp0[sort_idx]
    iso_n = iso_n[sort_idx]

    # Check Rouquerol criterion 1 only (n(1-P/P0) increasing).
    # Criterion 2 (monolayer P/P0 in range) depends on C and cannot
    # be fixed by narrowing — it's reported as a QC flag instead.
    if _rouquerol_criterion_1_satisfied(result.bet_x, iso_pp0, iso_n):
        return result

    # Adaptive narrowing: remove highest P/P0 points iteratively
    # until criterion 1 is satisfied or we hit the minimum point count.
    mask = (P_over_P0 >= p0_min) & (P_over_P0 <= p0_max) \
           & (n_adsorbed_mmol_g > 0)
    available_pp0 = np.sort(P_over_P0[mask])

    for i in range(1, len(available_pp0) - 2):
        # New upper bound: just below the i-th highest point
        new_max = float(available_pp0[-(i + 1)]) + 1e-10
        candidate = _fit_bet_window(P_over_P0, n_adsorbed_mmol_g,
                                     p0_min, new_max, cross_section_m2)

        if candidate.n_points < 3:
            break

        if _rouquerol_criterion_1_satisfied(candidate.bet_x, iso_pp0,
                                             iso_n):
            return candidate

    # Narrowing didn't help — return the original full-range fit
    return result


# ---------------------------------------------------------------------------
# Quality control
# ---------------------------------------------------------------------------

def bet_quality_checks(result: BETResult) -> BETQCFlags:
    """Run comprehensive quality checks on a BET result.

    Implements five categories of checks:

    1. **Fit quality** — number of points and R² of the BET linear fit.
    2. **BET constant** — C must be > 1 for a physically meaningful Type II
       isotherm; C > 100 is unusually high for physisorption.
    3. **Measurement sensitivity** — low V_N (< 1 mL) means net retention
       is only a few seconds, making peak detection imprecise.
    4. **Rouquerol criteria** (IUPAC 2015) — the standard consistency
       checks for validating that a BET fit is physically meaningful:
       (a) n(1 − P/P0) must be monotonically increasing in the fitted
       range; (b) the monolayer loading n_m should correspond to a P/P0
       that falls within the fitted range.
    5. **Isotherm shape** — the adsorption isotherm must be monotonically
       increasing (thermodynamic requirement for spontaneous adsorption).
    6. **V_N trend** — net retention volume should generally decrease with
       increasing injection concentration (transition from Henry's law to
       saturation).  An *increasing* V_N suggests an experimental problem.
    7. **Mass sensitivity** — flags cases where a ±5% error in sample
       mass would shift SSA by more than 10%.

    Parameters
    ----------
    result : BETResult
        A completed BET analysis result.

    Returns
    -------
    BETQCFlags
        Quality flags with descriptions.
    """
    qc = BETQCFlags()

    # --- 1. Fit quality ---
    if result.n_points < 5:
        qc.few_points = True
        qc.messages.append(
            f"Only {result.n_points} points in BET window (need ≥ 5 for "
            f"reliable fit)"
        )

    if not math.isnan(result.r_squared) and result.r_squared < 0.99:
        qc.low_r2 = True
        qc.messages.append(
            f"R² = {result.r_squared:.4f} < 0.99 — poor BET linearity"
        )

    # --- 2. BET constant quality tiers ---
    # C reflects the relative strength of adsorbate-surface vs adsorbate-
    # adsorbate interactions.  For Type II isotherms (BET applies), C > 2.
    # As C → 1, the isotherm approaches Type III (no monolayer inflection)
    # and the BET model becomes physically inappropriate.
    if not math.isnan(result.C_bet):
        if result.C_bet < 1.5:
            qc.c_below_1_5 = True
            qc.messages.append(
                f"C = {result.C_bet:.2f} < 1.5 — approaching Type III "
                f"isotherm; BET model is poorly applicable and SSA is "
                f"unreliable"
            )
        elif result.C_bet < 2.0:
            qc.c_below_2 = True
            qc.messages.append(
                f"C = {result.C_bet:.2f} < 2.0 — BET model is poorly "
                f"conditioned; monolayer capacity is poorly defined"
            )
        if result.C_bet > 100:
            qc.c_above_100 = True
            qc.messages.append(
                f"C = {result.C_bet:.2f} > 100 — unusually strong first-layer "
                f"interaction; may indicate chemisorption or fitting artifact"
            )

    # --- 3. Measurement sensitivity (V_N) ---
    if result.injections:
        min_vn = min(inj.V_N_mL for inj in result.injections)
        if min_vn < 0.5:
            qc.very_low_vn = True
            qc.messages.append(
                f"Min V_N = {min_vn:.3f} mL (< 0.5 mL) — very short "
                f"retention; BET result unreliable"
            )
        elif min_vn < 1.0:
            qc.low_vn = True
            qc.messages.append(
                f"Min V_N = {min_vn:.3f} mL (< 1.0 mL) — short retention "
                f"reduces peak detection precision"
            )

    # --- 4. Rouquerol criteria (IUPAC 2015 recommendations) ---
    if len(result.bet_x) >= 3 and len(result.isotherm) > 0:
        _check_rouquerol(result, qc)

    # --- 5. Isotherm shape ---
    if len(result.isotherm) >= 3:
        q_vals = [pt.n_adsorbed_mmol_g for pt in result.isotherm]
        for i in range(1, len(q_vals)):
            if q_vals[i] < q_vals[i - 1] - 1e-12:
                qc.isotherm_non_monotonic = True
                qc.messages.append(
                    f"Isotherm decreases at point {i}: "
                    f"q[{i-1}]={q_vals[i-1]:.6f} > q[{i}]={q_vals[i]:.6f} "
                    f"mMol/g — unphysical"
                )
                break

    # --- 6. V_N trend ---
    if len(result.injections) >= 4:
        _check_vn_trend(result, qc)

    # --- 7. Mass sensitivity ---
    if not math.isnan(result.ssa_m2_g) and result.sample_mass_mg > 0:
        _check_mass_sensitivity(result, qc)

    # --- 8. Net retention / dead-time integrity (probe-independent) ---
    _check_net_retention(result, qc)

    # --- 9. Detector saturation ---
    n_clipped = sum(1 for inj in result.injections if inj.peak_clipped)
    if n_clipped > 0:
        qc.peak_saturation = True
        qc.messages.append(
            f"{n_clipped} of {len(result.injections)} peaks appear clipped at "
            f"the digitiser ceiling — peak area and apex are underestimated"
        )

    # --- 10. Dead-time marker stability + signal-to-methane ---
    diag = result.diagnostics
    if diag is not None:
        m = diag.methane
        # Methane dead time should be highly reproducible.  Flag if its SD or
        # first→last drift exceeds ~2% of the mean dead time.
        if (m.n >= 2 and not math.isnan(m.mean_max_min) and m.mean_max_min > 0):
            rel_sd = m.sd_max_min / m.mean_max_min
            rel_drift = abs(m.drift_max_min) / m.mean_max_min
            if rel_sd > 0.02 or rel_drift > 0.02:
                qc.methane_unstable = True
                qc.messages.append(
                    f"Methane dead time unstable: SD={m.sd_max_min*60:.2f} s "
                    f"({rel_sd*100:.1f}%), drift={m.drift_max_min*60:+.2f} s "
                    f"({rel_drift*100:.1f}%) over {m.n} markers"
                )
        # If the probe's net retention is comparable to the methane timing
        # scatter, the isotherm rests on noise.
        if (not math.isnan(diag.signal_to_methane_ratio)
                and diag.signal_to_methane_ratio < 5.0):
            qc.low_signal_to_methane = True
            qc.messages.append(
                f"Median probe net retention is only "
                f"{diag.signal_to_methane_ratio:.1f}× the methane dead-time SD "
                f"— retention times poorly resolved from the dead-time marker"
            )

    return qc


def _check_net_retention(result: BETResult, qc: BETQCFlags) -> None:
    """Flag nonpositive net retention and probe CoM earlier than methane.

    A net retention time ``t_R − t0 ≤ 0`` means the probe eluted no later than
    the (unretained) methane marker — physically impossible for real
    adsorption and usually a sign the injection is below the detection floor.
    """
    if not result.injections:
        return
    n_neg = sum(1 for inj in result.injections
                if not math.isnan(inj.net_retention_time_min)
                and inj.net_retention_time_min <= 0)
    if n_neg > 0:
        qc.negative_net_retention = True
        qc.messages.append(
            f"{n_neg} of {len(result.injections)} injections have nonpositive "
            f"net retention (t_R ≤ t0) — excluded points; low-dose end of the "
            f"isotherm is unreliable"
        )
    # CoM earlier than methane CoM (only meaningful in CoM retention mode).
    diag = result.diagnostics
    if (result.retention_mode == "cofm" and diag is not None
            and not math.isnan(diag.methane.mean_cofm_min)):
        t0c = diag.methane.mean_cofm_min
        n_early = sum(1 for inj in result.injections
                      if inj.peak_cofm_time < t0c)
        if n_early > 0:
            qc.com_before_methane = True
            qc.messages.append(
                f"{n_early} probe center-of-mass times are earlier than the "
                f"methane CoM dead time ({t0c:.4f} min) — CoM unreliable"
            )


def _check_rouquerol(result: BETResult, qc: BETQCFlags) -> None:
    """Check Rouquerol consistency criteria for BET analysis.

    IUPAC 2015 recommendations (Thommes et al., Pure Appl. Chem. 87, 1051):

    Criterion 1: In the selected P/P0 range, the quantity n(1 - P/P0)
    must increase monotonically with P/P0.  If it decreases, the BET
    equation is being applied outside its valid range.

    Criterion 2: The monolayer capacity n_m should correspond to a P/P0
    value that falls within the fitted range.  Specifically, P/P0 at
    monolayer completion is approximately 1/(√C + 1).  If this falls
    outside the data range, the fit is extrapolating beyond the data.
    """
    # Build n_adsorbed array matching the BET x-values (P/P0 in fit range)
    # Match by finding isotherm points closest to each bet_x value
    iso_pp0 = np.array([pt.P_over_P0 for pt in result.isotherm])
    iso_n = np.array([pt.n_adsorbed_mmol_g for pt in result.isotherm])

    # For each BET x-value, find the closest isotherm point
    bet_n = np.interp(result.bet_x, iso_pp0, iso_n)

    # Criterion 1: n*(1 - P/P0) must be monotonically increasing
    consistency = bet_n * (1 - result.bet_x)
    diffs = np.diff(consistency)
    n_decreasing = int(np.sum(diffs < -1e-12))
    if n_decreasing > 0:
        qc.rouquerol_n_increasing = True
        qc.messages.append(
            f"Rouquerol criterion 1 violated: n·(1−P/P0) decreases at "
            f"{n_decreasing} of {len(diffs)} consecutive pairs in BET "
            f"window — P/P0 range may extend beyond BET validity"
        )

    # Criterion 2: P/P0 at monolayer completion ≈ 1/(√C + 1) should be
    # within the fitted range
    if not math.isnan(result.C_bet) and result.C_bet > 0:
        pp0_mono = 1.0 / (math.sqrt(result.C_bet) + 1.0)
        pp0_min_fit = float(result.bet_x[0])
        pp0_max_fit = float(result.bet_x[-1])
        if pp0_mono < pp0_min_fit or pp0_mono > pp0_max_fit:
            qc.rouquerol_nm_outside_range = True
            qc.messages.append(
                f"Rouquerol criterion 2: monolayer completion P/P0 = "
                f"{pp0_mono:.4f} (from C={result.C_bet:.2f}) falls "
                f"outside fitted range [{pp0_min_fit:.4f}, "
                f"{pp0_max_fit:.4f}]"
            )


def _check_vn_trend(result: BETResult, qc: BETQCFlags) -> None:
    """Check V_N trend within the BET fitting range.

    In the BET range, V_N can be roughly flat or gently declining —
    that's normal.  What matters is whether V_N is *increasing* through
    the BET range, which suggests multilayer onset is contaminating the
    fit.

    We compute a linear regression of V_N vs P/P0 for injections within
    the BET fitting range.  The flag triggers only if the slope is
    positive and the trend explains substantial variance (R² > 0.3),
    indicating a statistically meaningful increase rather than noise.
    """
    if len(result.bet_x) < 3:
        return

    pp0_min = float(result.bet_x[0])
    pp0_max = float(result.bet_x[-1])

    # Select injections within the BET fitting range
    bet_inj = [
        inj for inj in result.injections
        if pp0_min <= inj.P_over_P0 <= pp0_max
    ]
    if len(bet_inj) < 3:
        return

    pp0_vals = np.array([inj.P_over_P0 for inj in bet_inj])
    vn_vals = np.array([inj.V_N_mL for inj in bet_inj])

    # Linear regression: V_N = slope * P/P0 + intercept
    coeffs = np.polyfit(pp0_vals, vn_vals, 1)
    slope = float(coeffs[0])

    # R² of the linear fit
    vn_pred = np.polyval(coeffs, pp0_vals)
    ss_res = np.sum((vn_vals - vn_pred) ** 2)
    ss_tot = np.sum((vn_vals - vn_vals.mean()) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 1e-15 else 0.0

    # Flag if slope is positive AND explains meaningful variance
    if slope > 0 and r2 > 0.3:
        qc.vn_increasing_with_concentration = True
        qc.messages.append(
            f"V_N increases within BET range (P/P0 {pp0_min:.3f}–"
            f"{pp0_max:.3f}): slope = {slope:.3f} mL per unit P/P0, "
            f"R² = {r2:.3f} — possible multilayer onset in fit window"
        )


def _check_mass_sensitivity(result: BETResult, qc: BETQCFlags) -> None:
    """Flag if ±5% mass error would change SSA by > 10%.

    BET SSA depends on mass through the isotherm normalization (q = V_N/m).
    The effect is approximately linear: SSA ∝ 1/mass.  But the actual
    sensitivity depends on where the BET window falls relative to the
    isotherm — mass changes shift which P/P0 points map to which isotherm
    positions.

    For a purely linear (1/m) dependence, a 5% mass error gives a ~5%
    SSA change.  We flag cases where the effective sensitivity exceeds
    10%, which can happen when the BET window is near a transition in
    the isotherm shape.

    Since re-running the full pipeline at perturbed mass is expensive,
    we approximate using the 1/m scaling as a lower bound and flag based
    on additional indicators of sensitivity.
    """
    # For the linear approximation: delta_SSA/SSA ≈ delta_m/m
    # With 5% mass perturbation, SSA changes by ~5%.
    # That alone wouldn't trigger the 10% threshold.
    #
    # The real risk is when mass is very small (< 30 mg), where
    # weighing error is proportionally larger, AND when V_N is low
    # (compounding uncertainties).
    mass_mg = result.sample_mass_mg
    if mass_mg < 30.0:
        # At very low mass, a 1-2 mg weighing error (typical balance
        # precision) is > 5% of mass
        pct_from_1mg = 1.0 / mass_mg * 100
        if pct_from_1mg > 3.0:  # 1 mg is > 3% of mass
            qc.mass_sensitive = True
            qc.messages.append(
                f"Sample mass = {mass_mg:.1f} mg — a 1 mg weighing error "
                f"is {pct_from_1mg:.1f}% of mass, which propagates "
                f"directly to SSA"
            )


# ---------------------------------------------------------------------------
# Top-level BET analysis from raw exported data
# ---------------------------------------------------------------------------

def _is_clipped(signal: np.ndarray, min_run: int = 8) -> bool:
    """Heuristic flat-top (digitiser saturation) detector.

    A genuinely clipped peak rails at the ADC ceiling, producing a run of many
    *identical* maximum samples.  A normal peak apex — even a small one where
    quantisation puts a few samples near the top — never reaches an identical
    run of this length.  We therefore flag only when the longest run of
    samples exactly equal to the maximum (within a 1e-6 relative tolerance)
    reaches ``min_run`` (default 8), which keeps this free of false positives
    on ordinary or low-signal peaks while still catching true saturation.
    """
    if signal.size < min_run:
        return False
    smax = float(np.max(signal))
    if smax <= 0:
        return False
    at_ceiling = signal >= smax * (1.0 - 1e-6)
    best = run = 0
    for flag in at_ceiling:
        run = run + 1 if flag else 0
        best = max(best, run)
    return best >= min_run


@dataclass
class _InjectionRecord:
    """Raw per-injection data reused across retention/origin variants."""

    number: int
    target_coverage: float | None
    peak: dict
    n_injected_mol: float
    area: float
    temp_col_K: float
    flow_col_mL_min: float
    j_factor: float
    conditions_source: str
    clipped: bool


def _assemble_result(
    records: list[_InjectionRecord],
    *,
    mass_g: float,
    props: ProbeProperties,
    V_loop_m3: float,
    t0_max: float,
    t0_cofm: float,
    retention_mode: str,
    concentration_mode: str,
    origin: str,
    p0_min: float,
    p0_max: float,
) -> BETResult:
    """Build a :class:`BETResult` from processed injection records.

    Pure function of the already-extracted per-injection data plus the chosen
    conventions, so the main fit and the retention/origin sensitivity variants
    share exactly one implementation.
    """
    injections: list[InjectionResult] = []
    for rec in records:
        # Matched dead time: CoM retention uses the CoM methane dead time,
        # peak-max retention uses the peak-max methane dead time.  (The old
        # code subtracted the peak-max dead time from the probe CoM.)
        if retention_mode == "cofm":
            t0 = t0_cofm
            t_R = rec.peak["peak_cofm"]
        else:
            t0 = t0_max
            t_R = rec.peak["peak_max_time"]

        t_net = t_R - t0
        V_N_mL = rec.j_factor * t_net * rec.flow_col_mL_min

        P_sat = props.p_sat(rec.temp_col_K)
        if concentration_mode == "loop":
            c = rec.n_injected_mol / V_loop_m3
        else:
            c = eluted_peak_concentration(
                rec.peak["peak_max_value"], rec.n_injected_mol, rec.area,
                rec.flow_col_mL_min * 1e-6,
            )
        pp0 = c * R_GAS * rec.temp_col_K / P_sat

        injections.append(InjectionResult(
            injection_number=rec.number,
            target_coverage=rec.target_coverage,
            peak_max_time=rec.peak["peak_max_time"],
            peak_cofm_time=rec.peak["peak_cofm"],
            peak_area=rec.area,
            peak_height=rec.peak["peak_max_value"],
            n_injected_mol=rec.n_injected_mol,
            V_N_mL=V_N_mL,
            concentration_mol_m3=c,
            P_over_P0=pp0,
            net_retention_time_min=t_net,
            temp_col_K=rec.temp_col_K,
            flow_col_mL_min=rec.flow_col_mL_min,
            j_factor=rec.j_factor,
            conditions_source=rec.conditions_source,
            asymmetry_factor=rec.peak.get("asymmetry_factor", float("nan")),
            peak_clipped=rec.clipped,
        ))

    # Sort by concentration for isotherm integration.
    injections.sort(key=lambda inj: inj.concentration_mol_m3)

    # Only physically valid points (positive net retention → positive V_N)
    # enter the isotherm; nonpositive-retention injections are kept in the
    # result for QC but excluded from the cumulative integral.
    iso_inj = [inj for inj in injections if inj.V_N_mL > 0]

    concentrations = np.array([inj.concentration_mol_m3 for inj in iso_inj])
    V_N_m3_arr = np.array([inj.V_N_mL * 1e-6 for inj in iso_inj])
    pp0_arr = np.array([inj.P_over_P0 for inj in iso_inj])

    if len(iso_inj) >= 1:
        q_mmol_g = build_adsorption_isotherm(
            concentrations, V_N_m3_arr, mass_g, origin=origin) * 1000.0
    else:
        q_mmol_g = np.array([])

    isotherm = [
        IsothermPoint(P_over_P0=float(pp0), n_adsorbed_mmol_g=float(q),
                      V_N_mL=inj.V_N_mL, concentration_mol_m3=inj.concentration_mol_m3)
        for pp0, q, inj in zip(pp0_arr, q_mmol_g, iso_inj)
    ]

    if len(pp0_arr) >= 2:
        result = bet_linearization(pp0_arr, q_mmol_g, p0_min, p0_max,
                                   cross_section_m2=props.cross_section_m2)
    else:
        result = _fit_bet_window(pp0_arr, q_mmol_g, p0_min, p0_max,
                                 props.cross_section_m2)

    result.injections = injections
    result.isotherm = isotherm
    result.sample_mass_mg = mass_g * 1000.0
    result.probe = props.name
    result.cross_section_m2 = props.cross_section_m2
    result.antoine = props.antoine
    result.retention_mode = retention_mode
    result.concentration_mode = concentration_mode
    result.origin_strategy = origin
    # Representative (mean) temperature, flow and J across injections.
    if injections:
        result.temperature_K = float(np.mean([i.temp_col_K for i in injections]))
        result.james_martin_j = float(np.mean([i.j_factor for i in injections]))
        result.p_sat_Pa = props.p_sat(result.temperature_K)
        sources = {i.conditions_source for i in injections}
        result.conditions_source = (sources.pop() if len(sources) == 1
                                    else "mixed")
    return result
