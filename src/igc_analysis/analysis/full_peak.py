"""Full-peak nonlinear inverse chromatography for finite-concentration IGC.

Instead of reducing each injection to a peak-maximum or centre-of-mass
retention time, this module keeps the complete baseline-corrected FID trace and
fits a forward equilibrium-dispersive column model (see
:mod:`igc_analysis.analysis.column_model`) to **all** injections jointly, so an
adsorption isotherm is identified from peak shape.

Pipeline
--------
1. :func:`build_trace_dataset_from_neutral` — validated source-neutral,
   calibrated long-form traces (mass-conserving signal → molar flow → outlet
   concentration → P/P0).
2. :func:`characterize_blocks` — per-block methane transport (never pooled).
3. :func:`fit_model` — joint inverse fit with shared adsorption parameters.
4. :func:`compare_models` — AICc/BIC, held-out dose prediction, identifiability.
5. :func:`compute_ssa_if_identifiable` — the mandatory SSA guardrail.

Calibration, probe properties, chromatograms, and acquisition conditions are
read exclusively from a validated neutral bundle.

See ``docs/full_peak_architecture.md`` for the model equations and conventions.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from igc_analysis.constants import R_GAS
from igc_analysis.analysis.column_model import (
    ColumnGeometry, MethaneTransport, SolveResult, TransportParams,
    characterize_methane_transport, make_geometry, solve_column,
    peak_moments,
)
from igc_analysis.analysis.isotherm_models import (
    IsothermModel, get_model, is_cooperative,
)

# Default spatial discretisation.  Must be held FIXED between the methane
# calibration and the probe fit: the effective plate number absorbs the scheme's
# numerical dispersion, so changing n_cells changes what N means.
DEFAULT_N_CELLS = 160

# A parameter is called unidentifiable above this relative standard error.
IDENTIFIABILITY_RSE_LIMIT = 0.5

# Structural-confounding limits.  Both are computed on the *correlation* matrix
# rather than the raw Jacobian, so they are scale-free: a large J'J condition
# number usually only reflects parameters carrying different units (e.g. K_F
# ~1e-6 vs n ~1), which says nothing about identifiability.  Parameters
# correlated above CORRELATION_LIMIT — or a correlation matrix conditioned
# worse than CONDITION_LIMIT, which also catches 3-way collinearity that
# pairwise correlation misses — are not separately identifiable however small
# the residual happens to be.
CORRELATION_LIMIT = 0.99
CONDITION_LIMIT = 1e4


# ---------------------------------------------------------------------------
# Data containers
# ---------------------------------------------------------------------------

@dataclass
class Injection:
    """One calibrated probe injection ready for forward modelling."""

    block: str
    injection_number: int
    name: str
    probe: str
    time_min: np.ndarray
    signal_raw_uV: np.ndarray
    baseline_uV: np.ndarray
    signal_corrected_uV: np.ndarray   # may contain negative residuals
    n_injected_mol: float
    c_out_mol_m3: np.ndarray          # from the nonnegative-clipped signal
    pp0: np.ndarray
    flow_col_m3_min: float
    temp_col_K: float
    pressure_drop_torr: float
    fid_gain: float
    p_sat_Pa: float
    clipped: bool
    transport_mode: str = "fixed_block_mean"
    assigned_t0_min: float | None = None
    t0_interpolation_fraction: float | None = None
    pre_bracket_t0_min: float | None = None
    post_bracket_t0_min: float | None = None
    bracket_drift_min: float | None = None
    t0_assignment_basis: str = "block mean"
    source_methane_pre: tuple[str, ...] = ()
    source_methane_post: tuple[str, ...] = ()

    @property
    def peak_scale(self) -> float:
        """Amplitude scale used to normalise this injection's residuals."""
        m = float(np.max(self.c_out_mol_m3)) if self.c_out_mol_m3.size else 0.0
        return m if m > 0 else 1.0


@dataclass
class BlockData:
    """One completed measurement block: methane markers + probe injections."""

    block: str
    transport: MethaneTransport
    geometry: ColumnGeometry
    injections: list[Injection] = field(default_factory=list)
    methane_traces: list[tuple[np.ndarray, np.ndarray]] = field(default_factory=list)
    methane_names: list[str] = field(default_factory=list)
    schedule: list[str] = field(default_factory=list)


@dataclass
class FitResult:
    """Outcome of a joint inverse fit for one isotherm model."""

    model_name: str
    params: np.ndarray
    param_names: list[str]
    param_units: list[str]
    std_errors: np.ndarray
    rel_std_errors: np.ndarray
    identifiable: dict[str, bool]
    condition_number: float
    max_abs_correlation: float
    converged: bool
    n_starts: int
    n_successful_starts: int
    rss: float
    n_points: int
    n_params: int
    rmse_normalised: float
    r_squared: float
    aicc: float
    bic: float
    mass_balance_mean: float
    mass_balance_min: float
    per_injection_rmse: dict[str, float]
    cooperative: bool
    message: str = ""

    @property
    def all_identifiable(self) -> bool:
        return all(self.identifiable.values()) if self.identifiable else True


# ---------------------------------------------------------------------------
# 1. Trace dataset construction
# ---------------------------------------------------------------------------

def _neutral_condition_mean(
    conditions: pd.DataFrame,
    injection_id: str,
    quantity: str,
    *,
    value_role: str = "measured",
) -> float | None:
    rows = conditions[
        (conditions["injection_id"] == injection_id)
        & (conditions["quantity"] == quantity)
        & (conditions["value_role"] == value_role)
    ]
    if rows.empty:
        return None
    values = pd.to_numeric(rows["value"], errors="raise").to_numpy(dtype=float)
    return float(np.mean(values))


def _neutral_trace(
    traces: pd.DataFrame,
    injection_id: str,
    detector_channel: str,
) -> tuple[np.ndarray, np.ndarray, str]:
    rows = traces[
        (traces["injection_id"] == injection_id)
        & (traces["detector_channel"] == detector_channel)
    ].sort_values("point_index")
    if rows.empty:
        raise ValueError(f"neutral injection {injection_id!r} has no declared trace")
    units = set(rows["signal_unit"].astype(str))
    if len(units) != 1:
        raise ValueError(f"neutral injection {injection_id!r} has inconsistent signal units")
    return (
        pd.to_numeric(rows["time_s"], errors="raise").to_numpy(dtype=float) / 60.0,
        pd.to_numeric(rows["signal_raw"], errors="raise").to_numpy(dtype=float),
        units.pop(),
    )


def _neutral_calibrated_amount(
    time_min: np.ndarray,
    signal_corrected: np.ndarray,
    signal_unit: str,
    calibration: pd.Series,
) -> float:
    area_unit = str(calibration["area_unit"])
    if area_unit == f"{signal_unit}_min":
        area = float(np.trapezoid(np.maximum(signal_corrected, 0.0), time_min))
    elif area_unit == f"{signal_unit}_s":
        area = float(np.trapezoid(np.maximum(signal_corrected, 0.0), time_min * 60.0))
    else:
        raise ValueError(
            f"calibration area_unit {area_unit!r} is incompatible with "
            f"trace signal_unit {signal_unit!r}"
        )
    if area <= 0:
        return 0.0

    p0 = float(calibration["parameter_0"])
    p1 = float(calibration["parameter_1"])
    model = str(calibration["calibration_model"])
    if model == "linear":
        amount = p0 + p1 * area
    elif model == "quadratic":
        amount = p0 + p1 * area + float(calibration["parameter_2"]) * area**2
    elif model == "power_law":
        amount = p0 * area**p1
    else:
        raise ValueError(f"unsupported neutral calibration model {model!r}")
    if not math.isfinite(amount) or amount < 0:
        raise ValueError("neutral calibration produced a nonphysical amount")
    return amount


def build_trace_dataset_from_neutral(
    bundle_dirs: dict[str, str | Path],
    probe_override: str = "auto",
    n_cells: int = DEFAULT_N_CELLS,
    transport_mode: str = "fixed_block_mean",
    verbose: bool = True,
) -> list[BlockData]:
    """Build full-peak inputs from validated ``igc-neutral-data/0.2.0`` bundles.

    This is the public source-neutral path. It never reads an extraction
    directory, source filename convention, embedded method file, or protected
    schema. Each mapping entry is one independently characterized block.
    """

    from igc_analysis.analysis.peak_detection import detect_baseline, subtract_baseline
    from igc_analysis.io.neutral_data import read_neutral_bundle

    valid_modes = {"fixed_block_mean", "bracket_interpolated",
                   "bracket_pre", "bracket_post"}
    if transport_mode not in valid_modes:
        raise ValueError(f"unknown transport_mode {transport_mode!r}; "
                         f"choose from {sorted(valid_modes)}")

    blocks: list[BlockData] = []
    for label, bundle_dir in bundle_dirs.items():
        bundle = read_neutral_bundle(bundle_dir)
        experiment = bundle.table("experiment.csv")
        injections_table = bundle.table("injections.csv").sort_values("sequence_index")
        components = bundle.table("injection_components.csv")
        properties = bundle.table("probe_properties.csv")
        calibrations = bundle.table("calibration.csv")
        conditions = bundle.table("conditions.csv")
        traces = bundle.table("traces.csv")

        if len(experiment) != 1:
            raise ValueError(f"{label}: full-peak requires exactly one experiment")
        mass_g = float(experiment.iloc[0]["sample_mass_g"])
        block_ids = set(injections_table["block_id"].astype(str))
        if len(block_ids) != 1:
            raise ValueError(f"{label}: one neutral bundle must contain exactly one block")

        component_properties = components.merge(properties, on="probe_id", validate="many_to_one")
        analytes = component_properties[
            component_properties["component_role"] == "analyte"
        ]
        available_probe_ids = set(analytes["probe_id"].astype(str))
        if not available_probe_ids:
            raise ValueError(f"{label}: neutral bundle contains no analyte components")
        if probe_override == "auto":
            if len(available_probe_ids) != 1:
                raise ValueError(
                    f"{label}: full-peak requires one analyte probe per bundle; "
                    f"found {sorted(available_probe_ids)}"
                )
            selected_probe_id = next(iter(available_probe_ids))
        else:
            requested = probe_override.casefold()
            matching = analytes[
                analytes["probe_name"].astype(str).str.casefold() == requested
            ]
            matches = set(matching["probe_id"].astype(str))
            if len(matches) != 1:
                raise ValueError(f"{label}: probe {probe_override!r} is not uniquely available")
            selected_probe_id = next(iter(matches))
        probe_name = str(
            properties.loc[properties["probe_id"] == selected_probe_id, "probe_name"].iloc[0]
        ).upper()

        selected_components = components[
            (components["probe_id"] == selected_probe_id)
            & (components["component_role"] == "analyte")
        ]
        probe_ids = set(selected_components["injection_id"].astype(str))
        probe_rows = injections_table[
            injections_table["injection_id"].isin(probe_ids)
            & (injections_table["role"] == "probe")
        ]
        methane_rows = injections_table[injections_table["role"] == "dead_time"]
        if probe_rows.empty or methane_rows.empty:
            raise ValueError(f"{label}: full-peak requires probe and dead-time traces")

        schedule = list(injections_table["injection_id"].astype(str))
        probe_names = list(probe_rows["injection_id"].astype(str))
        methane_names = list(methane_rows["injection_id"].astype(str))

        probe_flow: dict[str, float] = {}
        probe_temp: dict[str, float] = {}
        probe_pressure_drop_torr: dict[str, float] = {}
        for injection_id in probe_names:
            temp_K = _neutral_condition_mean(conditions, injection_id, "column_temperature")
            if temp_K is None:
                temp_K = _neutral_condition_mean(
                    conditions, injection_id, "column_temperature", value_role="target"
                )
            flow_m3_s = _neutral_condition_mean(conditions, injection_id, "flow_column")
            if flow_m3_s is None:
                standard_flow = _neutral_condition_mean(
                    conditions, injection_id, "flow_standard"
                )
                if standard_flow is None:
                    standard_flow = _neutral_condition_mean(
                        conditions, injection_id, "flow_standard", value_role="target"
                    )
                if standard_flow is not None and temp_K is not None:
                    flow_m3_s = standard_flow * temp_K / 273.15
            if temp_K is None or flow_m3_s is None or flow_m3_s <= 0:
                raise ValueError(f"{label}/{injection_id}: temperature or column flow is missing")
            pressure_drop_pa = _neutral_condition_mean(
                conditions, injection_id, "pressure_drop"
            )
            probe_temp[injection_id] = temp_K
            probe_flow[injection_id] = flow_m3_s * 60.0
            probe_pressure_drop_torr[injection_id] = (
                pressure_drop_pa / 133.32236842105263
                if pressure_drop_pa is not None else float("nan")
            )

        flow_col_m3_min = float(np.mean(list(probe_flow.values())))

        methane_traces: list[tuple[np.ndarray, np.ndarray]] = []
        rough_t0: float | None = None
        for injection_id in methane_names:
            row = methane_rows[methane_rows["injection_id"] == injection_id].iloc[0]
            time_min, signal_raw, _ = _neutral_trace(
                traces, injection_id, str(row["detector_channel"])
            )
            if rough_t0 is None:
                rough_t0 = float(time_min[int(np.argmax(signal_raw))])
            intercept, gradient = detect_baseline(
                time_min, signal_raw, dead_time_min=rough_t0
            )
            methane_traces.append(
                (time_min, subtract_baseline(time_min, signal_raw, intercept, gradient))
            )

        transport = characterize_methane_transport(
            label, methane_traces, flow_col_m3_min, mass_g, n_cells=n_cells
        )
        geometry = make_geometry(mass_g, flow_col_m3_min, transport.t0_min)
        methane_moments = {
            name: peak_moments(time, signal)[1]
            for name, (time, signal) in zip(methane_names, methane_traces)
        }
        probe_positions = {name: schedule.index(name) + 1 for name in probe_names}
        first_probe = min(probe_positions.values())
        last_probe = max(probe_positions.values())
        pre_names = [name for name in methane_names if schedule.index(name) + 1 < first_probe]
        post_names = [name for name in methane_names if schedule.index(name) + 1 > last_probe]
        if transport_mode != "fixed_block_mean" and (not pre_names or not post_names):
            raise ValueError(
                f"{label}: {transport_mode} requires dead-time markers before and after "
                f"the probe block; schedule={schedule!r}"
            )
        pre_t0 = float(np.mean([methane_moments[n] for n in pre_names])) if pre_names else None
        post_t0 = float(np.mean([methane_moments[n] for n in post_names])) if post_names else None
        if pre_names and post_names:
            pre_anchor = float(np.mean([schedule.index(n) + 1 for n in pre_names]))
            post_anchor = float(np.mean([schedule.index(n) + 1 for n in post_names]))
        else:
            pre_anchor, post_anchor = 0.0, float(len(probe_names) + 1)

        if verbose:
            print(
                f"  [{label}] {probe_name}: {len(probe_names)} probe + "
                f"{len(methane_names)} dead-time injections from neutral bundle "
                f"{bundle.dataset_id}; t0={transport.t0_min:.4f} min"
            )

        result_injections: list[Injection] = []
        for index, injection_id in enumerate(probe_names, start=1):
            injection_row = probe_rows[probe_rows["injection_id"] == injection_id].iloc[0]
            component_rows = selected_components[
                selected_components["injection_id"] == injection_id
            ]
            if len(component_rows) != 1:
                raise ValueError(f"{label}/{injection_id}: expected one analyte component")
            component = component_rows.iloc[0]
            calibration_id = str(component["calibration_id"])
            calibration_rows = calibrations[
                calibrations["calibration_id"] == calibration_id
            ]
            if len(calibration_rows) != 1:
                raise ValueError(f"{label}/{injection_id}: calibration is missing or ambiguous")
            calibration = calibration_rows.iloc[0]
            if str(calibration["probe_id"]) != selected_probe_id:
                raise ValueError(f"{label}/{injection_id}: calibration probe does not match")

            time_min, signal_raw, signal_unit = _neutral_trace(
                traces, injection_id, str(injection_row["detector_channel"])
            )
            intercept, gradient = detect_baseline(
                time_min, signal_raw, dead_time_min=transport.t0_min
            )
            baseline = intercept + gradient * time_min
            signal_corrected = signal_raw - baseline
            area_min = float(np.trapezoid(np.maximum(signal_corrected, 0.0), time_min))
            amount_mol = _neutral_calibrated_amount(
                time_min, signal_corrected, signal_unit, calibration
            )
            flow_m3_min = probe_flow[injection_id]
            signal_positive = np.maximum(signal_corrected, 0.0)
            if area_min > 0 and flow_m3_min > 0:
                molar_flow = amount_mol * signal_positive / area_min
                concentration = molar_flow / flow_m3_min
            else:
                concentration = np.zeros_like(time_min)
            p_sat = float(component["saturation_vapor_pressure_Pa"])
            if not math.isfinite(p_sat) or p_sat <= 0:
                raise ValueError(f"{label}/{injection_id}: saturation pressure is required")
            temperature_K = probe_temp[injection_id]
            pp0 = concentration * R_GAS * temperature_K / p_sat

            clipping_value = str(injection_row["clipping_observed"]).casefold()
            if clipping_value in {"true", "false"}:
                clipped = clipping_value == "true"
            else:
                signal_max = float(np.max(signal_raw))
                clipped = bool(np.sum(signal_raw >= signal_max * (1 - 1e-6)) >= 8)

            result = Injection(
                block=label,
                injection_number=index,
                name=injection_id,
                probe=probe_name,
                time_min=time_min,
                signal_raw_uV=signal_raw,
                baseline_uV=baseline,
                signal_corrected_uV=signal_corrected,
                n_injected_mol=amount_mol,
                c_out_mol_m3=concentration,
                pp0=pp0,
                flow_col_m3_min=flow_m3_min,
                temp_col_K=temperature_K,
                pressure_drop_torr=probe_pressure_drop_torr[injection_id],
                fid_gain=(
                    float(injection_row["detector_gain"])
                    if str(injection_row["detector_gain"]) != "" else 1.0
                ),
                p_sat_Pa=p_sat,
                clipped=clipped,
            )
            result.transport_mode = transport_mode
            result.pre_bracket_t0_min = pre_t0
            result.post_bracket_t0_min = post_t0
            result.bracket_drift_min = (
                post_t0 - pre_t0 if pre_t0 is not None and post_t0 is not None else None
            )
            result.source_methane_pre = tuple(pre_names)
            result.source_methane_post = tuple(post_names)
            if transport_mode == "fixed_block_mean":
                result.assigned_t0_min = transport.t0_min
                result.t0_assignment_basis = "mean of all block dead-time first moments"
            elif transport_mode == "bracket_pre":
                result.assigned_t0_min = pre_t0
                result.t0_interpolation_fraction = 0.0
                result.t0_assignment_basis = "pre-block dead-time mean"
            elif transport_mode == "bracket_post":
                result.assigned_t0_min = post_t0
                result.t0_interpolation_fraction = 1.0
                result.t0_assignment_basis = "post-block dead-time mean"
            else:
                position = probe_positions[injection_id]
                fraction = (position - pre_anchor) / (post_anchor - pre_anchor)
                result.t0_interpolation_fraction = float(fraction)
                result.assigned_t0_min = float(pre_t0 + fraction * (post_t0 - pre_t0))
                result.t0_assignment_basis = (
                    "neutral sequence-position interpolation; acquisition timestamps absent"
                )
            result_injections.append(result)

        blocks.append(
            BlockData(
                block=label,
                transport=transport,
                geometry=geometry,
                injections=result_injections,
                methane_traces=methane_traces,
                methane_names=methane_names,
                schedule=schedule,
            )
        )

    if not blocks:
        raise ValueError("No neutral bundle contained full-peak-ready traces.")
    return blocks


def traces_to_dataframe(blocks: list[BlockData]) -> pd.DataFrame:
    """Flatten all blocks into the long-form ``full_peak_traces`` table."""
    rows = []
    for blk in blocks:
        for inj in blk.injections:
            n = len(inj.time_min)
            rows.append(pd.DataFrame({
                "block": inj.block,
                "injection": inj.injection_number,
                "injection_name": inj.name,
                "probe": inj.probe,
                "time_min": inj.time_min,
                "signal_raw_uV": inj.signal_raw_uV,
                "baseline_uV": inj.baseline_uV,
                "signal_corrected_uV": inj.signal_corrected_uV,
                "response_density_per_min": (
                    np.maximum(inj.signal_corrected_uV, 0.0)
                    / max(float(np.trapezoid(np.maximum(inj.signal_corrected_uV, 0.0),
                                             inj.time_min)), 1e-30)),
                "n_injected_mol": inj.n_injected_mol,
                "molar_flow_mol_min": inj.c_out_mol_m3 * inj.flow_col_m3_min,
                "c_out_mol_m3": inj.c_out_mol_m3,
                "pp0": inj.pp0,
                "flow_col_m3_min": inj.flow_col_m3_min,
                "temp_col_K": inj.temp_col_K,
                "pressure_drop_torr": inj.pressure_drop_torr,
                "fid_gain": inj.fid_gain,
                "p_sat_Pa": inj.p_sat_Pa,
                "clipped": inj.clipped,
                "transport_mode": inj.transport_mode,
                "assigned_t0_min": inj.assigned_t0_min,
                "t0_interpolation_fraction": inj.t0_interpolation_fraction,
                "pre_bracket_t0_min": inj.pre_bracket_t0_min,
                "post_bracket_t0_min": inj.post_bracket_t0_min,
                "bracket_drift_s": (inj.bracket_drift_min * 60.0
                                    if inj.bracket_drift_min is not None else np.nan),
                "t0_assignment_basis": inj.t0_assignment_basis,
                "source_methane_pre": ";".join(inj.source_methane_pre),
                "source_methane_post": ";".join(inj.source_methane_post),
            }))
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def transport_to_dataframe(blocks: list[BlockData]) -> pd.DataFrame:
    """Per-block methane transport summary table."""
    return pd.DataFrame([{
        "block": b.transport.block,
        "n_methane_markers": b.transport.n_markers,
        "t0_min": b.transport.t0_min,
        "t0_sd_s": b.transport.t0_sd_min * 60.0,
        "t0_range_s": b.transport.t0_range_min * 60.0,
        "sigma_mean_min": b.transport.sigma_mean_min,
        "plate_number_effective": b.transport.plate_number,
        "plate_number_moment_sd": b.transport.plate_number_sd,
        "t_inj_effective_min": b.transport.t_inj_min,
        "void_volume_m3": b.transport.void_volume_m3,
        "void_volume_mL": b.transport.void_volume_m3 * 1e6,
        "sample_mass_g": b.geometry.sample_mass_g,
        "phase_ratio_g_per_m3": b.geometry.phase_ratio_g_m3,
        "flow_col_mL_min": b.geometry.flow_col_m3_min * 1e6,
        "methane_fit_rmse_normalised": b.transport.fit_rmse,
    } for b in blocks])


# ---------------------------------------------------------------------------
# 2. Forward prediction and residuals
# ---------------------------------------------------------------------------

def predict_injection(inj: Injection, blk: BlockData, model: IsothermModel,
                      params: np.ndarray, n_cells: int = DEFAULT_N_CELLS
                      ) -> SolveResult:
    """Forward-solve one injection with the block's fixed transport."""
    if inj.assigned_t0_min is None or inj.transport_mode == "fixed_block_mean":
        transport = blk.transport.to_params()
        geometry = blk.geometry
    else:
        transport = TransportParams(
            t0_min=inj.assigned_t0_min,
            t_inj_min=blk.transport.t_inj_min,
            plate_number=blk.transport.plate_number,
        )
        geometry = make_geometry(
            blk.geometry.sample_mass_g, inj.flow_col_m3_min,
            inj.assigned_t0_min)
    return solve_column(inj.time_min, inj.n_injected_mol, transport,
                        geometry, model, params, n_cells=n_cells)


def bracket_assignment_to_dataframe(blocks: list[BlockData]) -> pd.DataFrame:
    """Return the auditable per-injection methane-bracket assignment."""
    rows = []
    for blk in blocks:
        for inj in blk.injections:
            _, probe_m1, _ = peak_moments(inj.time_min, inj.c_out_mol_m3)
            rows.append({
                "block": inj.block,
                "injection": inj.injection_number,
                "injection_name": inj.name,
                "n_injected_mol": inj.n_injected_mol,
                "transport_mode": inj.transport_mode,
                "assigned_t0_min": inj.assigned_t0_min,
                "probe_first_moment_min": probe_m1,
                "net_first_moment_min": probe_m1 - float(inj.assigned_t0_min),
                "t0_interpolation_fraction": inj.t0_interpolation_fraction,
                "pre_bracket_t0_min": inj.pre_bracket_t0_min,
                "post_bracket_t0_min": inj.post_bracket_t0_min,
                "bracket_drift_s": (inj.bracket_drift_min * 60.0
                                    if inj.bracket_drift_min is not None else np.nan),
                "assignment_basis": inj.t0_assignment_basis,
                "source_methane_pre": ";".join(inj.source_methane_pre),
                "source_methane_post": ";".join(inj.source_methane_post),
            })
    return pd.DataFrame(rows)


def _residuals(blocks: list[BlockData], model: IsothermModel,
               params: np.ndarray, n_cells: int,
               exclude: tuple[str, int] | None = None) -> np.ndarray:
    """Amplitude-normalised residual vector across all injections.

    Each injection's residual is divided by its own peak scale so that
    high-dose peaks cannot dominate the objective purely through amplitude.
    """
    out = []
    for blk in blocks:
        for inj in blk.injections:
            if exclude is not None and (inj.block, inj.injection_number) == exclude:
                continue
            try:
                res = predict_injection(inj, blk, model, params, n_cells)
                pred = res.c_out
            except Exception:
                pred = np.full_like(inj.c_out_mol_m3, np.inf)
            out.append((inj.c_out_mol_m3 - pred) / inj.peak_scale)
    return np.concatenate(out) if out else np.array([0.0])


# ---------------------------------------------------------------------------
# 3. Joint inverse fitting
# ---------------------------------------------------------------------------

def fit_model(
    blocks: list[BlockData],
    model_name: str,
    n_starts: int = 4,
    n_cells: int = DEFAULT_N_CELLS,
    seed: int = 0,
    exclude: tuple[str, int] | None = None,
    x0: np.ndarray | None = None,
    verbose: bool = True,
) -> FitResult:
    """Jointly fit one isotherm model to every probe injection.

    A single shared adsorption parameter set is used for all injections (one
    packed column); per-injection dose and conditions and per-block transport
    are fixed inputs, not free parameters.

    Parameters
    ----------
    blocks : list[BlockData]
    model_name : str
        Isotherm model key (``none``/``henry``/``langmuir``/``freundlich``).
    n_starts : int
        Number of starting points (first is the model default, rest random
        within bounds) — guards against local minima.
    n_cells : int
        Spatial discretisation; must match the methane calibration.
    seed : int
        RNG seed for reproducible multistart.
    exclude : (block, injection_number), optional
        Hold this injection out (used by leave-one-dose-out).
    x0 : np.ndarray, optional
        Warm start.  When given it replaces the multistart set — used by
        leave-one-dose-out, where the full-data optimum is an excellent
        starting point and re-running a full multistart per fold is wasteful.
    verbose : bool

    Returns
    -------
    FitResult
    """
    model = get_model(model_name)
    lo, hi = model.bounds

    def resid(p):
        return _residuals(blocks, model, np.asarray(p), n_cells, exclude)

    # --- Zero-parameter model: evaluate directly ---
    if model.n_params == 0:
        r = resid(np.array([]))
        return _package_fit(model, np.array([]), r, blocks, n_cells,
                            converged=True, n_starts=1, n_successful=1,
                            jac=None, exclude=exclude,
                            message="No free parameters (transport-only).")

    if x0 is not None:
        starts = [np.asarray(x0, dtype=float)]
    else:
        starts = _multistart_points(model, n_starts, seed)

    return _run_starts(model, starts, resid, lo, hi, blocks, n_cells, exclude,
                       verbose)


def _multistart_points(model, n_starts: int, seed: int) -> list[np.ndarray]:
    """Default start plus random draws within bounds (log-uniform if wide)."""
    lo, hi = model.bounds
    rng = np.random.default_rng(seed)
    starts = [np.array(model.initial, dtype=float)]
    for _ in range(max(0, n_starts - 1)):
        # Log-uniform sampling for positive-scale parameters.
        s = []
        for l, h, init in zip(lo, hi, model.initial):
            if l > 0 and h / max(l, 1e-30) > 100:
                s.append(float(np.exp(rng.uniform(np.log(l), np.log(h)))))
            else:
                s.append(float(rng.uniform(l, h)))
        starts.append(np.array(s))
    return starts


def _run_starts(model, starts, resid, lo, hi, blocks, n_cells, exclude,
                verbose):
    """Optimise from each start and package the best result."""
    from scipy.optimize import least_squares

    best = None
    n_ok = 0
    for i, x0 in enumerate(starts):
        x0 = np.clip(x0, lo, hi)
        try:
            sol = least_squares(resid, x0, bounds=(lo, hi), method="trf",
                                x_scale="jac", max_nfev=200 * (model.n_params + 1))
        except Exception as exc:  # pragma: no cover - optimiser safety net
            if verbose:
                print(f"      start {i}: failed ({exc})")
            continue
        n_ok += 1
        cost = float(sol.cost)
        if best is None or cost < best.cost:
            best = sol

    if best is None:
        raise RuntimeError(f"{model.name}: all optimiser starts failed")

    return _package_fit(model, np.asarray(best.x), np.asarray(best.fun), blocks,
                        n_cells, converged=bool(best.success), n_starts=len(starts),
                        n_successful=n_ok, jac=np.asarray(best.jac),
                        exclude=exclude, message=str(best.message))


def _package_fit(model, params, residual, blocks, n_cells, *, converged,
                 n_starts, n_successful, jac, exclude, message) -> FitResult:
    """Assemble statistics, uncertainty and identifiability for a fit."""
    resid = np.asarray(residual, dtype=float)
    rss = float(np.sum(resid ** 2))
    n = int(resid.size)
    k = int(model.n_params)

    # R² against the observed (normalised) traces.
    obs = []
    for blk in blocks:
        for inj in blk.injections:
            if exclude is not None and (inj.block, inj.injection_number) == exclude:
                continue
            obs.append(inj.c_out_mol_m3 / inj.peak_scale)
    obs_v = np.concatenate(obs) if obs else np.array([0.0])
    ss_tot = float(np.sum((obs_v - obs_v.mean()) ** 2))
    r2 = 1.0 - rss / ss_tot if ss_tot > 0 else float("nan")
    rmse = float(np.sqrt(rss / n)) if n else float("nan")

    # Information criteria.  NOTE: chromatographic residuals are strongly
    # autocorrelated, so these are comparative guides, not absolute evidence.
    if n > k + 1 and rss > 0:
        aic = n * math.log(rss / n) + 2 * k
        aicc = aic + (2 * k * (k + 1)) / (n - k - 1)
        bic = n * math.log(rss / n) + k * math.log(n)
    else:
        aicc = bic = float("nan")

    # --- Parameter uncertainty and identifiability -------------------------
    # Two complementary tests, because they fail in different situations:
    #
    #   (a) relative standard error — scales with the residual, so it catches
    #       "the data are too noisy to pin this down";
    #   (b) Jacobian conditioning and parameter correlation — scale-free, so it
    #       catches STRUCTURAL confounding (e.g. Langmuir q_s and K_L appearing
    #       only as their product in the near-linear regime) even when the fit
    #       is essentially exact and (a) would wrongly report high precision.
    se = np.full(k, np.nan)
    cond = float("nan")
    max_corr = float("nan")
    if k and jac is not None and n > k:
        try:
            JTJ = jac.T @ jac
            cov = np.linalg.inv(JTJ) * (rss / (n - k))
            se = np.sqrt(np.clip(np.diag(cov), 0.0, np.inf))
            if k > 1:
                denom = np.outer(se, se)
                with np.errstate(divide="ignore", invalid="ignore"):
                    corr = np.where(denom > 0, cov / denom, 0.0)
                off = corr[~np.eye(k, dtype=bool)]
                max_corr = float(np.nanmax(np.abs(off))) if off.size else 0.0
                # Scale-free conditioning: use the correlation matrix, not J'J,
                # so differing parameter units cannot masquerade as collinearity.
                cond = float(np.linalg.cond(corr))
            else:
                max_corr = 0.0
                cond = 1.0
        except np.linalg.LinAlgError:
            se = np.full(k, np.inf)      # singular → unidentifiable
            cond = float("inf")
    rse = np.array([s / abs(p) if (p != 0 and np.isfinite(s)) else np.inf
                    for s, p in zip(se, params)]) if k else np.array([])

    structurally_confounded = (
        (np.isfinite(max_corr) and max_corr > CORRELATION_LIMIT)
        or (np.isfinite(cond) and cond > CONDITION_LIMIT)
        or not np.isfinite(cond)
    )
    identifiable = {
        nm: bool(r < IDENTIFIABILITY_RSE_LIMIT and not structurally_confounded)
        for nm, r in zip(model.param_names, rse)
    }

    # Mass-balance diagnostics and per-injection RMSE.
    mbs, per_inj = [], {}
    for blk in blocks:
        for inj in blk.injections:
            if exclude is not None and (inj.block, inj.injection_number) == exclude:
                continue
            try:
                res = predict_injection(inj, blk, model, params, n_cells)
                mbs.append(res.mass_balance)
                d = (inj.c_out_mol_m3 - res.c_out) / inj.peak_scale
                per_inj[f"{inj.block}:{inj.injection_number}"] = float(
                    np.sqrt(np.mean(d ** 2)))
            except Exception:
                mbs.append(float("nan"))

    return FitResult(
        model_name=model.name, params=params, param_names=list(model.param_names),
        param_units=list(model.param_units), std_errors=se, rel_std_errors=rse,
        identifiable=identifiable, condition_number=cond,
        max_abs_correlation=max_corr, converged=converged, n_starts=n_starts,
        n_successful_starts=n_successful, rss=rss, n_points=n, n_params=k,
        rmse_normalised=rmse, r_squared=r2, aicc=aicc, bic=bic,
        mass_balance_mean=(float(np.mean(mbs))
                           if mbs and np.all(np.isfinite(mbs)) else float("nan")),
        mass_balance_min=(float(np.min(mbs))
                          if mbs and np.all(np.isfinite(mbs)) else float("nan")),
        per_injection_rmse=per_inj,
        cooperative=is_cooperative(model, params) if k else False,
        message=message,
    )


# ---------------------------------------------------------------------------
# 4. Model comparison
# ---------------------------------------------------------------------------

def leave_one_dose_out(
    blocks: list[BlockData], model_name: str, n_cells: int = DEFAULT_N_CELLS,
    warm_start: np.ndarray | None = None, verbose: bool = False,
) -> float:
    """Mean held-out normalised RMSE over leave-one-dose-out refits.

    Each injection is removed in turn, the model refitted on the remainder, and
    the omitted peak predicted.  This penalises models that fit in-sample only
    through flexibility.

    Each fold is warm-started from ``warm_start`` (normally the full-data
    optimum), which is both far cheaper than repeating a full multistart and
    more stable, since consecutive folds differ by only one injection.
    """
    model = get_model(model_name)
    errs = []
    for blk in blocks:
        for inj in blk.injections:
            key = (inj.block, inj.injection_number)
            try:
                fit = fit_model(blocks, model_name, n_cells=n_cells,
                                exclude=key, x0=warm_start, verbose=False)
                res = predict_injection(inj, blk, model, fit.params, n_cells)
                d = (inj.c_out_mol_m3 - res.c_out) / inj.peak_scale
                errs.append(float(np.sqrt(np.mean(d ** 2))))
            except Exception:
                errs.append(float("nan"))
    return float(np.nanmean(errs)) if errs else float("nan")


def transport_sensitivity(
    blocks: list[BlockData],
    model_name: str,
    base_params: np.ndarray,
    n_cells: int = DEFAULT_N_CELLS,
    verbose: bool = False,
) -> pd.DataFrame:
    """Refit the isotherm under perturbed transport and reporting assumptions.

    The adsorption parameters are only meaningful if they survive plausible
    changes in the *fixed* inputs.  Each scenario perturbs one assumption and
    refits (warm-started from ``base_params``):

    - ``t0 ± 1 SD``      — methane dead-time reproducibility, only when every
      block has at least two usable markers and a nonzero finite sample SD;
    - ``N × 0.8 / 1.25`` — the effective plate number (axial dispersion);
    - ``t_inj × 0.8/1.25`` — the inlet-pulse width assumption;
    - ``window 0–2 min`` — a narrower integration/peak window.

    Returns
    -------
    pd.DataFrame
        One row per scenario with the refitted parameters and the percentage
        change from the base fit.
    """
    import copy

    model = get_model(model_name)
    rows = []

    def _record(label, blks):
        try:
            fit = fit_model(blks, model_name, n_cells=n_cells,
                            x0=base_params, verbose=False)
            row = {"scenario": label, "rmse_normalised": fit.rmse_normalised}
            for nm, v, b in zip(model.param_names, fit.params, base_params):
                row[nm] = v
                row[f"{nm}_pct_change"] = (
                    100.0 * (v - b) / abs(b) if b else float("nan"))
            row["all_params_identifiable"] = fit.all_identifiable
            rows.append(row)
        except Exception as exc:                      # pragma: no cover
            rows.append({"scenario": label, "error": str(exc)})

    base_row = {"scenario": "base", "rmse_normalised": float("nan")}
    for nm, b in zip(model.param_names, base_params):
        base_row[nm] = b
        base_row[f"{nm}_pct_change"] = 0.0
    rows.append(base_row)

    def _perturb(fn):
        blks = []
        for b in blocks:
            nb = copy.copy(b)
            nb.transport = copy.copy(b.transport)
            nb.geometry = copy.copy(b.geometry)
            nb.injections = [copy.copy(i) for i in b.injections]
            fn(nb)
            blks.append(nb)
        return blks

    def _shift_t0(nb, delta):
        nb.transport.t0_min += delta
        nb.transport.void_volume_m3 = (
            nb.geometry.flow_col_m3_min * nb.transport.t0_min)
        nb.geometry = make_geometry(
            nb.geometry.sample_mass_g,
            nb.geometry.flow_col_m3_min,
            nb.transport.t0_min,
        )
        for inj in nb.injections:
            if inj.assigned_t0_min is not None:
                inj.assigned_t0_min += delta

    t0_sensitivity_available = all(
        np.isfinite(b.transport.t0_sd_min) and b.transport.t0_sd_min > 0
        for b in blocks
    )
    if t0_sensitivity_available:
        for sign in (+1.0, -1.0):
            lbl = f"t0 {'+' if sign > 0 else '-'}1 SD"
            _record(lbl, _perturb(
                lambda nb, s=sign: _shift_t0(
                    nb, s * nb.transport.t0_sd_min)))

    if any(i.transport_mode == "bracket_interpolated"
           for b in blocks for i in b.injections):
        def _use_bracket_endpoint(nb, endpoint):
            for inj in nb.injections:
                value = (inj.pre_bracket_t0_min if endpoint == "pre"
                         else inj.post_bracket_t0_min)
                if value is not None:
                    inj.assigned_t0_min = value
                    inj.transport_mode = f"bracket_{endpoint}"
        _record("all probes use pre-bracket t0", _perturb(
            lambda nb: _use_bracket_endpoint(nb, "pre")))
        _record("all probes use post-bracket t0", _perturb(
            lambda nb: _use_bracket_endpoint(nb, "post")))
    for f in (0.8, 1.25):
        _record(f"N x {f}", _perturb(
            lambda nb, ff=f: setattr(nb.transport, "plate_number",
                                     nb.transport.plate_number * ff)))
        _record(f"t_inj x {f}", _perturb(
            lambda nb, ff=f: setattr(nb.transport, "t_inj_min",
                                     nb.transport.t_inj_min * ff)))

    # Narrower peak window: truncate every trace at 2 minutes.
    def _windowed():
        blks = []
        for b in blocks:
            nb = copy.copy(b)
            nb.injections = []
            for inj in b.injections:
                m = inj.time_min <= 2.0
                ni = copy.copy(inj)
                ni.time_min = inj.time_min[m]
                ni.c_out_mol_m3 = inj.c_out_mol_m3[m]
                nb.injections.append(ni)
            blks.append(nb)
        return blks

    _record("window 0-2 min", _windowed())
    return pd.DataFrame(rows)


def compare_models(
    blocks: list[BlockData],
    model_names: tuple[str, ...] = ("none", "henry", "langmuir", "freundlich"),
    n_cells: int = DEFAULT_N_CELLS,
    n_starts: int = 4,
    do_lodo: bool = True,
    lodo_models: set[str] | None = None,
    verbose: bool = True,
) -> tuple[dict[str, FitResult], pd.DataFrame]:
    """Fit every candidate model and build the comparison table."""
    fits: dict[str, FitResult] = {}
    rows = []
    for name in model_names:
        if verbose:
            print(f"    fitting {name} ...")
        fit = fit_model(blocks, name, n_starts=n_starts, n_cells=n_cells,
                        verbose=verbose)
        fits[name] = fit
        run_lodo = do_lodo and (lodo_models is None or name in lodo_models)
        lodo = (leave_one_dose_out(blocks, name, n_cells=n_cells,
                                   warm_start=fit.params if fit.n_params else None)
                if run_lodo else float("nan"))
        rows.append({
            "model": name,
            "description": get_model(name).description,
            "n_params": fit.n_params,
            "rmse_normalised": fit.rmse_normalised,
            "r_squared": fit.r_squared,
            "aicc": fit.aicc,
            "bic": fit.bic,
            "lodo_rmse": lodo,
            "mass_balance_mean": fit.mass_balance_mean,
            "converged": fit.converged,
            "all_params_identifiable": fit.all_identifiable,
            "cooperative": fit.cooperative,
            "params": ", ".join(
                f"{n}={v:.4g}" for n, v in zip(fit.param_names, fit.params)),
        })
    table = pd.DataFrame(rows)
    if not table.empty:
        best = table.sort_values("aicc").iloc[0]["model"]
        table["best_by_aicc"] = table["model"] == best
    return fits, table


# ---------------------------------------------------------------------------
# 5. SSA guardrail (mandatory)
# ---------------------------------------------------------------------------

@dataclass
class SSAVerdict:
    """Outcome of the SSA guardrail."""

    ssa_m2_g: float | None
    reportable: bool
    reason: str


def compute_ssa_if_identifiable(
    fit: FitResult,
    blocks: list[BlockData],
    cross_section_m2: float,
    saturation_fraction_required: float = 0.2,
) -> SSAVerdict:
    """Return an SSA **only** when it is structurally and numerically earned.

    All three conditions must hold:

    1. the fitted model has a structural finite capacity (a monolayer exists);
    2. that capacity parameter is identifiable (finite, positive, relative
       standard error below the limit);
    3. the measured concentration range actually approaches saturation — at
       least ``saturation_fraction_required`` of ``q_s`` is reached.

    Otherwise no number is produced.  The measured/predicted P/P0 ratio is
    never converted into an SSA correction.
    """
    from igc_analysis.constants import N_AVOGADRO

    model = get_model(fit.model_name)

    if not model.has_finite_capacity:
        return SSAVerdict(None, False,
                          f"Model '{fit.model_name}' has no finite monolayer "
                          f"capacity by construction — a geometric SSA is not "
                          f"defined and is not reported.")

    cap = model.capacity_param
    if cap not in fit.param_names:
        return SSAVerdict(None, False,
                          f"Capacity parameter {cap!r} absent from the fit.")
    i = fit.param_names.index(cap)
    q_s = float(fit.params[i])

    if not (np.isfinite(q_s) and q_s > 0):
        return SSAVerdict(None, False, f"Fitted {cap} is not positive/finite.")
    if not fit.identifiable.get(cap, False):
        rse = fit.rel_std_errors[i] if i < len(fit.rel_std_errors) else float("nan")
        return SSAVerdict(None, False,
                          f"{cap} is not identifiable over the measured range "
                          f"(relative SE {rse:.0%} > "
                          f"{IDENTIFIABILITY_RSE_LIMIT:.0%}).")

    c_max = max((float(np.max(inj.c_out_mol_m3))
                 for blk in blocks for inj in blk.injections), default=0.0)
    q_at_max = float(model.q(np.array([c_max]), fit.params)[0])
    frac = q_at_max / q_s if q_s > 0 else 0.0
    if frac < saturation_fraction_required:
        return SSAVerdict(None, False,
                          f"Measured range reaches only {frac:.1%} of the "
                          f"fitted monolayer capacity (need "
                          f"{saturation_fraction_required:.0%}); {cap} is "
                          f"extrapolated, so SSA is not reported.")

    ssa = q_s * N_AVOGADRO * cross_section_m2
    return SSAVerdict(float(ssa), True,
                      f"{cap} identifiable and {frac:.0%} of capacity reached.")


def recovered_isotherm(fit: FitResult, blocks: list[BlockData],
                       n_points: int = 200) -> pd.DataFrame:
    """Tabulate q(c) and dq/dc over the measured concentration range."""
    model = get_model(fit.model_name)
    c_max = max((float(np.max(inj.c_out_mol_m3))
                 for blk in blocks for inj in blk.injections), default=1.0)
    c = np.linspace(0.0, c_max, n_points)
    # A representative P/P0 axis using the mean measured conditions.
    T = float(np.mean([inj.temp_col_K for blk in blocks for inj in blk.injections]))
    Ps = float(np.mean([inj.p_sat_Pa for blk in blocks for inj in blk.injections]))
    return pd.DataFrame({
        "c_mol_m3": c,
        "pp0": c * R_GAS * T / Ps if Ps > 0 else np.zeros_like(c),
        "q_mol_g": model.q(c, fit.params),
        "dqdc_m3_g": model.dqdc(c, fit.params),
        "model": fit.model_name,
    })
