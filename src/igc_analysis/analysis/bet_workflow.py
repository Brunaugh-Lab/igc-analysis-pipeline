"""End-to-end BET workflow for validated neutral-data bundles.

The calculation uses eluted peak-apex concentration, per-injection measured
column flow, James--Martin pressure correction, matched dead-time conventions,
sensitivity checks, and strict reportability gating.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from igc_analysis.analysis.bet import (
    BETResult,
    MethaneStats,
    ProbeSelection,
    _InjectionRecord,
    _assemble_result,
    _is_clipped,
    _pct_delta,
    bet_quality_checks,
    classify_isotherm,
    compute_bet_diagnostics,
    james_martin_j,
)
from igc_analysis.analysis.full_peak import (
    _neutral_calibrated_amount,
    _neutral_condition_mean,
    _neutral_trace,
)
from igc_analysis.analysis.peak_detection import process_chromatogram
from igc_analysis.constants import T_STANDARD_K
from igc_analysis.io.neutral_data import NeutralBundle, read_neutral_bundle


@dataclass(frozen=True)
class _EvaluatedProbeProperties:
    """Probe identity for bundles that carry evaluated P_sat per injection."""

    name: str
    cross_section_m2: float

    @property
    def antoine(self) -> dict[str, float]:
        return {}

    def p_sat(self, temperature_K: float) -> float:
        raise RuntimeError(
            "neutral BET must use each component's declared evaluated "
            "saturation vapour pressure"
        )


def _condition_with_role(
    conditions: pd.DataFrame,
    injection_id: str,
    quantity: str,
) -> tuple[float | None, str | None]:
    for role in ("measured", "target"):
        value = _neutral_condition_mean(
            conditions, injection_id, quantity, value_role=role
        )
        if value is not None:
            return value, role
    return None, None


def _column_conditions(
    conditions: pd.DataFrame,
    injection_id: str,
) -> tuple[float, float, float, tuple[str, ...], str]:
    """Return temperature K, column flow mL/min, standard flow, provenance."""

    temperature_K, temperature_role = _condition_with_role(
        conditions, injection_id, "column_temperature"
    )
    if temperature_K is None or temperature_K <= 0:
        raise ValueError(f"{injection_id}: column temperature is missing")

    flow_column_m3_s, flow_role = _condition_with_role(
        conditions, injection_id, "flow_column"
    )
    flow_standard_m3_s: float | None = None
    if flow_column_m3_s is None:
        flow_standard_m3_s, flow_role = _condition_with_role(
            conditions, injection_id, "flow_standard"
        )
        if flow_standard_m3_s is not None:
            flow_column_m3_s = flow_standard_m3_s * temperature_K / T_STANDARD_K
        flow_quantity = "flow_standard"
    else:
        flow_standard_m3_s = flow_column_m3_s * T_STANDARD_K / temperature_K
        flow_quantity = "flow_column"

    if flow_column_m3_s is None or flow_column_m3_s <= 0:
        raise ValueError(f"{injection_id}: column or standard flow is missing")
    if flow_standard_m3_s is None:
        raise ValueError(f"{injection_id}: standard flow could not be resolved")

    roles = tuple(sorted({str(temperature_role), str(flow_role)}))
    flow_rows = conditions[
        (conditions["injection_id"] == injection_id)
        & (conditions["quantity"] == flow_quantity)
        & (conditions["value_role"] == flow_role)
    ]
    channels = sorted({
        str(value).strip() for value in flow_rows["source_channel"]
        if str(value).strip()
    })
    flow_channel = channels[0] if channels else "declared-target"
    return (
        temperature_K,
        flow_column_m3_s * 60.0 * 1e6,
        flow_standard_m3_s * 60.0 * 1e6,
        roles,
        flow_channel,
    )


@dataclass(frozen=True)
class _PressureResolution:
    factor: float
    roles: tuple[str, ...]
    basis: str


def _pressure_factor(
    conditions: pd.DataFrame,
    injection_id: str,
    *,
    enabled: bool,
    ambient_pressure_pa: float,
) -> _PressureResolution:
    if not enabled:
        return _PressureResolution(1.0, (), "disabled")

    inlet, inlet_role = _condition_with_role(
        conditions, injection_id, "pressure_inlet"
    )
    outlet, outlet_role = _condition_with_role(
        conditions, injection_id, "pressure_outlet"
    )
    pressure_drop, drop_role = _condition_with_role(
        conditions, injection_id, "pressure_drop"
    )
    if inlet is not None and inlet <= 0:
        raise ValueError(f"{injection_id}: pressure_inlet must be positive")
    if outlet is not None and outlet <= 0:
        raise ValueError(f"{injection_id}: pressure_outlet must be positive")
    if pressure_drop is not None and pressure_drop < 0:
        raise ValueError(f"{injection_id}: pressure_drop cannot be negative")
    used_ambient_outlet = outlet is None
    if used_ambient_outlet:
        outlet = ambient_pressure_pa
    if inlet is None and pressure_drop is not None:
        inlet = outlet + pressure_drop
    if inlet is None:
        raise ValueError(
            f"{injection_id}: pressure correction requires pressure_inlet "
            "or pressure_drop"
        )
    if inlet < outlet:
        raise ValueError(
            f"{injection_id}: pressure_inlet cannot be below pressure_outlet"
        )
    if pressure_drop is not None and not math.isclose(
        inlet - outlet, pressure_drop, rel_tol=0.02, abs_tol=100.0
    ):
        raise ValueError(
            f"{injection_id}: pressure_inlet, pressure_outlet, and "
            "pressure_drop are inconsistent"
        )
    roles = tuple(sorted({
        str(role) for role in (inlet_role, outlet_role, drop_role)
        if role is not None
    }))
    if inlet_role is not None and not used_ambient_outlet:
        basis = "declared_absolute_inlet_outlet"
    elif inlet_role is not None:
        basis = "declared_absolute_inlet_plus_ambient_outlet"
    elif not used_ambient_outlet:
        basis = "declared_drop_plus_declared_absolute_outlet"
    else:
        basis = "declared_drop_plus_ambient_absolute_outlet"
    return _PressureResolution(james_martin_j(inlet, outlet), roles, basis)


def _condition_source(roles: set[str]) -> str:
    if roles == {"measured"}:
        return "measured"
    if roles == {"target"}:
        return "method_target"
    return "mixed"


def run_bet_from_neutral(
    bundle_dir: str | Path | NeutralBundle,
    *,
    p0_min: float = 0.05,
    p0_max: float = 0.35,
    probe: str = "auto",
    retention_mode: str = "peak_max",
    concentration_mode: str = "eluted",
    origin: str = "legacy",
    pressure_correction: bool = True,
    ambient_pressure_pa: float = 101325.0,
    sensitivity: bool = True,
) -> BETResult:
    """Run corrected BET analysis from one validated neutral bundle."""

    if retention_mode not in {"peak_max", "cofm"}:
        raise ValueError("retention_mode must be 'peak_max' or 'cofm'")
    if concentration_mode not in {"eluted", "loop"}:
        raise ValueError("concentration_mode must be 'eluted' or 'loop'")
    if origin not in {"legacy", "rectangular", "linear"}:
        raise ValueError("origin must be 'legacy', 'rectangular', or 'linear'")
    if not (0 <= p0_min < p0_max < 1):
        raise ValueError("BET bounds must satisfy 0 <= p0_min < p0_max < 1")
    if ambient_pressure_pa <= 0:
        raise ValueError("ambient_pressure_pa must be positive")

    bundle = (
        bundle_dir
        if isinstance(bundle_dir, NeutralBundle)
        else read_neutral_bundle(bundle_dir)
    )
    experiment = bundle.table("experiment.csv")
    injections = bundle.table("injections.csv").sort_values("sequence_index")
    components = bundle.table("injection_components.csv")
    properties = bundle.table("probe_properties.csv")
    calibrations = bundle.table("calibration.csv")
    conditions = bundle.table("conditions.csv")
    traces = bundle.table("traces.csv")

    if len(experiment) != 1:
        raise ValueError("BET requires exactly one experiment per neutral bundle")
    if len(set(injections["block_id"].astype(str))) != 1:
        raise ValueError("BET requires exactly one acquisition block per bundle")

    experiment_row = experiment.iloc[0]
    mass_g = float(experiment_row["sample_mass_g"])
    loop_raw = str(experiment_row["injection_loop_volume_m3"])
    loop_volume_m3 = float(loop_raw) if loop_raw else float("nan")
    if concentration_mode == "loop" and (
        not math.isfinite(loop_volume_m3) or loop_volume_m3 <= 0
    ):
        raise ValueError("loop concentration requires injection_loop_volume_m3")

    merged = components.merge(properties, on="probe_id", validate="many_to_one")
    analytes = merged[merged["component_role"] == "analyte"]
    probe_rows = injections[injections["role"] == "probe"]
    dead_time_rows = injections[injections["role"] == "dead_time"]
    if probe_rows.empty or dead_time_rows.empty:
        raise ValueError("BET requires probe and dead-time injections")

    available = analytes[["probe_id", "probe_name"]].drop_duplicates()
    if probe == "auto":
        if len(available) != 1:
            raise ValueError(
                "BET requires exactly one analyte probe per bundle; found "
                f"{sorted(available['probe_name'].astype(str))}"
            )
        selected = available.iloc[0]
    else:
        matches = available[
            available["probe_name"].astype(str).str.casefold() == probe.casefold()
        ]
        if len(matches) != 1:
            raise ValueError(
                f"requested probe {probe!r} does not uniquely match the bundle"
            )
        selected = matches.iloc[0]
    selected_probe_id = str(selected["probe_id"])
    selected_probe_name = str(selected["probe_name"])
    property_row = properties[properties["probe_id"] == selected_probe_id]
    if len(property_row) != 1:
        raise ValueError("selected probe properties are missing or ambiguous")
    cross_section_m2 = float(property_row.iloc[0]["cross_section_m2"])
    if not math.isfinite(cross_section_m2) or cross_section_m2 <= 0:
        raise ValueError("BET requires a positive declared molecular cross-section")
    props = _EvaluatedProbeProperties(selected_probe_name, cross_section_m2)

    selected_components = components[
        (components["probe_id"] == selected_probe_id)
        & (components["component_role"] == "analyte")
    ]
    selected_ids = set(selected_components["injection_id"].astype(str))
    probe_rows = probe_rows[probe_rows["injection_id"].isin(selected_ids)]
    if probe_rows.empty:
        raise ValueError("selected probe has no probe-role injections")

    dead_max: list[float] = []
    dead_cofm: list[float] = []
    for row in dead_time_rows.itertuples(index=False):
        time_min, signal_raw, _ = _neutral_trace(
            traces, str(row.injection_id), str(row.detector_channel)
        )
        peak = process_chromatogram(time_min, signal_raw)
        dead_max.append(float(peak["peak_max_time"]))
        dead_cofm.append(float(peak["peak_cofm"]))
    t0_max = float(np.mean(dead_max))
    t0_cofm = float(np.mean(dead_cofm))
    methane = MethaneStats(
        n=len(dead_max),
        mean_max_min=t0_max,
        mean_cofm_min=t0_cofm,
        sd_max_min=float(np.std(dead_max)) if len(dead_max) > 1 else 0.0,
        range_max_min=float(np.max(dead_max) - np.min(dead_max)),
        drift_max_min=dead_max[-1] - dead_max[0] if len(dead_max) >= 2 else 0.0,
    )

    records: list[_InjectionRecord] = []
    standard_flows_mL_min: list[float] = []
    flow_source_channels: set[str] = set()
    pressure_bases: set[str] = set()
    pressure_roles: set[str] = set()
    for number, row in enumerate(probe_rows.itertuples(index=False), start=1):
        injection_id = str(row.injection_id)
        component_rows = selected_components[
            selected_components["injection_id"] == injection_id
        ]
        if len(component_rows) != 1:
            raise ValueError(f"{injection_id}: expected one selected analyte component")
        component = component_rows.iloc[0]
        p_sat_raw = str(component["saturation_vapor_pressure_Pa"])
        if not p_sat_raw:
            raise ValueError(f"{injection_id}: saturation vapour pressure is required")
        p_sat_Pa = float(p_sat_raw)
        if not math.isfinite(p_sat_Pa) or p_sat_Pa <= 0:
            raise ValueError(f"{injection_id}: saturation vapour pressure is required")

        calibration_id = str(component["calibration_id"])
        calibration_rows = calibrations[
            calibrations["calibration_id"] == calibration_id
        ]
        if len(calibration_rows) != 1:
            raise ValueError(f"{injection_id}: calibration is missing or ambiguous")
        calibration = calibration_rows.iloc[0]
        if str(calibration["probe_id"]) != selected_probe_id:
            raise ValueError(f"{injection_id}: calibration probe does not match")

        time_min, signal_raw, signal_unit = _neutral_trace(
            traces, injection_id, str(row.detector_channel)
        )
        peak = process_chromatogram(
            time_min, signal_raw, dead_time_min=t0_max
        )
        baseline = peak["baseline_intercept"] + peak["baseline_gradient"] * time_min
        corrected = signal_raw - baseline
        amount_mol = _neutral_calibrated_amount(
            time_min, corrected, signal_unit, calibration
        )
        (temperature_K, flow_column_mL_min, flow_standard_mL_min,
         condition_roles, flow_source_channel) = (
            _column_conditions(conditions, injection_id)
        )
        standard_flows_mL_min.append(flow_standard_mL_min)
        flow_source_channels.add(flow_source_channel)
        pressure = _pressure_factor(
            conditions,
            injection_id,
            enabled=pressure_correction,
            ambient_pressure_pa=ambient_pressure_pa,
        )
        pressure_bases.add(pressure.basis)
        pressure_roles.update(pressure.roles)
        source = _condition_source(set(condition_roles) | set(pressure.roles))

        clipping_value = str(row.clipping_observed).casefold()
        clipped = clipping_value == "true" or _is_clipped(signal_raw)
        target_raw = str(row.target_coverage_fraction)
        records.append(_InjectionRecord(
            number=number,
            target_coverage=float(target_raw) if target_raw else None,
            peak=peak,
            n_injected_mol=amount_mol,
            area=float(peak["peak_area"]),
            temp_col_K=temperature_K,
            flow_col_mL_min=flow_column_mL_min,
            j_factor=pressure.factor,
            conditions_source=source,
            clipped=clipped,
            p_sat_Pa=p_sat_Pa,
            injection_id=injection_id,
        ))

    main_result = _assemble_result(
        records,
        mass_g=mass_g,
        props=props,
        V_loop_m3=loop_volume_m3,
        t0_max=t0_max,
        t0_cofm=t0_cofm,
        retention_mode=retention_mode,
        concentration_mode=concentration_mode,
        origin=origin,
        p0_min=p0_min,
        p0_max=p0_max,
    )
    main_result.sample_name = str(experiment_row["sample_id"])
    main_result.probe_selection = ProbeSelection(
        probe=selected_probe_name, source="neutral_bundle"
    )
    main_result.flow_sccm = float(np.mean(standard_flows_mL_min))
    main_result.flow_source_channels = tuple(sorted(flow_source_channels))
    main_result.pressure_source = (
        "disabled" if not pressure_correction else _condition_source(pressure_roles)
    )
    main_result.pressure_basis = tuple(sorted(pressure_bases))
    diagnostics = compute_bet_diagnostics(
        main_result, methane, p0_min, p0_max
    )
    main_result.diagnostics = diagnostics

    if sensitivity:
        alternative_retention = (
            "cofm" if retention_mode == "peak_max" else "peak_max"
        )
        retention_result = _assemble_result(
            records,
            mass_g=mass_g,
            props=props,
            V_loop_m3=loop_volume_m3,
            t0_max=t0_max,
            t0_cofm=t0_cofm,
            retention_mode=alternative_retention,
            concentration_mode=concentration_mode,
            origin=origin,
            p0_min=p0_min,
            p0_max=p0_max,
        )
        diagnostics.ssa_alt_retention = retention_result.ssa_m2_g
        diagnostics.alt_retention_mode = alternative_retention

        alternative_origin = "rectangular" if origin == "legacy" else "legacy"
        origin_result = _assemble_result(
            records,
            mass_g=mass_g,
            props=props,
            V_loop_m3=loop_volume_m3,
            t0_max=t0_max,
            t0_cofm=t0_cofm,
            retention_mode=retention_mode,
            concentration_mode=concentration_mode,
            origin=alternative_origin,
            p0_min=p0_min,
            p0_max=p0_max,
        )
        diagnostics.ssa_alt_origin = origin_result.ssa_m2_g
        diagnostics.alt_origin_strategy = alternative_origin

    main_result.qc = bet_quality_checks(main_result)
    if sensitivity and not math.isnan(main_result.ssa_m2_g):
        retention_delta = _pct_delta(
            diagnostics.ssa_alt_retention, main_result.ssa_m2_g
        )
        if not math.isnan(retention_delta) and retention_delta > 10.0:
            main_result.qc.retention_convention_sensitive = True
            main_result.qc.messages.append(
                f"SSA shifts {retention_delta:.0f}% between peak-max and CoM "
                f"({main_result.ssa_m2_g:.3f} vs "
                f"{diagnostics.ssa_alt_retention:.3f} m2/g)"
            )
        origin_delta = _pct_delta(
            diagnostics.ssa_alt_origin, main_result.ssa_m2_g
        )
        if not math.isnan(origin_delta) and origin_delta > 2.0:
            main_result.qc.origin_sensitive = True
            main_result.qc.messages.append(
                f"SSA shifts {origin_delta:.1f}% between '{origin}' and "
                f"'{diagnostics.alt_origin_strategy}' origin treatments "
                f"({main_result.ssa_m2_g:.3f} vs "
                f"{diagnostics.ssa_alt_origin:.3f} m2/g)"
            )

    main_result.classification = classify_isotherm(main_result)
    return main_result
