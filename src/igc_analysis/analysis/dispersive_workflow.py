"""Source-neutral Dorris--Gray dispersive surface-energy workflow.

The workflow preserves the final validated pre-split conventions while
removing source-specific ingestion: calibrated actual coverage, declared SSA
and probe-property provenance, measured per-injection conditions,
James--Martin correction, center-of-mass retention as primary, peak maximum as
a sensitivity result, and explicit interpolation/extrapolation diagnostics.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from igc_analysis.analysis.bet import _is_clipped
from igc_analysis.analysis.bet_workflow import _column_conditions, _pressure_factor
from igc_analysis.analysis.calibration import monolayer_capacity
from igc_analysis.analysis.dispersive import dorris_gray_gamma_d
from igc_analysis.analysis.full_peak import (
    _neutral_calibrated_amount,
    _neutral_condition_mean,
    _neutral_trace,
)
from igc_analysis.analysis.peak_detection import process_chromatogram
from igc_analysis.analysis.quality import run_qc_checks
from igc_analysis.io.neutral_data import NeutralBundle, read_neutral_bundle


GAMMA_D_WARNING_BOUNDS_MJ_M2 = (15.0, 80.0)
GAMMA_D_CRITICAL_BOUNDS_MJ_M2 = (5.0, 150.0)


@dataclass(frozen=True)
class DispersiveResult:
    """Complete calculation, diagnostics, and provenance for one bundle."""

    dataset_id: str
    sample_id: str
    specific_surface_area_m2_g: float
    surface_area_source: str
    temperature_K: float
    dead_time_cofm_min: float
    dead_time_peak_max_min: float
    injections: pd.DataFrame
    interpolated: pd.DataFrame
    additional_injections: pd.DataFrame
    additional_interpolated: pd.DataFrame
    gamma_d: pd.DataFrame
    qc: dict
    reportable: bool
    flow_source_channels: tuple[str, ...]
    pressure_basis: tuple[str, ...]
    pressure_roles: tuple[str, ...]
    properties_sources: tuple[str, ...]
    calibration_sources: tuple[str, ...]
    detector_gains: tuple[float, ...]


def _required_float(row: pd.Series, field: str, context: str) -> float:
    raw = str(row[field]).strip()
    if not raw:
        raise ValueError(f"{context}: {field} is required")
    value = float(raw)
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f"{context}: {field} must be finite and positive")
    return value


def _linear_value(x: np.ndarray, y: np.ndarray, target: float) -> tuple[float, str]:
    """Evaluate a piecewise-linear curve and label interpolation status."""

    if target < x[0]:
        i0, i1, status = 0, 1, "extrapolated_low"
    elif target > x[-1]:
        i0, i1, status = len(x) - 2, len(x) - 1, "extrapolated_high"
    else:
        return float(np.interp(target, x, y)), "interpolated"
    slope = (y[i1] - y[i0]) / (x[i1] - x[i0])
    return float(y[i0] + slope * (target - x[i0])), status


def _interpolate_retention(
    injections: pd.DataFrame,
    targets: list[float],
    *,
    vn_column: str,
    retention_mode: str,
    extrapolate: bool,
    allow_invalid_retention: bool = False,
) -> pd.DataFrame:
    rows: list[dict] = []
    for probe_id, group in injections.groupby("probe_id", sort=False):
        g = group.sort_values("actual_coverage").copy()
        invalid_retention_count = 0
        if allow_invalid_retention:
            valid = (
                np.isfinite(pd.to_numeric(g["actual_coverage"], errors="coerce"))
                & np.isfinite(pd.to_numeric(g[vn_column], errors="coerce"))
                & (pd.to_numeric(g["actual_coverage"], errors="coerce") > 0)
                & (pd.to_numeric(g[vn_column], errors="coerce") > 0)
            )
            invalid_retention_count = int((~valid).sum())
            g = g[valid].copy()
        if len(g) < 2:
            first = group.iloc[0]
            carbon_number = pd.to_numeric(
                pd.Series([first["carbon_number"]]), errors="coerce"
            ).iloc[0]
            for target in targets:
                rows.append({
                    "retention_mode": retention_mode,
                    "probe_id": str(probe_id),
                    "probe_name": str(first["probe_name"]),
                    "carbon_number": (
                        int(carbon_number) if pd.notna(carbon_number) else pd.NA
                    ),
                    "target_coverage": target,
                    "VN_mL_g": float("nan"),
                    "actual_coverage_min": (
                        float(g["actual_coverage"].min()) if len(g) else float("nan")
                    ),
                    "actual_coverage_max": (
                        float(g["actual_coverage"].max()) if len(g) else float("nan")
                    ),
                    "interpolation_status": "insufficient_positive_retention",
                    "invalid_retention_count": invalid_retention_count,
                })
            continue
        x = g["actual_coverage"].to_numpy(dtype=float)
        y = g[vn_column].to_numpy(dtype=float)
        if len(x) < 2:
            raise ValueError(
                f"probe {probe_id!r} requires at least two coverage points"
            )
        if len(np.unique(x)) != len(x):
            raise ValueError(
                f"probe {probe_id!r} has duplicate calibrated coverage values"
            )
        if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
            raise ValueError(f"probe {probe_id!r} has non-finite coverage or V_N")
        if np.any(x <= 0) or np.any(y <= 0):
            raise ValueError(f"probe {probe_id!r} has nonpositive coverage or V_N")

        first = g.iloc[0]
        for target in targets:
            below = target < x[0] and not math.isclose(
                target, x[0], rel_tol=1e-9, abs_tol=1e-12
            )
            above = target > x[-1] and not math.isclose(
                target, x[-1], rel_tol=1e-9, abs_tol=1e-12
            )
            outside = below or above
            evaluation_target = (
                float(x[0]) if math.isclose(target, x[0], rel_tol=1e-9, abs_tol=1e-12)
                else float(x[-1])
                if math.isclose(target, x[-1], rel_tol=1e-9, abs_tol=1e-12)
                else target
            )
            if outside and not extrapolate:
                value, status = float("nan"), "outside_measured_range"
            else:
                value, status = _linear_value(x, y, evaluation_target)
                if not math.isfinite(value) or value <= 0:
                    value, status = float("nan"), "extrapolated_nonphysical"
            carbon_number = pd.to_numeric(
                pd.Series([first["carbon_number"]]), errors="coerce"
            ).iloc[0]
            rows.append({
                "retention_mode": retention_mode,
                "probe_id": str(probe_id),
                "probe_name": str(first["probe_name"]),
                "carbon_number": (
                    int(carbon_number) if pd.notna(carbon_number) else pd.NA
                ),
                "target_coverage": target,
                "VN_mL_g": value,
                "actual_coverage_min": float(x[0]),
                "actual_coverage_max": float(x[-1]),
                "interpolation_status": status,
                "invalid_retention_count": invalid_retention_count,
            })
    return pd.DataFrame(rows)


def _fit_profile(interpolated: pd.DataFrame, temperature_K: float) -> pd.DataFrame:
    coverages = pd.DataFrame({
        "coverage": sorted(interpolated["target_coverage"].unique())
    })
    fit_input = interpolated.rename(columns={
        "target_coverage": "coverage",
        "VN_mL_g": "VN",
    })[["coverage", "carbon_number", "VN"]].dropna(subset=["VN"])
    fit_input["carbon_number"] = pd.to_numeric(
        fit_input["carbon_number"], errors="raise"
    )
    if fit_input.empty:
        fitted = pd.DataFrame(columns=[
            "coverage", "gamma_d_mJm2", "r_squared", "slope_Jmol",
            "intercept", "n_alkanes",
        ])
    else:
        fitted = dorris_gray_gamma_d(fit_input, temperature_K)
    profile = coverages.merge(fitted, on="coverage", how="left")
    profile["n_alkanes"] = profile["n_alkanes"].fillna(0).astype(int)
    return profile


def _append_flag(qc: dict, flag: dict) -> None:
    qc["flags"].append(flag)


def _finish_qc(qc: dict) -> None:
    n_warn = sum(flag["severity"] == "warning" for flag in qc["flags"])
    n_crit = sum(flag["severity"] == "critical" for flag in qc["flags"])
    qc["pass"] = n_crit == 0
    status = "PASS" if qc["pass"] else "FAIL"
    qc["summary"] = (
        f"QC {status}: {n_warn} warning(s), {n_crit} critical; "
        f"profile: {qc['profile_shape']}"
    )


def run_dispersive_from_neutral(
    bundle_dir: str | Path | NeutralBundle,
    *,
    homologous_probe_ids: tuple[str, ...] | list[str] | None = None,
    additional_probe_ids: tuple[str, ...] | list[str] | None = None,
    pressure_correction: bool = True,
    ambient_pressure_pa: float = 101325.0,
    extrapolate: bool = True,
    max_temperature_span_K: float = 1.0,
    max_flow_relative_span: float = 0.05,
) -> DispersiveResult:
    """Run coverage-resolved Dorris--Gray analysis on one neutral bundle."""

    if ambient_pressure_pa <= 0:
        raise ValueError("ambient_pressure_pa must be positive")
    if max_temperature_span_K < 0:
        raise ValueError("max_temperature_span_K cannot be negative")
    if max_flow_relative_span < 0:
        raise ValueError("max_flow_relative_span cannot be negative")

    bundle = (
        bundle_dir
        if isinstance(bundle_dir, NeutralBundle)
        else read_neutral_bundle(bundle_dir)
    )
    experiment = bundle.table("experiment.csv")
    injections = bundle.table("injections.csv").sort_values("sequence_index")
    components = bundle.table("injection_components.csv")
    traces = bundle.table("traces.csv")
    properties = bundle.table("probe_properties.csv")
    calibrations = bundle.table("calibration.csv")
    conditions = bundle.table("conditions.csv")

    if len(experiment) != 1:
        raise ValueError("dispersive analysis requires exactly one experiment row")
    if injections["block_id"].nunique() != 1:
        raise ValueError("dispersive analysis requires exactly one acquisition block")

    experiment_row = experiment.iloc[0]
    mass_g = _required_float(experiment_row, "sample_mass_g", "experiment")
    ssa_m2_g = _required_float(
        experiment_row, "specific_surface_area_m2_g", "experiment"
    )
    ssa_source = str(experiment_row["surface_area_source"]).strip()
    if not ssa_source:
        raise ValueError("experiment: surface_area_source is required")

    merged = components.merge(properties, on="probe_id", validate="many_to_one")
    analytes = merged[merged["component_role"] == "analyte"].copy()
    analytes["carbon_number"] = pd.to_numeric(
        analytes["carbon_number"], errors="coerce"
    )
    carbon_declared = analytes[analytes["carbon_number"].notna()].copy()
    if homologous_probe_ids is not None:
        selected_ids = tuple(str(value).strip() for value in homologous_probe_ids)
        if any(not value for value in selected_ids):
            raise ValueError("homologous probe IDs cannot be empty")
        if len(set(selected_ids)) != len(selected_ids):
            raise ValueError("homologous probe IDs must be unique")
        available_ids = set(analytes["probe_id"].astype(str))
        missing_ids = sorted(set(selected_ids) - available_ids)
        if missing_ids:
            raise ValueError(
                "selected homologous probe IDs are not analytes in this bundle: "
                + ", ".join(missing_ids)
            )
        homologous = carbon_declared[
            carbon_declared["probe_id"].astype(str).isin(selected_ids)
        ].copy()
        selected_without_carbon = sorted(
            set(selected_ids) - set(homologous["probe_id"].astype(str))
        )
        if selected_without_carbon:
            raise ValueError(
                "selected homologous probes require declared carbon_number: "
                + ", ".join(selected_without_carbon)
            )
    else:
        if carbon_declared["probe_id"].nunique() > 3:
            raise ValueError(
                "bundles with more than three carbon-numbered analytes require "
                "explicit homologous_probe_ids/--homologous-probe-id selection"
            )
        homologous = carbon_declared
    if homologous.empty:
        raise ValueError("dispersive analysis requires a homologous analyte series")

    probe_defs = homologous[[
        "probe_id", "probe_name", "carbon_number", "cross_section_m2",
        "properties_source",
    ]].drop_duplicates()
    if probe_defs["probe_id"].nunique() < 3:
        raise ValueError("dispersive analysis requires at least three homologous probes")
    if probe_defs["carbon_number"].nunique() != len(probe_defs):
        raise ValueError(
            "each homologous probe must have a unique carbon number; use "
            "homologous_probe_ids/--homologous-probe-id when the bundle also "
            "contains other carbon-numbered analytes"
        )
    for row in probe_defs.itertuples(index=False):
        if not str(row.properties_source).strip():
            raise ValueError(f"probe {row.probe_id!r} is missing properties_source")
        _required_float(
            pd.Series(row._asdict()),
            "cross_section_m2",
            f"probe {row.probe_id!r}",
        )
        carbon_number = float(row.carbon_number)
        if not carbon_number.is_integer():
            raise ValueError(
                f"probe {row.probe_id!r} requires an integer carbon_number"
        )

    additional = analytes.iloc[0:0].copy()
    if additional_probe_ids is not None:
        additional_ids = tuple(str(value).strip() for value in additional_probe_ids)
        if any(not value for value in additional_ids):
            raise ValueError("additional probe IDs cannot be empty")
        if len(set(additional_ids)) != len(additional_ids):
            raise ValueError("additional probe IDs must be unique")
        homologous_id_values = set(probe_defs["probe_id"].astype(str))
        overlap = sorted(set(additional_ids) & homologous_id_values)
        if overlap:
            raise ValueError(
                "additional and homologous probe selections must be disjoint: "
                + ", ".join(overlap)
            )
        available_ids = set(analytes["probe_id"].astype(str))
        missing_ids = sorted(set(additional_ids) - available_ids)
        if missing_ids:
            raise ValueError(
                "selected additional probe IDs are not analytes in this bundle: "
                + ", ".join(missing_ids)
            )
        additional = analytes[
            analytes["probe_id"].astype(str).isin(additional_ids)
        ].copy()
        additional_defs = additional[[
            "probe_id", "probe_name", "carbon_number", "cross_section_m2",
            "properties_source",
        ]].drop_duplicates()
        if additional_defs["probe_id"].nunique() != len(additional_ids):
            raise ValueError("each additional probe must have one property record")
        for row in additional_defs.itertuples(index=False):
            if not str(row.properties_source).strip():
                raise ValueError(f"probe {row.probe_id!r} is missing properties_source")
            _required_float(
                pd.Series(row._asdict()), "cross_section_m2",
                f"probe {row.probe_id!r}",
            )

    probe_rows = injections[injections["role"] == "probe"].copy()
    dead_rows = injections[injections["role"] == "dead_time"].copy()
    selected_components = pd.concat([homologous, additional], ignore_index=True)
    selected_injection_ids = set(selected_components["injection_id"].astype(str))
    probe_rows = probe_rows[
        probe_rows["injection_id"].isin(selected_injection_ids)
    ]
    if probe_rows.empty or dead_rows.empty:
        raise ValueError("dispersive analysis requires probe and dead-time injections")
    required_rows = pd.concat([probe_rows, dead_rows], ignore_index=True)
    detector_gains = np.array([
        _required_float(row, "detector_gain", str(row["injection_id"]))
        for _, row in required_rows.iterrows()
    ])
    selected_probe_defs = selected_components[["probe_id"]].drop_duplicates()
    selected_probe_components = selected_components[
        selected_components["injection_id"].isin(probe_rows["injection_id"])
    ]
    coverage_counts = selected_probe_components.groupby("probe_id")[
        "injection_id"
    ].nunique().reindex(selected_probe_defs["probe_id"].astype(str), fill_value=0)
    sparse_probes = coverage_counts[coverage_counts < 3]
    if not sparse_probes.empty:
        raise ValueError(
            "analysis requires at least three calibrated coverage points for "
            "every selected probe; insufficient: "
            + ", ".join(str(value) for value in sparse_probes.index)
        )

    dead_max: list[float] = []
    dead_cofm: list[float] = []
    dead_clipped: list[str] = []
    all_temperatures: list[float] = []
    all_flows: list[float] = []
    flow_channels: set[str] = set()
    pressure_bases: set[str] = set()
    pressure_roles: set[str] = set()
    property_sources: set[str] = set()
    calibration_sources: set[str] = set()
    for row in dead_rows.itertuples(index=False):
        injection_id = str(row.injection_id)
        (dead_temperature_K, dead_flow_mL_min, _, condition_roles,
         dead_flow_channel) = _column_conditions(conditions, injection_id)
        if set(condition_roles) != {"measured"}:
            raise ValueError(
                f"{injection_id}: dispersive analysis requires measured "
                "column temperature and flow"
            )
        all_temperatures.append(dead_temperature_K)
        all_flows.append(dead_flow_mL_min)
        flow_channels.add(dead_flow_channel)
        time_min, signal_raw, _ = _neutral_trace(
            traces, injection_id, str(row.detector_channel)
        )
        peak = process_chromatogram(time_min, signal_raw)
        dead_max.append(float(peak["peak_max_time"]))
        dead_cofm.append(float(peak["peak_cofm"]))
        declared = str(row.clipping_observed).casefold() == "true"
        if declared or _is_clipped(signal_raw):
            dead_clipped.append(injection_id)
    t0_max = float(np.mean(dead_max))
    t0_cofm = float(np.mean(dead_cofm))

    records: list[dict] = []
    probe_temperatures: list[float] = []
    clipped_probe_ids: list[str] = []

    selected_by_injection = selected_components.set_index(
        "injection_id", drop=False
    )
    for row in probe_rows.itertuples(index=False):
        injection_id = str(row.injection_id)
        analyte_components = components[
            (components["injection_id"] == injection_id)
            & (components["component_role"] == "analyte")
        ]
        if len(analyte_components) != 1:
            raise ValueError(
                f"{injection_id}: dispersive analysis requires exactly one "
                "analyte component"
            )
        component = selected_by_injection.loc[injection_id]
        if isinstance(component, pd.DataFrame):
            raise ValueError(f"{injection_id}: expected one selected analyte")
        target_raw = str(row.target_coverage_fraction).strip()
        if not target_raw:
            raise ValueError(f"{injection_id}: target_coverage_fraction is required")
        target_coverage = float(target_raw)
        if not math.isfinite(target_coverage) or target_coverage <= 0:
            raise ValueError(f"{injection_id}: target coverage must be positive")

        calibration_id = str(component["calibration_id"]).strip()
        calibration_rows = calibrations[
            calibrations["calibration_id"] == calibration_id
        ]
        if len(calibration_rows) != 1:
            raise ValueError(f"{injection_id}: calibration is missing or ambiguous")
        calibration = calibration_rows.iloc[0]
        if str(calibration["probe_id"]) != str(component["probe_id"]):
            raise ValueError(f"{injection_id}: calibration probe does not match")
        calibration_source = str(calibration["calibration_source"]).strip()
        if not calibration_source:
            raise ValueError(f"{injection_id}: calibration_source is required")

        time_min, signal_raw, signal_unit = _neutral_trace(
            traces, injection_id, str(row.detector_channel)
        )
        peak = process_chromatogram(time_min, signal_raw, dead_time_min=t0_max)
        baseline = peak["baseline_intercept"] + peak["baseline_gradient"] * time_min
        amount_mol = _neutral_calibrated_amount(
            time_min, signal_raw - baseline, signal_unit, calibration
        )
        if amount_mol <= 0:
            raise ValueError(f"{injection_id}: calibrated amount must be positive")

        cross_section_m2 = _required_float(
            component, "cross_section_m2", f"probe {component['probe_id']!r}"
        )
        capacity_mol = monolayer_capacity(ssa_m2_g, mass_g, cross_section_m2)
        actual_coverage = amount_mol / capacity_mol
        if not math.isfinite(actual_coverage) or actual_coverage <= 0:
            raise ValueError(f"{injection_id}: calibrated coverage is nonphysical")

        (temperature_K, flow_column_mL_min, _, condition_roles,
         flow_channel) = _column_conditions(conditions, injection_id)
        if set(condition_roles) != {"measured"}:
            raise ValueError(
                f"{injection_id}: dispersive analysis requires measured "
                "column temperature and flow"
            )
        pressure_rows = conditions[
            ~(
                conditions["quantity"].isin(
                    ["pressure_inlet", "pressure_outlet", "pressure_drop"]
                )
                & (conditions["value_role"] != "measured")
            )
        ]
        if pressure_correction:
            measured_inlet = _neutral_condition_mean(
                conditions,
                injection_id,
                "pressure_inlet",
                value_role="measured",
            )
            measured_drop = _neutral_condition_mean(
                conditions,
                injection_id,
                "pressure_drop",
                value_role="measured",
            )
            if measured_inlet is None and measured_drop is None:
                raise ValueError(
                    f"{injection_id}: pressure correction requires measured "
                    "pressure_inlet or pressure_drop; outlet may be measured "
                    "or use the configured ambient pressure"
                )
        pressure = _pressure_factor(
            pressure_rows,
            injection_id,
            enabled=pressure_correction,
            ambient_pressure_pa=ambient_pressure_pa,
        )
        net_max = float(peak["peak_max_time"]) - t0_max
        net_cofm = float(peak["peak_cofm"]) - t0_cofm
        vn_max = (
            net_max * flow_column_mL_min * pressure.factor / mass_g
            if net_max > 0 else float("nan")
        )
        vn_cofm = (
            net_cofm * flow_column_mL_min * pressure.factor / mass_g
            if net_cofm > 0 else float("nan")
        )

        declared_clip = str(row.clipping_observed).casefold() == "true"
        clipped = declared_clip or _is_clipped(signal_raw)
        if clipped:
            clipped_probe_ids.append(injection_id)
        probe_temperatures.append(temperature_K)
        all_temperatures.append(temperature_K)
        all_flows.append(flow_column_mL_min)
        flow_channels.add(flow_channel)
        pressure_bases.add(pressure.basis)
        pressure_roles.update(pressure.roles)
        property_sources.add(str(component["properties_source"]).strip())
        calibration_sources.add(calibration_source)
        records.append({
            "injection_id": injection_id,
            "sequence_index": int(row.sequence_index),
            "probe_id": str(component["probe_id"]),
            "probe_name": str(component["probe_name"]),
            "carbon_number": pd.to_numeric(
                pd.Series([component["carbon_number"]]), errors="coerce"
            ).iloc[0],
            "target_coverage": target_coverage,
            "actual_coverage": actual_coverage,
            "amount_mol": amount_mol,
            "peak_area": float(peak["peak_area"]),
            "peak_max_time_min": float(peak["peak_max_time"]),
            "peak_cofm_time_min": float(peak["peak_cofm"]),
            "asymmetry_factor": float(peak["asymmetry_factor"]),
            "tailing_factor": float(peak["tailing_factor"]),
            "cofm_peak_max_divergence_min": float(
                peak["com_max_divergence_min"]
            ),
            "cofm_peak_max_divergence_fraction": float(
                peak["com_max_divergence_frac"]
            ),
            "net_retention_peak_max_min": net_max,
            "net_retention_cofm_min": net_cofm,
            "VN_peak_max_mL_g": vn_max,
            "VN_cofm_mL_g": vn_cofm,
            "temperature_K": temperature_K,
            "flow_column_mL_min": flow_column_mL_min,
            "james_martin_factor": pressure.factor,
            "pressure_basis": pressure.basis,
            "flow_source_channel": flow_channel,
            "clipped": clipped,
            "calibration_id": calibration_id,
            "calibration_source": calibration_source,
            "cross_section_m2": cross_section_m2,
            "properties_source": str(component["properties_source"]).strip(),
        })

    injection_table = pd.DataFrame(records).sort_values("sequence_index")
    temperature_span = float(np.ptp(all_temperatures))
    mean_flow = float(np.mean(all_flows))
    flow_relative_span = float(np.ptp(all_flows) / mean_flow)
    if temperature_span > max_temperature_span_K:
        raise ValueError(
            "column temperature is not stable across selected injections: "
            f"span {temperature_span:.3g} K exceeds {max_temperature_span_K:.3g} K"
        )
    if flow_relative_span > max_flow_relative_span:
        raise ValueError(
            "column flow is not stable across selected injections: relative "
            f"span {flow_relative_span:.3g} exceeds {max_flow_relative_span:.3g}"
        )
    temperature_K = float(np.mean(probe_temperatures))

    homologous_id_values = set(probe_defs["probe_id"].astype(str))
    additional_id_values = set(additional["probe_id"].astype(str))
    homologous_injections = injection_table[
        injection_table["probe_id"].isin(homologous_id_values)
    ].copy()
    additional_injections = injection_table[
        injection_table["probe_id"].isin(additional_id_values)
    ].copy()
    targets = sorted(set(homologous_injections["target_coverage"].astype(float)))
    cofm = _interpolate_retention(
        homologous_injections,
        targets,
        vn_column="VN_cofm_mL_g",
        retention_mode="cofm",
        extrapolate=extrapolate,
    )
    peak_max = _interpolate_retention(
        homologous_injections,
        targets,
        vn_column="VN_peak_max_mL_g",
        retention_mode="peak_max",
        extrapolate=extrapolate,
    )
    interpolated = pd.concat([cofm, peak_max], ignore_index=True)
    if additional_injections.empty:
        additional_interpolated = pd.DataFrame(columns=interpolated.columns)
    else:
        additional_interpolated = pd.concat([
            _interpolate_retention(
                additional_injections, targets,
                vn_column="VN_cofm_mL_g", retention_mode="cofm",
                extrapolate=extrapolate, allow_invalid_retention=True,
            ),
            _interpolate_retention(
                additional_injections, targets,
                vn_column="VN_peak_max_mL_g", retention_mode="peak_max",
                extrapolate=extrapolate, allow_invalid_retention=True,
            ),
        ], ignore_index=True)
    injection_table = homologous_injections

    primary = _fit_profile(cofm, temperature_K)
    secondary = _fit_profile(peak_max, temperature_K)
    secondary = secondary[[
        "coverage", "gamma_d_mJm2", "r_squared", "slope_Jmol"
    ]].rename(columns={
        "gamma_d_mJm2": "gamma_d_pm_mJm2",
        "r_squared": "r_squared_pm",
        "slope_Jmol": "slope_pm_Jmol",
    })
    profile = primary.merge(secondary, on="coverage", how="left")
    profile["delta_cofm_pm"] = (
        profile["gamma_d_mJm2"] - profile["gamma_d_pm_mJm2"]
    )
    profile["W_cohesion_d_mJm2"] = 2.0 * profile["gamma_d_mJm2"]
    profile["W_cohesion_d_pm_mJm2"] = 2.0 * profile["gamma_d_pm_mJm2"]

    qc_input = cofm.rename(columns={
        "target_coverage": "target_coverage",
        "VN_mL_g": "VN_cofm",
    })
    qc = run_qc_checks(
        profile,
        qc_input,
        bounds=GAMMA_D_WARNING_BOUNDS_MJ_M2,
        critical_bounds=GAMMA_D_CRITICAL_BOUNDS_MJ_M2,
    )
    extrapolated = cofm[cofm["interpolation_status"].str.startswith("extrapolated")]
    if not extrapolated.empty:
        _append_flag(qc, {
            "check": "coverage_extrapolation",
            "severity": "warning",
            "coverage": None,
            "message": (
                f"{len(extrapolated)} homologous-series values were linearly "
                "extrapolated beyond their probe-specific measured coverage range"
            ),
            "value": int(len(extrapolated)),
        })
    outside_range = cofm[
        cofm["interpolation_status"] == "outside_measured_range"
    ]
    if not outside_range.empty:
        _append_flag(qc, {
            "check": "coverage_outside_measured_range",
            "severity": "warning",
            "coverage": None,
            "message": (
                f"{len(outside_range)} homologous-series values were left "
                "undefined because target coverage was outside the measured range"
            ),
            "value": int(len(outside_range)),
        })
    distinct_gains = sorted(set(detector_gains.tolist()))
    gain_variation = len(distinct_gains) > 1
    if gain_variation:
        _append_flag(qc, {
            "check": "detector_gain_variation",
            "severity": "warning",
            "coverage": None,
            "message": (
                "required injections use multiple detector gains; confirm the "
                "declared area-to-amount calibrations are valid across them; "
                "the profile is non-reportable until that review is resolved"
            ),
            "value": distinct_gains,
        })
    if dead_clipped or clipped_probe_ids:
        _append_flag(qc, {
            "check": "detector_clipping",
            "severity": "critical",
            "coverage": None,
            "message": "clipping was observed in one or more required traces",
            "value": sorted(dead_clipped + clipped_probe_ids),
        })
    _finish_qc(qc)

    complete_three_probe_fits = bool(
        len(profile)
        and (profile["n_alkanes"] >= 3).all()
        and np.isfinite(profile["gamma_d_mJm2"]).all()
        and (profile["slope_Jmol"] > 0).all()
    )
    no_extrapolation = bool(
        (interpolated["interpolation_status"] == "interpolated").all()
    )
    reportable = bool(
        qc["pass"]
        and complete_three_probe_fits
        and no_extrapolation
        and not gain_variation
    )

    return DispersiveResult(
        dataset_id=bundle.dataset_id,
        sample_id=str(experiment_row["sample_id"]),
        specific_surface_area_m2_g=ssa_m2_g,
        surface_area_source=ssa_source,
        temperature_K=temperature_K,
        dead_time_cofm_min=t0_cofm,
        dead_time_peak_max_min=t0_max,
        injections=injection_table,
        interpolated=interpolated,
        additional_injections=additional_injections,
        additional_interpolated=additional_interpolated,
        gamma_d=profile,
        qc=qc,
        reportable=reportable,
        flow_source_channels=tuple(sorted(flow_channels)),
        pressure_basis=tuple(sorted(pressure_bases)),
        pressure_roles=tuple(sorted(pressure_roles)),
        properties_sources=tuple(sorted(property_sources)),
        calibration_sources=tuple(sorted(calibration_sources)),
        detector_gains=tuple(sorted(set(detector_gains.tolist()))),
    )
