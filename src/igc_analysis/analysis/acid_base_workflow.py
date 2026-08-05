"""Coverage-resolved Gutmann acid/base workflow for neutral-data bundles.

Probe roles and properties are declared by the caller and the neutral bundle.
No probe identity or property value is inferred from a chemical name.
Center-of-mass retention is primary; peak maximum is a sensitivity result.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from igc_analysis.analysis.acid_base import gutmann_ka_kb
from igc_analysis.analysis.dispersive_workflow import (
    GAMMA_D_CRITICAL_BOUNDS_MJ_M2,
    GAMMA_D_WARNING_BOUNDS_MJ_M2,
    DispersiveResult,
    run_dispersive_from_neutral,
)
from igc_analysis.constants import N_AVOGADRO, R_GAS
from igc_analysis.io.neutral_data import NeutralBundle, read_neutral_bundle


GAMMA_L_D_PLAUSIBLE_BOUNDS_MJ_M2 = (1.0, 100.0)
SCHULTZ_R2_WARNING = 0.98
SCHULTZ_R2_CRITICAL = 0.95
DELTA_G_COLUMNS = (
    "retention_mode", "coverage", "probe", "probe_id", "probe_name",
    "VN_mL_g", "RT_ln_VN_Jmol", "x_schultz", "alkane_predicted_Jmol",
    "delta_g_sp_Jmol", "delta_g_sp_kJmol", "dn", "an_star",
    "properties_source",
)


@dataclass(frozen=True)
class AcidBaseResult:
    """Calculation tables, QC, reportability, and declared provenance."""

    dataset_id: str
    sample_id: str
    temperature_K: float
    profile: pd.DataFrame
    delta_g_sp: pd.DataFrame
    schultz_lines: pd.DataFrame
    injections: pd.DataFrame
    interpolated: pd.DataFrame
    qc: dict
    reportable: bool
    properties_sources: tuple[str, ...]
    calibration_sources: tuple[str, ...]
    flow_source_channels: tuple[str, ...]
    pressure_basis: tuple[str, ...]
    pressure_roles: tuple[str, ...]
    detector_gains: tuple[float, ...]
    dispersive: DispersiveResult


def _finite_positive(row: pd.Series, field: str, probe_id: str) -> float:
    raw = str(row[field]).strip()
    if not raw:
        raise ValueError(f"probe {probe_id!r}: {field} is required")
    value = float(raw)
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f"probe {probe_id!r}: {field} must be finite and positive")
    return value


def _finite_nonnegative(row: pd.Series, field: str, probe_id: str) -> float:
    raw = str(row[field]).strip()
    if not raw:
        raise ValueError(f"probe {probe_id!r}: {field} is required")
    value = float(raw)
    if not math.isfinite(value) or value < 0:
        raise ValueError(
            f"probe {probe_id!r}: {field} must be finite and nonnegative"
        )
    return value


def _property_records(
    bundle: NeutralBundle,
    homologous_ids: tuple[str, ...],
    polar_ids: tuple[str, ...],
) -> pd.DataFrame:
    properties = bundle.table("probe_properties.csv").copy()
    properties["probe_id"] = properties["probe_id"].astype(str)
    selected = properties[
        properties["probe_id"].isin(set(homologous_ids) | set(polar_ids))
    ].copy()
    if selected["probe_id"].nunique() != len(set(homologous_ids) | set(polar_ids)):
        raise ValueError("every selected probe requires one declared property record")
    if selected["probe_id"].duplicated().any():
        raise ValueError("selected probe property records must be unique")
    selected = selected.set_index("probe_id", drop=False)
    for probe_id in homologous_ids + polar_ids:
        row = selected.loc[probe_id]
        if not str(row["properties_source"]).strip():
            raise ValueError(f"probe {probe_id!r}: properties_source is required")
        _finite_positive(row, "cross_section_m2", probe_id)
        gamma_l_d = _finite_positive(row, "gamma_l_d_mJ_m2", probe_id)
        if not (
            GAMMA_L_D_PLAUSIBLE_BOUNDS_MJ_M2[0]
            <= gamma_l_d
            <= GAMMA_L_D_PLAUSIBLE_BOUNDS_MJ_M2[1]
        ):
            raise ValueError(
                f"probe {probe_id!r}: gamma_l_d_mJ_m2={gamma_l_d:g} is outside "
                "the 1 to 100 mJ/m2 unit/plausibility gate"
            )
    for probe_id in polar_ids:
        row = selected.loc[probe_id]
        _finite_nonnegative(row, "donor_number_kJ_mol", probe_id)
        _finite_positive(row, "acceptor_number_kJ_mol", probe_id)
    return selected


def _mode_profile(
    interpolated: pd.DataFrame,
    properties: pd.DataFrame,
    homologous_ids: tuple[str, ...],
    polar_ids: tuple[str, ...],
    temperature_K: float,
    retention_mode: str,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    mode = interpolated[
        interpolated["retention_mode"] == retention_mode
    ].copy()
    profile_rows: list[dict] = []
    dg_rows: list[dict] = []
    line_rows: list[dict] = []
    for coverage in sorted(mode["target_coverage"].unique()):
        at_coverage = mode[np.isclose(mode["target_coverage"], coverage)]
        alkane = at_coverage[
            at_coverage["probe_id"].isin(homologous_ids)
        ].dropna(subset=["VN_mL_g"])
        points: list[tuple[float, float]] = []
        for row in alkane.itertuples(index=False):
            prop = properties.loc[str(row.probe_id)]
            x_value = float(prop["cross_section_m2"]) * math.sqrt(
                float(prop["gamma_l_d_mJ_m2"]) * 1e-3
            )
            y_value = R_GAS * temperature_K * math.log(float(row.VN_mL_g))
            points.append((x_value, y_value))
        if len(points) >= 2:
            x_values = np.array([value[0] for value in points])
            y_values = np.array([value[1] for value in points])
            fit = stats.linregress(x_values, y_values)
            slope = float(fit.slope)
            intercept = float(fit.intercept)
            line_r2 = float(fit.rvalue**2)
            gamma_d = (slope / (2.0 * N_AVOGADRO)) ** 2 * 1e3
        else:
            slope = intercept = line_r2 = gamma_d = float("nan")
        line_rows.append({
            "retention_mode": retention_mode,
            "coverage": float(coverage),
            "slope": slope,
            "intercept": intercept,
            "r_squared": line_r2,
            "n_homologs": len(points),
            "gamma_d_schultz_mJm2": gamma_d,
        })

        polar = at_coverage[
            at_coverage["probe_id"].isin(polar_ids)
        ].dropna(subset=["VN_mL_g"])
        coverage_dg: list[dict] = []
        if math.isfinite(slope):
            for row in polar.itertuples(index=False):
                prop = properties.loc[str(row.probe_id)]
                x_value = float(prop["cross_section_m2"]) * math.sqrt(
                    float(prop["gamma_l_d_mJ_m2"]) * 1e-3
                )
                observed = R_GAS * temperature_K * math.log(float(row.VN_mL_g))
                predicted = slope * x_value + intercept
                delta = observed - predicted
                record = {
                    "retention_mode": retention_mode,
                    "coverage": float(coverage),
                    "probe": str(row.probe_id),
                    "probe_id": str(row.probe_id),
                    "probe_name": str(prop["probe_name"]),
                    "VN_mL_g": float(row.VN_mL_g),
                    "RT_ln_VN_Jmol": observed,
                    "x_schultz": x_value,
                    "alkane_predicted_Jmol": predicted,
                    "delta_g_sp_Jmol": delta,
                    "delta_g_sp_kJmol": delta / 1000.0,
                    "dn": float(prop["donor_number_kJ_mol"]),
                    "an_star": float(prop["acceptor_number_kJ_mol"]),
                    "properties_source": str(prop["properties_source"]),
                }
                coverage_dg.append(record)
                dg_rows.append(record)
        dg_frame = pd.DataFrame(coverage_dg)
        gutmann = gutmann_ka_kb(dg_frame)
        profile_rows.append({
            "retention_mode": retention_mode,
            "coverage": float(coverage),
            "Ka": gutmann["Ka"],
            "Kb": gutmann["Kb"],
            "Kb_Ka_ratio": gutmann["Kb_Ka_ratio"],
            "r_squared": gutmann["r_squared"],
            "fit_method": gutmann["fit_method"],
            "n_polar_probes": gutmann["n_probes"],
            "probe_ids_used": ", ".join(gutmann["probes_used"]),
            "schultz_r_squared": line_r2,
            "n_homologs": len(points),
            "gamma_d_schultz_mJm2": gamma_d,
        })
    return (
        pd.DataFrame(profile_rows),
        pd.DataFrame(dg_rows, columns=DELTA_G_COLUMNS),
        pd.DataFrame(line_rows),
    )


def _quality(
    primary: pd.DataFrame,
    delta_g_sp: pd.DataFrame,
    dispersive: DispersiveResult,
    interpolated: pd.DataFrame,
    injections: pd.DataFrame,
) -> tuple[dict, bool]:
    flags = [dict(flag) for flag in dispersive.qc["flags"]]
    for row in primary.itertuples(index=False):
        if row.n_homologs < 3:
            flags.append({"check": "schultz_probe_count", "severity": "critical",
                          "coverage": row.coverage,
                          "message": "fewer than three homologs define the reference line"})
        if math.isfinite(row.schultz_r_squared):
            if row.schultz_r_squared < SCHULTZ_R2_CRITICAL:
                flags.append({
                    "check": "schultz_r_squared", "severity": "critical",
                    "coverage": row.coverage,
                    "message": (
                        f"Schultz reference R2={row.schultz_r_squared:.3f} is "
                        f"below {SCHULTZ_R2_CRITICAL:.2f}"
                    ),
                })
            elif row.schultz_r_squared < SCHULTZ_R2_WARNING:
                flags.append({
                    "check": "schultz_r_squared", "severity": "warning",
                    "coverage": row.coverage,
                    "message": (
                        f"Schultz reference R2={row.schultz_r_squared:.3f} is "
                        f"below {SCHULTZ_R2_WARNING:.2f}"
                    ),
                })
        if math.isfinite(row.gamma_d_schultz_mJm2):
            if not (
                GAMMA_D_CRITICAL_BOUNDS_MJ_M2[0]
                <= row.gamma_d_schultz_mJm2
                <= GAMMA_D_CRITICAL_BOUNDS_MJ_M2[1]
            ):
                flags.append({
                    "check": "schultz_gamma_d_bounds", "severity": "critical",
                    "coverage": row.coverage,
                    "message": "Schultz-derived gamma_d is outside critical bounds",
                })
            elif not (
                GAMMA_D_WARNING_BOUNDS_MJ_M2[0]
                <= row.gamma_d_schultz_mJm2
                <= GAMMA_D_WARNING_BOUNDS_MJ_M2[1]
            ):
                flags.append({
                    "check": "schultz_gamma_d_bounds", "severity": "warning",
                    "coverage": row.coverage,
                    "message": "Schultz-derived gamma_d is outside expected bounds",
                })
        if row.n_polar_probes < 3 or row.fit_method != "regression":
            flags.append({"check": "gutmann_probe_count", "severity": "critical",
                          "coverage": row.coverage,
                          "message": "at least three polar probes are required for regression"})
        if math.isfinite(row.r_squared) and row.r_squared < 0.3:
            flags.append({"check": "gutmann_r_squared", "severity": "critical",
                          "coverage": row.coverage,
                          "message": f"Gutmann R2={row.r_squared:.3f} is below 0.3"})
        elif math.isfinite(row.r_squared) and row.r_squared < 0.5:
            flags.append({"check": "gutmann_r_squared", "severity": "warning",
                          "coverage": row.coverage,
                          "message": f"Gutmann R2={row.r_squared:.3f} is below 0.5"})
        if math.isfinite(row.Ka) and row.Ka < 0:
            flags.append({"check": "negative_Ka", "severity": "warning",
                          "coverage": row.coverage, "message": "Ka is negative"})
        if math.isfinite(row.Kb) and row.Kb < 0:
            flags.append({"check": "negative_Kb", "severity": "warning",
                          "coverage": row.coverage, "message": "Kb is negative"})
    for row in delta_g_sp[
        (delta_g_sp["retention_mode"] == "cofm")
        & (delta_g_sp["delta_g_sp_kJmol"] < 0)
    ].itertuples(index=False):
        flags.append({"check": "negative_delta_g_sp", "severity": "warning",
                      "coverage": row.coverage,
                      "message": f"probe {row.probe_id!r} lies below the reference line"})
    outside = interpolated[
        interpolated["interpolation_status"] != "interpolated"
    ]
    invalid_injections = int(
        (
            ~np.isfinite(injections["VN_cofm_mL_g"])
            | ~np.isfinite(injections["VN_peak_max_mL_g"])
        ).sum()
    )
    if invalid_injections:
        flags.append({
            "check": "positive_net_retention",
            "severity": "critical",
            "coverage": None,
            "message": (
                f"{invalid_injections} selected-probe injection(s) have "
                "nonpositive net retention"
            ),
        })
    if not outside.empty:
        flags.append({"check": "coverage_mapping", "severity": "warning",
                      "coverage": None,
                      "message": f"{len(outside)} selected-probe values are outside measured coverage ranges"})
    n_critical = sum(flag["severity"] == "critical" for flag in flags)
    n_warning = sum(flag["severity"] == "warning" for flag in flags)
    complete = bool(
        len(primary)
        and (primary["n_homologs"] >= 3).all()
        and (primary["n_polar_probes"] >= 3).all()
        and (primary["fit_method"] == "regression").all()
        and np.isfinite(primary[["Ka", "Kb", "r_squared"]]).all().all()
    )
    reportable = bool(
        n_critical == 0
        and complete
        and outside.empty
        and len(dispersive.detector_gains) == 1
    )
    status = "PASS" if n_critical == 0 else "FAIL"
    return {
        "pass": n_critical == 0,
        "flags": flags,
        "summary": f"QC {status}: {n_warning} warning(s), {n_critical} critical",
    }, reportable


def run_acid_base_from_neutral(
    bundle_dir: str | Path | NeutralBundle,
    *,
    homologous_probe_ids: tuple[str, ...] | list[str],
    polar_probe_ids: tuple[str, ...] | list[str],
    pressure_correction: bool = True,
    ambient_pressure_pa: float = 101325.0,
    extrapolate: bool = True,
    max_temperature_span_K: float = 1.0,
    max_flow_relative_span: float = 0.05,
) -> AcidBaseResult:
    """Run the declared Schultz--Gutmann workflow on a neutral bundle."""

    homologous_ids = tuple(str(value).strip() for value in homologous_probe_ids)
    polar_ids = tuple(str(value).strip() for value in polar_probe_ids)
    if len(homologous_ids) < 3 or len(set(homologous_ids)) != len(homologous_ids):
        raise ValueError("declare at least three unique homologous probe IDs")
    if len(polar_ids) < 3 or len(set(polar_ids)) != len(polar_ids):
        raise ValueError("declare at least three unique polar probe IDs")
    if any(not value for value in homologous_ids + polar_ids):
        raise ValueError("selected probe IDs cannot be empty")
    if set(homologous_ids) & set(polar_ids):
        raise ValueError("homologous and polar probe selections must be disjoint")
    bundle = bundle_dir if isinstance(bundle_dir, NeutralBundle) else read_neutral_bundle(bundle_dir)
    properties = _property_records(bundle, homologous_ids, polar_ids)
    dispersive = run_dispersive_from_neutral(
        bundle,
        homologous_probe_ids=homologous_ids,
        additional_probe_ids=polar_ids,
        pressure_correction=pressure_correction,
        ambient_pressure_pa=ambient_pressure_pa,
        extrapolate=extrapolate,
        max_temperature_span_K=max_temperature_span_K,
        max_flow_relative_span=max_flow_relative_span,
    )
    interpolated = pd.concat([
        dispersive.interpolated, dispersive.additional_interpolated
    ], ignore_index=True)
    profiles = []
    delta_tables = []
    line_tables = []
    for mode in ("cofm", "peak_max"):
        profile, delta, lines = _mode_profile(
            interpolated, properties, homologous_ids, polar_ids,
            dispersive.temperature_K, mode,
        )
        profiles.append(profile)
        delta_tables.append(delta)
        line_tables.append(lines)
    profile = pd.concat(profiles, ignore_index=True)
    delta_g_sp = pd.concat(delta_tables, ignore_index=True)
    schultz_lines = pd.concat(line_tables, ignore_index=True)
    primary = profile[profile["retention_mode"] == "cofm"]
    all_injections = pd.concat([
        dispersive.injections, dispersive.additional_injections
    ], ignore_index=True).sort_values("sequence_index")
    qc, reportable = _quality(
        primary, delta_g_sp, dispersive, interpolated, all_injections
    )
    return AcidBaseResult(
        dataset_id=bundle.dataset_id,
        sample_id=dispersive.sample_id,
        temperature_K=dispersive.temperature_K,
        profile=profile,
        delta_g_sp=delta_g_sp,
        schultz_lines=schultz_lines,
        injections=all_injections,
        interpolated=interpolated,
        qc=qc,
        reportable=reportable,
        properties_sources=dispersive.properties_sources,
        calibration_sources=dispersive.calibration_sources,
        flow_source_channels=dispersive.flow_source_channels,
        pressure_basis=dispersive.pressure_basis,
        pressure_roles=dispersive.pressure_roles,
        detector_gains=dispersive.detector_gains,
        dispersive=dispersive,
    )
