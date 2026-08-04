"""Generate the deterministic closed-form synthetic BET contract fixture."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = (
    ROOT / "src" / "igc_analysis" / "contracts" / "igc-neutral-data" / "0.2.0"
    / "examples" / "synthetic_bet_isotherm"
)

DATASET_ID = "synthetic-bet-isotherm-001"
EXPERIMENT_ID = "experiment-bet-001"
SAMPLE_ID = "sample-bet-synthetic-001"
COLUMN_ID = "column-bet-001"
BLOCK_ID = "block-bet-001"
PROBE_ID = "probe-bet-synthetic"
DEAD_TIME_PROBE_ID = "probe-dead-time-synthetic"
CALIBRATION_ID = "calibration-bet-synthetic"
DETECTOR_CHANNEL = "detector-primary"

TEMPERATURE_K = 303.15
STANDARD_TEMPERATURE_K = 273.15
FLOW_COLUMN_M3_S = 1.85e-7
FLOW_COLUMN_ML_MIN = FLOW_COLUMN_M3_S * 60.0 * 1e6
PRESSURE_PA = 101325.0
P_SAT_PA = 12000.0
SAMPLE_MASS_G = 0.25
CROSS_SECTION_M2 = 6.3e-19
CALIBRATION_MOL_PER_AREA = 1.0e-6
R_GAS = 8.314462618
N_AVOGADRO = 6.02214076e23
T0_MIN = 0.4
TRACE_END_MIN = 2.5
TRACE_STEP_MIN = 0.002
PEAK_SIGMA_MIN = 0.025
DEAD_TIME_SIGMA_MIN = 0.015

TARGET_N_MONOLAYER_MMOL_G = 0.01
TARGET_C_BET = 20.0
TARGET_SSA_M2_G = (
    TARGET_N_MONOLAYER_MMOL_G * 1e-3 * N_AVOGADRO * CROSS_SECTION_M2
)
RELATIVE_PRESSURES = (
    0.001, 0.02, 0.04, 0.06, 0.08, 0.10, 0.13,
    0.16, 0.20, 0.25, 0.30, 0.35, 0.40,
)

FILE_ORDER = (
    "experiment.csv",
    "columns.csv",
    "conditions.csv",
    "injections.csv",
    "injection_components.csv",
    "traces.csv",
    "probe_properties.csv",
    "calibration.csv",
)


def _write_csv(
    output: Path, filename: str, fieldnames: list[str], rows: list[dict]
) -> None:
    path = output / filename
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _bet_q_mol_g(x: float) -> float:
    n_m = TARGET_N_MONOLAYER_MMOL_G * 1e-3
    denominator = (1.0 - x) * (1.0 + (TARGET_C_BET - 1.0) * x)
    return n_m * TARGET_C_BET * x / denominator


def _bet_dqdx_mol_g(x: float) -> float:
    n_m = TARGET_N_MONOLAYER_MMOL_G * 1e-3
    denominator = 1.0 + (TARGET_C_BET - 2.0) * x - (TARGET_C_BET - 1.0) * x**2
    derivative = (TARGET_C_BET - 2.0) - 2.0 * (TARGET_C_BET - 1.0) * x
    return (
        n_m * TARGET_C_BET * (denominator - x * derivative) / denominator**2
    )


def _probe_spec(relative_pressure: float) -> tuple[float, float, float]:
    concentration = relative_pressure * P_SAT_PA / (R_GAS * TEMPERATURE_K)
    apex = concentration * FLOW_COLUMN_M3_S * 60.0 / CALIBRATION_MOL_PER_AREA
    dc_dx = P_SAT_PA / (R_GAS * TEMPERATURE_K)
    retention_volume_m3 = SAMPLE_MASS_G * _bet_dqdx_mol_g(relative_pressure) / dc_dx
    retention_time_min = T0_MIN + retention_volume_m3 * 1e6 / FLOW_COLUMN_ML_MIN
    retention_time_min = round(retention_time_min / TRACE_STEP_MIN) * TRACE_STEP_MIN
    return concentration, apex, retention_time_min


def _gaussian(time_min: float, center_min: float, amplitude: float, sigma: float) -> float:
    z = (time_min - center_min) / sigma
    return amplitude * math.exp(-0.5 * z * z)


def generate(output: Path = DEFAULT_OUTPUT) -> None:
    output.mkdir(parents=True, exist_ok=True)

    injection_rows: list[dict] = []
    component_rows: list[dict] = []
    condition_rows: list[dict] = []
    trace_rows: list[dict] = []
    trace_specs: list[tuple[str, float, float, float]] = []

    ids = ["dead-time-000"] + [
        f"probe-injection-{index:03d}"
        for index in range(1, len(RELATIVE_PRESSURES) + 1)
    ] + ["dead-time-999"]

    for sequence_index, injection_id in enumerate(ids):
        is_dead_time = injection_id.startswith("dead-time")
        acquired_minute = sequence_index * 5
        injection_rows.append({
            "experiment_id": EXPERIMENT_ID,
            "injection_id": injection_id,
            "block_id": BLOCK_ID,
            "sequence_index": sequence_index,
            "acquired_at": f"2026-01-20T{9 + acquired_minute // 60:02d}:{acquired_minute % 60:02d}:00-05:00",
            "role": "dead_time" if is_dead_time else "probe",
            "target_coverage_fraction": "" if is_dead_time else RELATIVE_PRESSURES[sequence_index - 1],
            "detector_gain": 1.0,
            "detector_channel": DETECTOR_CHANNEL,
            "clipping_observed": "false",
        })

        if is_dead_time:
            component_rows.append({
                "injection_id": injection_id,
                "component_index": 0,
                "probe_id": DEAD_TIME_PROBE_ID,
                "component_role": "dead_time_marker",
                "target_amount_mol": 1e-8,
                "calibration_id": "",
                "saturation_vapor_pressure_Pa": "",
                "vapor_pressure_source": "",
                "vapor_pressure_model_id": "",
            })
            trace_specs.append((injection_id, T0_MIN, 5.0, DEAD_TIME_SIGMA_MIN))
        else:
            relative_pressure = RELATIVE_PRESSURES[sequence_index - 1]
            _, apex, retention_time_min = _probe_spec(relative_pressure)
            target_amount = (
                apex * math.sqrt(2.0 * math.pi) * PEAK_SIGMA_MIN
                * CALIBRATION_MOL_PER_AREA
            )
            component_rows.append({
                "injection_id": injection_id,
                "component_index": 0,
                "probe_id": PROBE_ID,
                "component_role": "analyte",
                "target_amount_mol": f"{target_amount:.12g}",
                "calibration_id": CALIBRATION_ID,
                "saturation_vapor_pressure_Pa": P_SAT_PA,
                "vapor_pressure_source": "synthetic-closed-form-bet",
                "vapor_pressure_model_id": "synthetic-constant-v1",
            })
            trace_specs.append((
                injection_id, retention_time_min, apex, PEAK_SIGMA_MIN
            ))

        for quantity, value, unit, source_channel in (
            ("column_temperature", TEMPERATURE_K, "K", ""),
            ("flow_column", FLOW_COLUMN_M3_S, "m3_s", "flow-channel-1"),
            ("pressure_inlet", PRESSURE_PA, "Pa", ""),
            ("pressure_outlet", PRESSURE_PA, "Pa", ""),
            ("pressure_drop", 0.0, "Pa", ""),
        ):
            condition_rows.append({
                "condition_id": f"condition-{len(condition_rows):05d}",
                "injection_id": injection_id,
                "quantity": quantity,
                "value": value,
                "unit": unit,
                "value_role": "measured",
                "measured_at": injection_rows[-1]["acquired_at"],
                "measurement_basis": "direct",
                "source_channel": source_channel,
            })

    n_points = int(round(TRACE_END_MIN / TRACE_STEP_MIN)) + 1
    for injection_id, center, amplitude, sigma in trace_specs:
        for point_index in range(n_points):
            time_min = point_index * TRACE_STEP_MIN
            signal = _gaussian(time_min, center, amplitude, sigma)
            trace_rows.append({
                "injection_id": injection_id,
                "point_index": point_index,
                "time_s": f"{time_min * 60.0:.6f}",
                "detector_channel": DETECTOR_CHANNEL,
                "signal_raw": f"{signal:.12g}",
                "signal_unit": "arbitrary_unit",
                "signal_corrected": "",
                "preprocessing_method": "",
                "preprocessing_version": "",
            })

    _write_csv(output, "experiment.csv", [
        "dataset_id", "experiment_id", "sample_id", "sample_mass_g",
        "injection_loop_volume_m3", "column_id", "acquisition_started_at",
        "specific_surface_area_m2_g", "surface_area_source",
    ], [{
        "dataset_id": DATASET_ID,
        "experiment_id": EXPERIMENT_ID,
        "sample_id": SAMPLE_ID,
        "sample_mass_g": SAMPLE_MASS_G,
        "injection_loop_volume_m3": 2.0e-6,
        "column_id": COLUMN_ID,
        "acquisition_started_at": "2026-01-20T09:00:00-05:00",
        "specific_surface_area_m2_g": "",
        "surface_area_source": "",
    }])
    _write_csv(output, "columns.csv", [
        "column_id", "packing_replicate_id", "column_role",
        "internal_diameter_m", "packed_bed_length_m",
        "conditioning_description", "sample_batch_id", "density_kg_m3",
        "density_basis",
    ], [{
        "column_id": COLUMN_ID,
        "packing_replicate_id": "packing-bet-001",
        "column_role": "sample",
        "internal_diameter_m": 0.004,
        "packed_bed_length_m": 0.3,
        "conditioning_description": "synthetic-protocol",
        "sample_batch_id": "batch-bet-synthetic-001",
        "density_kg_m3": 450,
        "density_basis": "bulk_packed",
    }])
    _write_csv(output, "conditions.csv", [
        "condition_id", "injection_id", "quantity", "value", "unit",
        "value_role", "measured_at", "measurement_basis", "source_channel",
    ], condition_rows)
    _write_csv(output, "injections.csv", [
        "experiment_id", "injection_id", "block_id", "sequence_index",
        "acquired_at", "role", "target_coverage_fraction", "detector_gain",
        "detector_channel", "clipping_observed",
    ], injection_rows)
    _write_csv(output, "injection_components.csv", [
        "injection_id", "component_index", "probe_id", "component_role",
        "target_amount_mol", "calibration_id", "saturation_vapor_pressure_Pa",
        "vapor_pressure_source", "vapor_pressure_model_id",
    ], component_rows)
    _write_csv(output, "traces.csv", [
        "injection_id", "point_index", "time_s", "detector_channel",
        "signal_raw", "signal_unit", "signal_corrected",
        "preprocessing_method", "preprocessing_version",
    ], trace_rows)
    _write_csv(output, "probe_properties.csv", [
        "probe_id", "probe_name", "molar_mass_g_mol", "cross_section_m2",
        "gamma_l_d_mJ_m2", "donor_number_kJ_mol",
        "acceptor_number_kJ_mol", "carbon_number", "properties_source",
    ], [
        {
            "probe_id": DEAD_TIME_PROBE_ID,
            "probe_name": "synthetic dead-time marker",
            "molar_mass_g_mol": 16.0,
            "cross_section_m2": 2.0e-19,
            "gamma_l_d_mJ_m2": "",
            "donor_number_kJ_mol": "",
            "acceptor_number_kJ_mol": "",
            "carbon_number": "",
            "properties_source": "synthetic-closed-form-bet",
        },
        {
            "probe_id": PROBE_ID,
            "probe_name": "synthetic BET probe",
            "molar_mass_g_mol": 114.0,
            "cross_section_m2": CROSS_SECTION_M2,
            "gamma_l_d_mJ_m2": 21.0,
            "donor_number_kJ_mol": "",
            "acceptor_number_kJ_mol": "",
            "carbon_number": 8,
            "properties_source": "synthetic-closed-form-bet",
        },
    ])
    _write_csv(output, "calibration.csv", [
        "calibration_id", "probe_id", "calibration_model", "parameter_0",
        "parameter_1", "parameter_2", "area_unit", "amount_unit",
        "calibration_source",
    ], [{
        "calibration_id": CALIBRATION_ID,
        "probe_id": PROBE_ID,
        "calibration_model": "linear",
        "parameter_0": 0,
        "parameter_1": CALIBRATION_MOL_PER_AREA,
        "parameter_2": "",
        "area_unit": "arbitrary_unit_min",
        "amount_unit": "mol",
        "calibration_source": "synthetic-closed-form-bet",
    }])

    files = {}
    for filename in FILE_ORDER:
        path = output / filename
        with path.open(encoding="utf-8", newline="") as handle:
            row_count = sum(1 for _ in handle) - 1
        files[filename] = {
            "row_count": row_count,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        }
    manifest = {
        "contract_name": "igc-neutral-data",
        "contract_version": "0.2.0",
        "profile": "trace-core",
        "dataset_id": DATASET_ID,
        "created_at": "2026-01-20T12:00:00-05:00",
        "adapter_version": "synthetic-bet-generator-1.0.0",
        "source_fingerprint": hashlib.sha256(
            b"synthetic-closed-form-bet-v1"
        ).hexdigest(),
        "files": files,
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="fixture destination (default: packaged contract example)",
    )
    generate(parser.parse_args().output)
