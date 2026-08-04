"""Generate a deterministic closed-form dispersive-energy contract fixture."""

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
    / "examples" / "synthetic_dispersive_profile"
)

DATASET_ID = "synthetic-dispersive-profile-001"
EXPERIMENT_ID = "experiment-dispersive-001"
SAMPLE_ID = "sample-dispersive-synthetic-001"
COLUMN_ID = "column-dispersive-001"
BLOCK_ID = "block-dispersive-001"
DEAD_PROBE_ID = "probe-dead-time-synthetic"
DETECTOR_CHANNEL = "detector-primary"

TEMPERATURE_K = 303.15
FLOW_COLUMN_M3_S = 1.66666666666667e-7
FLOW_COLUMN_ML_MIN = FLOW_COLUMN_M3_S * 60.0 * 1e6
PRESSURE_OUTLET_PA = 101325.0
PRESSURE_INLET_PA = 105000.0
PRESSURE_DROP_PA = PRESSURE_INLET_PA - PRESSURE_OUTLET_PA
SAMPLE_MASS_G = 0.10
SSA_M2_G = 10.0
CALIBRATION_MOL_PER_AREA = 1.0e-8
T0_MIN = 0.40
TRACE_END_MIN = 3.0
TRACE_STEP_MIN = 0.002
PEAK_SIGMA_MIN = 0.020
DEAD_SIGMA_MIN = 0.012
TAIL_SIGMA_MIN = 0.050
TAIL_DELAY_MIN = 0.040
TAIL_FRACTION = 0.20
N_AVOGADRO = 6.02214076e23
R_GAS = 8.314462618
A_CH2_M2 = 6.0e-20
GAMMA_CH2_J_M2 = (35.6 - 0.058 * (TEMPERATURE_K - 293.15)) * 1e-3
TARGET_GAMMA_D_MJ_M2 = 40.0
TARGET_COVERAGES = (0.01, 0.02, 0.03, 0.04)

PROBES = (
    ("probe-homolog-08", "synthetic homolog C8", 8, 5.1e-19, -0.001),
    ("probe-homolog-09", "synthetic homolog C9", 9, 5.7e-19, 0.000),
    ("probe-homolog-10", "synthetic homolog C10", 10, 6.3e-19, 0.001),
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
    with (output / filename).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _gaussian(time_min: float, center_min: float, amplitude: float, sigma: float) -> float:
    z = (time_min - center_min) / sigma
    return amplitude * math.exp(-0.5 * z * z)


def _probe_signal(time_min: float, center_min: float, amplitude: float) -> float:
    return (
        (1.0 - TAIL_FRACTION)
        * _gaussian(time_min, center_min, amplitude, PEAK_SIGMA_MIN)
        + TAIL_FRACTION
        * _gaussian(
            time_min, center_min + TAIL_DELAY_MIN, amplitude, TAIL_SIGMA_MIN
        )
    )


def _probe_area_factor() -> float:
    return math.sqrt(2.0 * math.pi) * (
        (1.0 - TAIL_FRACTION) * PEAK_SIGMA_MIN
        + TAIL_FRACTION * TAIL_SIGMA_MIN
    )


def _probe_mean_offset() -> float:
    tail_area = TAIL_FRACTION * TAIL_SIGMA_MIN
    total_area = (
        (1.0 - TAIL_FRACTION) * PEAK_SIGMA_MIN + tail_area
    )
    return TAIL_DELAY_MIN * tail_area / total_area


def _james_martin_factor() -> float:
    ratio = PRESSURE_INLET_PA / PRESSURE_OUTLET_PA
    return 1.5 * (ratio**2 - 1.0) / (ratio**3 - 1.0)


def _dorris_slope_j_mol() -> float:
    gamma_d_j_m2 = TARGET_GAMMA_D_MJ_M2 * 1e-3
    return math.sqrt(
        gamma_d_j_m2
        * 4.0
        * N_AVOGADRO**2
        * A_CH2_M2**2
        * GAMMA_CH2_J_M2
    )


def _retention_volume_mL_g(carbon_number: int, coverage: float) -> float:
    # The coverage factor is linear, so the workflow's declared piecewise-linear
    # interpolation recovers the same homologous-series slope exactly.
    coverage_factor = 0.0045 - 0.010 * coverage
    return coverage_factor * math.exp(
        _dorris_slope_j_mol() * carbon_number / (R_GAS * TEMPERATURE_K)
    )


def generate(output: Path = DEFAULT_OUTPUT) -> None:
    output.mkdir(parents=True, exist_ok=True)
    injection_rows: list[dict] = []
    component_rows: list[dict] = []
    condition_rows: list[dict] = []
    trace_specs: list[tuple[str, float, float, float, bool]] = []

    def add_injection(
        injection_id: str,
        sequence_index: int,
        *,
        role: str,
        probe_id: str,
        target_coverage: float | None,
        actual_coverage: float | None = None,
        carbon_number: int | None = None,
        cross_section_m2: float | None = None,
    ) -> None:
        acquired_minute = sequence_index * 4
        acquired_at = (
            f"2026-01-21T{9 + acquired_minute // 60:02d}:"
            f"{acquired_minute % 60:02d}:00-05:00"
        )
        injection_rows.append({
            "experiment_id": EXPERIMENT_ID,
            "injection_id": injection_id,
            "block_id": BLOCK_ID,
            "sequence_index": sequence_index,
            "acquired_at": acquired_at,
            "role": role,
            "target_coverage_fraction": "" if target_coverage is None else target_coverage,
            "detector_gain": 1.0,
            "detector_channel": DETECTOR_CHANNEL,
            "clipping_observed": "false",
        })
        is_dead = role == "dead_time"
        calibration_id = "" if is_dead else f"calibration-{probe_id}"
        component_rows.append({
            "injection_id": injection_id,
            "component_index": 0,
            "probe_id": probe_id,
            "component_role": "dead_time_marker" if is_dead else "analyte",
            "target_amount_mol": 1e-9 if is_dead else "",
            "calibration_id": calibration_id,
            "saturation_vapor_pressure_Pa": "",
            "vapor_pressure_source": "",
            "vapor_pressure_model_id": "",
        })

        if is_dead:
            trace_specs.append((injection_id, T0_MIN, 8.0, DEAD_SIGMA_MIN, False))
        else:
            assert actual_coverage is not None
            assert carbon_number is not None
            assert cross_section_m2 is not None
            capacity_mol = (
                SSA_M2_G * SAMPLE_MASS_G
                / (N_AVOGADRO * cross_section_m2)
            )
            amount_mol = actual_coverage * capacity_mol
            area = amount_mol / CALIBRATION_MOL_PER_AREA
            amplitude = area / _probe_area_factor()
            vn_mL_g = _retention_volume_mL_g(carbon_number, actual_coverage)
            cofm_retention_min = (
                T0_MIN
                + vn_mL_g * SAMPLE_MASS_G
                / (FLOW_COLUMN_ML_MIN * _james_martin_factor())
            )
            center_min = cofm_retention_min - _probe_mean_offset()
            center_min = round(center_min / TRACE_STEP_MIN) * TRACE_STEP_MIN
            trace_specs.append((
                injection_id, center_min, amplitude, PEAK_SIGMA_MIN, True
            ))

        for quantity, value, unit, source_channel in (
            ("column_temperature", TEMPERATURE_K, "K", ""),
            ("flow_column", FLOW_COLUMN_M3_S, "m3_s", "flow-channel-synthetic"),
            ("pressure_inlet", PRESSURE_INLET_PA, "Pa", ""),
            ("pressure_outlet", PRESSURE_OUTLET_PA, "Pa", ""),
            ("pressure_drop", PRESSURE_DROP_PA, "Pa", ""),
        ):
            condition_rows.append({
                "condition_id": f"condition-{len(condition_rows):05d}",
                "injection_id": injection_id,
                "quantity": quantity,
                "value": value,
                "unit": unit,
                "value_role": "measured",
                "measured_at": acquired_at,
                "measurement_basis": "direct",
                "source_channel": source_channel,
            })

    add_injection(
        "dead-time-000", 0, role="dead_time", probe_id=DEAD_PROBE_ID,
        target_coverage=None,
    )
    sequence = 1
    for probe_id, _, carbon_number, cross_section_m2, offset in PROBES:
        for target in TARGET_COVERAGES:
            add_injection(
                f"injection-{carbon_number:02d}-{sequence:03d}",
                sequence,
                role="probe",
                probe_id=probe_id,
                target_coverage=target,
                actual_coverage=(
                    target if target in (TARGET_COVERAGES[0], TARGET_COVERAGES[-1])
                    else target + offset
                ),
                carbon_number=carbon_number,
                cross_section_m2=cross_section_m2,
            )
            sequence += 1
    add_injection(
        "dead-time-999", sequence, role="dead_time", probe_id=DEAD_PROBE_ID,
        target_coverage=None,
    )

    trace_rows: list[dict] = []
    n_points = int(round(TRACE_END_MIN / TRACE_STEP_MIN)) + 1
    for injection_id, center, amplitude, sigma, tailed in trace_specs:
        for point_index in range(n_points):
            time_min = point_index * TRACE_STEP_MIN
            signal = (
                _probe_signal(time_min, center, amplitude)
                if tailed
                else _gaussian(time_min, center, amplitude, sigma)
            )
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
        "acquisition_started_at": "2026-01-21T09:00:00-05:00",
        "specific_surface_area_m2_g": SSA_M2_G,
        "surface_area_source": "synthetic-declared-ssa-v1",
    }])
    _write_csv(output, "columns.csv", [
        "column_id", "packing_replicate_id", "column_role",
        "internal_diameter_m", "packed_bed_length_m",
        "conditioning_description", "sample_batch_id", "density_kg_m3",
        "density_basis",
    ], [{
        "column_id": COLUMN_ID,
        "packing_replicate_id": "packing-dispersive-001",
        "column_role": "sample",
        "internal_diameter_m": 0.004,
        "packed_bed_length_m": 0.3,
        "conditioning_description": "synthetic-protocol",
        "sample_batch_id": "batch-dispersive-synthetic-001",
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
    property_rows = [{
        "probe_id": DEAD_PROBE_ID,
        "probe_name": "synthetic dead-time marker",
        "molar_mass_g_mol": 16.0,
        "cross_section_m2": 2.0e-19,
        "gamma_l_d_mJ_m2": "",
        "donor_number_kJ_mol": "",
        "acceptor_number_kJ_mol": "",
        "carbon_number": "",
        "properties_source": "synthetic-reference-properties-v1",
    }]
    calibration_rows = []
    for probe_id, probe_name, carbon_number, cross_section_m2, _ in PROBES:
        property_rows.append({
            "probe_id": probe_id,
            "probe_name": probe_name,
            "molar_mass_g_mol": 14.0 * carbon_number + 2.0,
            "cross_section_m2": cross_section_m2,
            "gamma_l_d_mJ_m2": "",
            "donor_number_kJ_mol": "",
            "acceptor_number_kJ_mol": "",
            "carbon_number": carbon_number,
            "properties_source": "synthetic-reference-properties-v1",
        })
        calibration_rows.append({
            "calibration_id": f"calibration-{probe_id}",
            "probe_id": probe_id,
            "calibration_model": "linear",
            "parameter_0": 0,
            "parameter_1": CALIBRATION_MOL_PER_AREA,
            "parameter_2": "",
            "area_unit": "arbitrary_unit_min",
            "amount_unit": "mol",
            "calibration_source": "synthetic-calibration-v1",
        })
    _write_csv(output, "probe_properties.csv", [
        "probe_id", "probe_name", "molar_mass_g_mol", "cross_section_m2",
        "gamma_l_d_mJ_m2", "donor_number_kJ_mol",
        "acceptor_number_kJ_mol", "carbon_number", "properties_source",
    ], property_rows)
    _write_csv(output, "calibration.csv", [
        "calibration_id", "probe_id", "calibration_model", "parameter_0",
        "parameter_1", "parameter_2", "area_unit", "amount_unit",
        "calibration_source",
    ], calibration_rows)

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
        "created_at": "2026-01-21T12:00:00-05:00",
        "adapter_version": "synthetic-dispersive-generator-1.0.0",
        "source_fingerprint": hashlib.sha256(
            b"synthetic-closed-form-dispersive-v1"
        ).hexdigest(),
        "files": files,
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    generate(parser.parse_args().output)
