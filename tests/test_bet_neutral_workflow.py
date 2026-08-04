"""Regression tests for the source-neutral corrected BET workflow."""

from __future__ import annotations

import json
import hashlib
import shutil
from pathlib import Path

import pandas as pd
import numpy as np
import pytest

from igc_analysis.analysis.bet import james_martin_j
import igc_analysis.analysis.bet_workflow as bet_workflow
from igc_analysis.analysis.bet_workflow import _pressure_factor, run_bet_from_neutral
from igc_analysis.analysis.full_peak import build_trace_dataset_from_neutral
from igc_analysis.cli.bet import main
from igc_analysis import __version__
from igc_analysis.constants import R_GAS
from igc_analysis.io.neutral_data import bundled_contract_path


@pytest.fixture()
def synthetic_bundle() -> Path:
    return bundled_contract_path() / "examples" / "synthetic_peak_shape"


@pytest.fixture()
def synthetic_bet_bundle() -> Path:
    return bundled_contract_path() / "examples" / "synthetic_bet_isotherm"


def _copy_bundle(synthetic_bundle: Path, tmp_path: Path) -> Path:
    copied = tmp_path / "bundle"
    shutil.copytree(synthetic_bundle, copied)
    return copied


def _write_table_and_refresh_manifest(bundle: Path, filename: str, table: pd.DataFrame):
    path = bundle / filename
    table.to_csv(path, index=False)
    manifest_path = bundle / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["files"][filename]["row_count"] = len(table)
    manifest["files"][filename]["sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


def test_neutral_bet_preserves_corrected_concentration_flow_and_pressure(
    synthetic_bundle: Path,
):
    result = run_bet_from_neutral(synthetic_bundle)
    by_number = {item.injection_number: item for item in result.injections}
    first = by_number[1]
    assert first.injection_id == "injection-001"

    # Calibration declares area in signal-unit seconds, while the peak core
    # reports signal-unit minutes. The conversion must be explicit.
    assert first.n_injected_mol == pytest.approx(first.peak_area * 60.0e-9)

    # Corrected concentration is the eluted peak-apex concentration, never the
    # pre-injection loop concentration.
    expected_c = (
        first.peak_height
        * first.n_injected_mol
        / first.peak_area
        / (first.flow_col_mL_min * 1e-6)
    )
    assert first.concentration_mol_m3 == pytest.approx(expected_c)
    assert first.P_over_P0 == pytest.approx(
        expected_c * R_GAS * first.temp_col_K / 12000.0
    )

    # The declared measured standard flow is converted to column temperature;
    # it is not replaced by a default or a different channel.
    expected_column_flow = 1.67e-7 * first.temp_col_K / 273.15 * 60.0 * 1e6
    assert first.flow_col_mL_min == pytest.approx(expected_column_flow)

    # Direct declared absolute pressures drive the James--Martin factor.
    assert first.j_factor == pytest.approx(james_martin_j(104845.0, 101325.0))
    assert first.p_sat_Pa == pytest.approx(12000.0)


def test_pressure_correction_can_be_disabled_for_sensitivity_only(
    synthetic_bundle: Path,
):
    corrected = run_bet_from_neutral(synthetic_bundle)
    uncorrected = run_bet_from_neutral(
        synthetic_bundle, pressure_correction=False
    )
    assert all(item.j_factor < 1 for item in corrected.injections)
    assert all(item.j_factor == 1 for item in uncorrected.injections)
    assert all(
        uncorrected_item.V_N_mL > corrected_item.V_N_mL
        for corrected_item, uncorrected_item in zip(
            corrected.injections, uncorrected.injections
        )
    )


def test_pressure_correction_refuses_missing_or_inverted_pressure():
    missing = pd.DataFrame(columns=[
        "injection_id", "quantity", "value", "value_role"
    ])
    with pytest.raises(ValueError, match="requires pressure_inlet or pressure_drop"):
        _pressure_factor(
            missing, "injection-001", enabled=True, ambient_pressure_pa=101325.0
        )

    inverted = pd.DataFrame([
        {"injection_id": "injection-001", "quantity": "pressure_inlet",
         "value": 100000.0, "value_role": "measured"},
        {"injection_id": "injection-001", "quantity": "pressure_outlet",
         "value": 101325.0, "value_role": "measured"},
    ])
    with pytest.raises(ValueError, match="cannot be below"):
        _pressure_factor(
            inverted, "injection-001", enabled=True,
            ambient_pressure_pa=101325.0
        )


def test_pressure_correction_refuses_inconsistent_redundant_values():
    inconsistent = pd.DataFrame([
        {"injection_id": "injection-001", "quantity": "pressure_inlet",
         "value": 105000.0, "value_role": "measured"},
        {"injection_id": "injection-001", "quantity": "pressure_outlet",
         "value": 101325.0, "value_role": "measured"},
        {"injection_id": "injection-001", "quantity": "pressure_drop",
         "value": 1000.0, "value_role": "measured"},
    ])
    with pytest.raises(ValueError, match="are inconsistent"):
        _pressure_factor(
            inconsistent, "injection-001", enabled=True,
            ambient_pressure_pa=101325.0
        )


def test_pressure_drop_only_uses_explicit_ambient_pressure(
    synthetic_bundle: Path, tmp_path: Path,
):
    copied = _copy_bundle(synthetic_bundle, tmp_path)
    conditions = pd.read_csv(copied / "conditions.csv", keep_default_na=False)
    conditions = conditions[
        ~conditions["quantity"].isin(["pressure_inlet", "pressure_outlet"])
    ]
    _write_table_and_refresh_manifest(copied, "conditions.csv", conditions)

    default = run_bet_from_neutral(copied)
    higher_ambient = run_bet_from_neutral(copied, ambient_pressure_pa=110000.0)
    first_default = min(default.injections, key=lambda item: item.injection_number)
    first_higher = min(higher_ambient.injections, key=lambda item: item.injection_number)
    assert first_default.j_factor == pytest.approx(
        james_martin_j(101325.0 + 3520.0, 101325.0)
    )
    assert first_higher.j_factor > first_default.j_factor
    assert default.pressure_basis == ("declared_drop_plus_ambient_absolute_outlet",)


def test_mixed_condition_roles_are_reported(
    synthetic_bundle: Path, tmp_path: Path,
):
    copied = _copy_bundle(synthetic_bundle, tmp_path)
    conditions = pd.read_csv(copied / "conditions.csv", keep_default_na=False)
    probe_ids = {"injection-001", "injection-002"}
    drop_measured_temp = (
        conditions["injection_id"].isin(probe_ids)
        & (conditions["quantity"] == "column_temperature")
        & (conditions["value_role"] == "measured")
    )
    conditions = conditions[~drop_measured_temp]
    _write_table_and_refresh_manifest(copied, "conditions.csv", conditions)
    result = run_bet_from_neutral(copied)
    assert result.conditions_source == "mixed"


def test_ambiguous_measured_flow_channels_are_rejected(
    synthetic_bundle: Path, tmp_path: Path,
):
    copied = _copy_bundle(synthetic_bundle, tmp_path)
    conditions = pd.read_csv(copied / "conditions.csv", keep_default_na=False)
    source = conditions[
        (conditions["injection_id"] == "injection-001")
        & (conditions["quantity"] == "flow_standard")
        & (conditions["value_role"] == "measured")
    ].iloc[0].copy()
    source["condition_id"] = "condition-0012-alt"
    source["measurement_basis"] = "converted"
    source["source_channel"] = "flow-channel-2"
    conditions = pd.concat([conditions, pd.DataFrame([source])], ignore_index=True)
    _write_table_and_refresh_manifest(copied, "conditions.csv", conditions)
    with pytest.raises(ValueError, match="ambiguous measured flow_standard"):
        run_bet_from_neutral(copied)


def test_missing_saturation_pressure_is_rejected(
    synthetic_bundle: Path, tmp_path: Path,
):
    copied = _copy_bundle(synthetic_bundle, tmp_path)
    components = pd.read_csv(
        copied / "injection_components.csv", keep_default_na=False
    )
    mask = components["component_role"] == "analyte"
    components.loc[mask, [
        "saturation_vapor_pressure_Pa", "vapor_pressure_source",
        "vapor_pressure_model_id",
    ]] = ""
    _write_table_and_refresh_manifest(
        copied, "injection_components.csv", components
    )
    with pytest.raises(ValueError, match="saturation vapour pressure is required"):
        run_bet_from_neutral(copied)


def test_multiple_blocks_are_rejected(
    synthetic_bundle: Path, tmp_path: Path,
):
    copied = _copy_bundle(synthetic_bundle, tmp_path)
    injections = pd.read_csv(copied / "injections.csv", keep_default_na=False)
    injections.loc[injections["injection_id"] == "injection-002", "block_id"] = (
        "block-002"
    )
    _write_table_and_refresh_manifest(copied, "injections.csv", injections)
    with pytest.raises(ValueError, match="exactly one acquisition block"):
        run_bet_from_neutral(copied)


def test_probe_override_cannot_relabel_bundle(synthetic_bundle: Path):
    with pytest.raises(ValueError, match="does not uniquely match"):
        run_bet_from_neutral(synthetic_bundle, probe="different probe")


def test_small_fixture_is_numerical_but_not_reportable(synthetic_bundle: Path):
    result = run_bet_from_neutral(synthetic_bundle)
    assert result.qc.few_points
    assert result.classification is not None
    assert result.classification.bet_applicable is False


def test_closed_form_type_ii_bundle_recovers_known_bet_result(
    synthetic_bet_bundle: Path,
):
    result = run_bet_from_neutral(synthetic_bet_bundle)
    expected_ssa = 0.01 * 1e-3 * 6.02214076e23 * 6.3e-19
    assert result.n_monolayer_mmol_g == pytest.approx(0.01, rel=0.01)
    assert result.C_bet == pytest.approx(20.0, rel=0.02)
    assert result.ssa_m2_g == pytest.approx(expected_ssa, rel=0.01)
    assert result.r_squared > 0.9999
    assert result.n_points >= 8
    assert result.qc.passed
    assert result.classification.isotherm_type == "II"
    assert result.classification.bet_applicable is True


def test_synthetic_rising_retention_variant_is_rejected_as_type_iii(
    synthetic_bet_bundle: Path, tmp_path: Path,
):
    copied = _copy_bundle(synthetic_bet_bundle, tmp_path)
    traces = pd.read_csv(copied / "traces.csv", keep_default_na=False)
    injections = pd.read_csv(copied / "injections.csv", keep_default_na=False)
    probe_rows = injections[injections["role"] == "probe"]
    flow_mL_min = 1.85e-7 * 60.0 * 1e6
    for row in probe_rows.itertuples(index=False):
        mask = traces["injection_id"] == row.injection_id
        time_min = traces.loc[mask, "time_s"].to_numpy(dtype=float) / 60.0
        amplitude = float(traces.loc[mask, "signal_raw"].max())
        retention_volume_mL = 1.0 + 8.0 * float(row.target_coverage_fraction)
        center = 0.4 + retention_volume_mL / flow_mL_min
        signal = amplitude * np.exp(-0.5 * ((time_min - center) / 0.025) ** 2)
        traces.loc[mask, "signal_raw"] = signal
    _write_table_and_refresh_manifest(copied, "traces.csv", traces)

    result = run_bet_from_neutral(copied)
    assert result.qc.vn_increasing_with_concentration
    assert result.classification.isotherm_type in {"III", "II/III borderline"}
    assert result.classification.bet_applicable is False


def test_cofm_uses_matched_neutral_dead_time(synthetic_bundle: Path):
    result = run_bet_from_neutral(synthetic_bundle, retention_mode="cofm")
    first = min(result.injections, key=lambda item: item.injection_number)
    assert first.net_retention_time_min == pytest.approx(
        first.peak_cofm_time - result.diagnostics.methane.mean_cofm_min
    )
    assert first.net_retention_time_min != pytest.approx(
        first.peak_cofm_time - result.diagnostics.methane.mean_max_min
    )


def test_declared_false_does_not_suppress_clipping_detection(
    synthetic_bundle: Path, monkeypatch: pytest.MonkeyPatch,
):
    monkeypatch.setattr(bet_workflow, "_is_clipped", lambda signal: True)
    result = run_bet_from_neutral(synthetic_bundle)
    assert result.qc.peak_saturation
    assert all(item.peak_clipped for item in result.injections)


def test_loop_comparison_handles_degenerate_fit_without_crashing(
    synthetic_bundle: Path, tmp_path: Path,
):
    output = tmp_path / "loop-output"
    main([
        "--neutral-bundle", str(synthetic_bundle),
        "--concentration", "loop",
        "--output", str(output),
    ])
    record = json.loads((output / "bet_run.json").read_text(encoding="utf-8"))
    assert record["settings"]["concentration_mode"] == "loop"
    assert record["result"]["fit_range"] is None


def test_generated_readme_reports_actual_nondefault_conventions(
    synthetic_bundle: Path, tmp_path: Path,
):
    output = tmp_path / "nondefault-output"
    main([
        "--neutral-bundle", str(synthetic_bundle),
        "--retention", "cofm",
        "--origin", "rectangular",
        "--no-pressure-correction",
        "--output", str(output),
    ])
    readme = (output / "README.md").read_text(encoding="utf-8")
    assert "Retention: `cofm`" in readme
    assert "Isotherm origin: `rectangular`" in readme
    assert "Pressure correction: DISABLED" in readme
    assert "Pressure correction: applied" not in readme


def test_cli_error_is_concise_and_does_not_create_output(
    synthetic_bundle: Path, tmp_path: Path,
):
    output = tmp_path / "should-not-exist"
    with pytest.raises(SystemExit, match="does not uniquely match"):
        main([
            "--neutral-bundle", str(synthetic_bundle),
            "--probe", "not-declared",
            "--output", str(output),
        ])
    assert not output.exists()


def test_bet_cli_writes_auditable_outputs_without_input_path(
    synthetic_bundle: Path,
    tmp_path: Path,
):
    output = tmp_path / "bet-output"
    main([
        "--neutral-bundle", str(synthetic_bundle),
        "--output", str(output),
    ])
    expected = {
        "README.md",
        "bet_diagnostics.pdf",
        "bet_diagnostics.png",
        "bet_injections.csv",
        "bet_isotherm.csv",
        "bet_linearization.csv",
        "bet_run.json",
    }
    assert {path.name for path in output.iterdir()} == expected
    record_text = (output / "bet_run.json").read_text(encoding="utf-8")
    record = json.loads(record_text)
    assert record["settings"]["package_version"] == __version__
    assert record["settings"]["concentration_mode"] == "eluted"
    assert record["settings"]["pressure_correction"] is True
    assert record["input"]["dataset_id"] == "synthetic-peak-shape-001"
    assert len(record["input"]["manifest_digest"]) == 64
    assert record["reportability"]["bet_applicable"] is False
    assert record["declared_provenance"]["vapor_pressure_model_ids"] == [
        "synthetic-constant-v1"
    ]
    assert record["input"]["dead_time_injection_ids"] == [
        "injection-000", "injection-003"
    ]
    assert record["result"]["flow_source_channels"] == ["flow-channel-1"]
    assert record["result"]["pressure_source"] == "measured"
    assert record["result"]["sample_id"] == "sample-synthetic-001"
    assert str(synthetic_bundle) not in record_text


def test_cli_runs_packaged_closed_form_example(tmp_path: Path):
    output = tmp_path / "packaged-example"
    main(["--synthetic-example", "--output", str(output)])
    record = json.loads((output / "bet_run.json").read_text(encoding="utf-8"))
    assert record["reportability"]["bet_applicable"] is True
    assert record["result"]["ssa_m2_g"] == pytest.approx(3.794, rel=0.01)


def test_full_peak_declared_false_does_not_suppress_clipping_detection(
    synthetic_bundle: Path, tmp_path: Path,
):
    copied = _copy_bundle(synthetic_bundle, tmp_path)
    traces = pd.read_csv(copied / "traces.csv", keep_default_na=False)
    mask = traces["injection_id"] == "injection-001"
    peak_indices = traces.loc[mask, "signal_raw"].astype(float).nlargest(8).index
    traces.loc[peak_indices, "signal_raw"] = float(
        traces.loc[peak_indices, "signal_raw"].astype(float).max()
    )
    _write_table_and_refresh_manifest(copied, "traces.csv", traces)

    blocks = build_trace_dataset_from_neutral(
        {"synthetic": copied}, n_cells=40, verbose=False
    )
    first = next(
        injection for injection in blocks[0].injections
        if injection.name == "injection-001"
    )
    assert first.clipped is True
