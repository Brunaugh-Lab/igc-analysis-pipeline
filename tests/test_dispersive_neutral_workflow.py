"""End-to-end tests for the source-neutral dispersive workflow."""

from __future__ import annotations

import hashlib
import json
import shutil
import stat
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

import igc_analysis.analysis.dispersive_workflow as dispersive_workflow
from igc_analysis import __version__
from igc_analysis.analysis.dispersive_workflow import run_dispersive_from_neutral
from igc_analysis.analysis.peak_detection import find_peak_cofm
from igc_analysis.cli.dispersive import main
from igc_analysis.io.neutral_data import bundled_contract_path, read_neutral_bundle


@pytest.fixture()
def synthetic_dispersive_bundle() -> Path:
    return bundled_contract_path() / "examples" / "synthetic_dispersive_profile"


def _copy_bundle(source: Path, tmp_path: Path) -> Path:
    copied = tmp_path / "bundle"
    shutil.copytree(source, copied)
    return copied


def _write_table_and_refresh_manifest(
    bundle: Path, filename: str, table: pd.DataFrame
) -> None:
    path = bundle / filename
    table.to_csv(path, index=False)
    manifest_path = bundle / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["files"][filename]["row_count"] = len(table)
    manifest["files"][filename]["sha256"] = hashlib.sha256(
        path.read_bytes()
    ).hexdigest()
    manifest_path.write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )


def test_closed_form_bundle_recovers_known_dispersive_profile(
    synthetic_dispersive_bundle: Path,
):
    result = run_dispersive_from_neutral(synthetic_dispersive_bundle)

    assert set(result.injections["probe_id"]) == {
        "probe-homolog-08",
        "probe-homolog-09",
        "probe-homolog-10",
    }
    assert result.reportable is True
    assert result.qc["pass"] is True
    assert result.qc["profile_shape"] == "flat"
    assert len(result.gamma_d) == 4
    assert (result.gamma_d["n_alkanes"] == 3).all()
    assert result.gamma_d["gamma_d_mJm2"].tolist() == pytest.approx(
        [40.0] * 4, abs=0.20
    )
    assert result.gamma_d["r_squared"].min() > 0.9999
    assert (result.injections["james_martin_factor"] < 1.0).all()
    assert result.gamma_d["delta_cofm_pm"].abs().min() > 1.0
    assert result.gamma_d["W_cohesion_d_mJm2"].to_numpy() == pytest.approx(
        2.0 * result.gamma_d["gamma_d_mJm2"].to_numpy()
    )


def test_actual_coverage_and_matched_dead_time_are_used(
    synthetic_dispersive_bundle: Path,
):
    result = run_dispersive_from_neutral(synthetic_dispersive_bundle)
    first = result.injections.iloc[1]
    assert first["actual_coverage"] != pytest.approx(first["target_coverage"])
    assert first["net_retention_cofm_min"] == pytest.approx(
        first["peak_cofm_time_min"] - result.dead_time_cofm_min
    )
    assert first["net_retention_peak_max_min"] == pytest.approx(
        first["peak_max_time_min"] - result.dead_time_peak_max_min
    )
    assert set(result.interpolated["interpolation_status"]) == {"interpolated"}


def test_declared_provenance_is_preserved(
    synthetic_dispersive_bundle: Path,
):
    result = run_dispersive_from_neutral(synthetic_dispersive_bundle)
    assert result.surface_area_source == "synthetic-declared-ssa-v1"
    assert result.properties_sources == ("synthetic-reference-properties-v1",)
    assert result.calibration_sources == ("synthetic-calibration-v1",)
    assert result.flow_source_channels == ("flow-channel-synthetic",)
    assert result.pressure_basis == ("declared_absolute_inlet_outlet",)
    assert result.pressure_roles == ("measured",)


@pytest.mark.parametrize(
    ("column", "message"),
    [
        ("specific_surface_area_m2_g", "SSA and surface_area_source must occur together"),
        ("surface_area_source", "SSA and surface_area_source must occur together"),
    ],
)
def test_missing_surface_area_or_provenance_is_rejected(
    synthetic_dispersive_bundle: Path,
    tmp_path: Path,
    column: str,
    message: str,
):
    copied = _copy_bundle(synthetic_dispersive_bundle, tmp_path)
    experiment = pd.read_csv(copied / "experiment.csv", keep_default_na=False)
    experiment[column] = experiment[column].astype(object)
    experiment.loc[0, column] = ""
    _write_table_and_refresh_manifest(copied, "experiment.csv", experiment)
    with pytest.raises(ValueError, match=message):
        run_dispersive_from_neutral(copied)


def test_missing_probe_property_provenance_is_rejected(
    synthetic_dispersive_bundle: Path, tmp_path: Path,
):
    copied = _copy_bundle(synthetic_dispersive_bundle, tmp_path)
    properties = pd.read_csv(copied / "probe_properties.csv", keep_default_na=False)
    properties.loc[properties["carbon_number"] != "", "properties_source"] = ""
    _write_table_and_refresh_manifest(copied, "probe_properties.csv", properties)
    with pytest.raises(ValueError, match="properties_source: required value is empty"):
        run_dispersive_from_neutral(copied)


def test_missing_calibration_provenance_is_rejected(
    synthetic_dispersive_bundle: Path, tmp_path: Path,
):
    copied = _copy_bundle(synthetic_dispersive_bundle, tmp_path)
    calibration = pd.read_csv(copied / "calibration.csv", keep_default_na=False)
    calibration.loc[0, "calibration_source"] = ""
    _write_table_and_refresh_manifest(copied, "calibration.csv", calibration)
    with pytest.raises(ValueError, match="calibration_source: required value is empty"):
        run_dispersive_from_neutral(copied)


def test_co_injected_analytes_are_rejected(
    synthetic_dispersive_bundle: Path, tmp_path: Path,
):
    copied = _copy_bundle(synthetic_dispersive_bundle, tmp_path)
    components = pd.read_csv(
        copied / "injection_components.csv", keep_default_na=False
    )
    duplicate = components[
        components["injection_id"] == "injection-08-001"
    ].iloc[0].copy()
    duplicate["component_index"] = 1
    components = pd.concat(
        [components, pd.DataFrame([duplicate])], ignore_index=True
    ).sort_values(["injection_id", "component_index"])
    _write_table_and_refresh_manifest(
        copied, "injection_components.csv", components
    )
    with pytest.raises(ValueError, match="exactly one analyte component"):
        run_dispersive_from_neutral(copied)


def test_varying_detector_gain_is_recorded_and_flagged(
    synthetic_dispersive_bundle: Path, tmp_path: Path,
):
    copied = _copy_bundle(synthetic_dispersive_bundle, tmp_path)
    injections = pd.read_csv(copied / "injections.csv", keep_default_na=False)
    injections.loc[
        injections["injection_id"] == "injection-08-001", "detector_gain"
    ] = 2.0
    _write_table_and_refresh_manifest(copied, "injections.csv", injections)
    result = run_dispersive_from_neutral(copied)
    assert result.detector_gains == (1.0, 2.0)
    assert result.reportable is False
    assert any(
        flag["check"] == "detector_gain_variation"
        and flag["severity"] == "warning"
        for flag in result.qc["flags"]
    )


def test_explicit_homologous_probe_selection(
    synthetic_dispersive_bundle: Path,
):
    selected = [
        "probe-homolog-08",
        "probe-homolog-09",
        "probe-homolog-10",
    ]
    result = run_dispersive_from_neutral(
        synthetic_dispersive_bundle,
        homologous_probe_ids=selected,
    )
    assert set(result.injections["probe_id"]) == set(selected)

    with pytest.raises(ValueError, match="not analytes in this bundle"):
        run_dispersive_from_neutral(
            synthetic_dispersive_bundle,
            homologous_probe_ids=[*selected[:2], "probe-missing"],
        )


def test_mixed_carbon_numbered_bundle_requires_explicit_selection(
    synthetic_dispersive_bundle: Path, tmp_path: Path,
):
    copied = _copy_bundle(synthetic_dispersive_bundle, tmp_path)
    properties = pd.read_csv(copied / "probe_properties.csv", keep_default_na=False)
    extra_property = properties[
        properties["probe_id"] == "probe-homolog-08"
    ].iloc[0].copy()
    extra_property["probe_id"] = "probe-extra-carbon-07"
    extra_property["probe_name"] = "synthetic extra carbon analyte"
    extra_property["carbon_number"] = "7"
    properties = pd.concat(
        [properties, pd.DataFrame([extra_property])], ignore_index=True
    )
    _write_table_and_refresh_manifest(copied, "probe_properties.csv", properties)
    calibration = pd.read_csv(copied / "calibration.csv", keep_default_na=False)
    extra_calibration = calibration[
        calibration["probe_id"] == "probe-homolog-08"
    ].iloc[0].copy()
    extra_calibration["calibration_id"] = "calibration-probe-extra-carbon-07"
    extra_calibration["probe_id"] = "probe-extra-carbon-07"
    calibration = pd.concat(
        [calibration, pd.DataFrame([extra_calibration])], ignore_index=True
    )
    _write_table_and_refresh_manifest(copied, "calibration.csv", calibration)
    components = pd.read_csv(
        copied / "injection_components.csv", keep_default_na=False
    )
    components.loc[
        components["injection_id"] == "injection-08-001",
        ["probe_id", "calibration_id"],
    ] = ["probe-extra-carbon-07", "calibration-probe-extra-carbon-07"]
    _write_table_and_refresh_manifest(
        copied, "injection_components.csv", components
    )

    with pytest.raises(ValueError, match="require explicit homologous_probe_ids"):
        run_dispersive_from_neutral(copied)

    result = run_dispersive_from_neutral(
        copied,
        homologous_probe_ids=[
            "probe-homolog-08",
            "probe-homolog-09",
            "probe-homolog-10",
        ],
    )
    assert set(result.injections["probe_id"]) == {
        "probe-homolog-08",
        "probe-homolog-09",
        "probe-homolog-10",
    }


def test_missing_cross_section_has_a_located_error(
    synthetic_dispersive_bundle: Path, tmp_path: Path,
):
    copied = _copy_bundle(synthetic_dispersive_bundle, tmp_path)
    properties = pd.read_csv(copied / "probe_properties.csv", keep_default_na=False)
    properties["cross_section_m2"] = properties["cross_section_m2"].astype(object)
    properties.loc[
        properties["probe_id"] == "probe-homolog-10", "cross_section_m2"
    ] = ""
    _write_table_and_refresh_manifest(copied, "probe_properties.csv", properties)
    with pytest.raises(ValueError, match="cross_section_m2 is required"):
        run_dispersive_from_neutral(copied)


def test_unstable_measured_flow_is_rejected(
    synthetic_dispersive_bundle: Path, tmp_path: Path,
):
    copied = _copy_bundle(synthetic_dispersive_bundle, tmp_path)
    conditions = pd.read_csv(copied / "conditions.csv", keep_default_na=False)
    mask = (
        (conditions["injection_id"] == "injection-08-001")
        & (conditions["quantity"] == "flow_column")
    )
    conditions.loc[mask, "value"] = 2.5e-7
    _write_table_and_refresh_manifest(copied, "conditions.csv", conditions)
    with pytest.raises(ValueError, match="column flow is not stable"):
        run_dispersive_from_neutral(copied)


def test_target_only_conditions_are_rejected(
    synthetic_dispersive_bundle: Path, tmp_path: Path,
):
    copied = _copy_bundle(synthetic_dispersive_bundle, tmp_path)
    conditions = pd.read_csv(copied / "conditions.csv", keep_default_na=False)
    mask = (
        (conditions["injection_id"] == "injection-08-001")
        & conditions["quantity"].isin(["column_temperature", "flow_column"])
    )
    conditions.loc[mask, "value_role"] = "target"
    conditions.loc[mask, "measurement_basis"] = "declared_target"
    _write_table_and_refresh_manifest(copied, "conditions.csv", conditions)
    with pytest.raises(ValueError, match="requires measured column temperature"):
        run_dispersive_from_neutral(copied)


def test_pressure_correction_requires_measured_inlet_or_drop(
    synthetic_dispersive_bundle: Path, tmp_path: Path,
):
    copied = _copy_bundle(synthetic_dispersive_bundle, tmp_path)
    conditions = pd.read_csv(copied / "conditions.csv", keep_default_na=False)
    mask = conditions["quantity"].isin(["pressure_inlet", "pressure_drop"])
    conditions.loc[mask, "value_role"] = "target"
    conditions.loc[mask, "measurement_basis"] = "declared_target"
    _write_table_and_refresh_manifest(copied, "conditions.csv", conditions)
    with pytest.raises(
        ValueError, match="requires measured pressure_inlet or pressure_drop"
    ):
        run_dispersive_from_neutral(copied)


def test_measured_drop_takes_precedence_over_target_absolute_pressures(
    synthetic_dispersive_bundle: Path, tmp_path: Path,
):
    copied = _copy_bundle(synthetic_dispersive_bundle, tmp_path)
    conditions = pd.read_csv(copied / "conditions.csv", keep_default_na=False)
    mask = conditions["quantity"].isin(["pressure_inlet", "pressure_outlet"])
    conditions.loc[mask, "value_role"] = "target"
    conditions.loc[mask, "measurement_basis"] = "declared_target"
    _write_table_and_refresh_manifest(copied, "conditions.csv", conditions)
    result = run_dispersive_from_neutral(copied)
    assert result.pressure_roles == ("measured",)
    assert result.pressure_basis == (
        "declared_drop_plus_ambient_absolute_outlet",
    )


def test_dead_time_conditions_participate_in_stability_gate(
    synthetic_dispersive_bundle: Path, tmp_path: Path,
):
    copied = _copy_bundle(synthetic_dispersive_bundle, tmp_path)
    conditions = pd.read_csv(copied / "conditions.csv", keep_default_na=False)
    mask = (
        (conditions["injection_id"] == "dead-time-000")
        & (conditions["quantity"] == "flow_column")
    )
    conditions.loc[mask, "value"] = 2.5e-7
    _write_table_and_refresh_manifest(copied, "conditions.csv", conditions)
    with pytest.raises(ValueError, match="column flow is not stable"):
        run_dispersive_from_neutral(copied)


def test_fractional_carbon_number_is_rejected(
    synthetic_dispersive_bundle: Path, tmp_path: Path,
):
    copied = _copy_bundle(synthetic_dispersive_bundle, tmp_path)
    properties = pd.read_csv(copied / "probe_properties.csv", keep_default_na=False)
    properties.loc[
        properties["probe_id"] == "probe-homolog-10", "carbon_number"
    ] = "10.5"
    _write_table_and_refresh_manifest(copied, "probe_properties.csv", properties)
    with pytest.raises(ValueError, match="carbon_number: '10.5' is not a valid integer"):
        run_dispersive_from_neutral(copied)


def test_sparse_probe_coverage_is_rejected(
    synthetic_dispersive_bundle: Path, tmp_path: Path,
):
    copied = _copy_bundle(synthetic_dispersive_bundle, tmp_path)
    injections = pd.read_csv(copied / "injections.csv", keep_default_na=False)
    mask = injections["injection_id"].isin(["injection-10-011", "injection-10-012"])
    injections.loc[mask, "role"] = "reference"
    _write_table_and_refresh_manifest(copied, "injections.csv", injections)
    with pytest.raises(ValueError, match="at least three calibrated coverage points"):
        run_dispersive_from_neutral(copied)


def test_clipping_is_critical_even_when_declared_false(
    synthetic_dispersive_bundle: Path, monkeypatch: pytest.MonkeyPatch,
):
    monkeypatch.setattr(dispersive_workflow, "_is_clipped", lambda signal: True)
    result = run_dispersive_from_neutral(synthetic_dispersive_bundle)
    assert result.reportable is False
    assert result.qc["pass"] is False
    assert any(
        flag["check"] == "detector_clipping"
        and flag["severity"] == "critical"
        for flag in result.qc["flags"]
    )


def test_center_of_mass_weights_nonuniform_time_spacing():
    time = pd.Series([0.0, 1.0, 3.0]).to_numpy()
    signal = pd.Series([1.0, 1.0, 0.0]).to_numpy()
    assert find_peak_cofm(time, signal) == pytest.approx(0.75)


def test_disabling_extrapolation_fails_closed_at_incomplete_edges(
    synthetic_dispersive_bundle: Path, tmp_path: Path,
):
    copied = _copy_bundle(synthetic_dispersive_bundle, tmp_path)
    injections = pd.read_csv(copied / "injections.csv", keep_default_na=False)
    injections.loc[
        injections["injection_id"] == "injection-08-001",
        "target_coverage_fraction",
    ] = "0.005"
    _write_table_and_refresh_manifest(copied, "injections.csv", injections)

    historical_mode = run_dispersive_from_neutral(copied, extrapolate=True)
    assert historical_mode.reportable is False
    assert any(
        flag["check"] == "coverage_extrapolation"
        for flag in historical_mode.qc["flags"]
    )

    fail_closed = run_dispersive_from_neutral(copied, extrapolate=False)
    assert fail_closed.reportable is False
    assert (fail_closed.gamma_d["n_alkanes"] < 3).any()
    assert "outside_measured_range" in set(
        fail_closed.interpolated["interpolation_status"]
    )
    assert any(
        flag["check"] == "coverage_outside_measured_range"
        for flag in fail_closed.qc["flags"]
    )


def test_cli_writes_auditable_outputs_without_local_path(
    synthetic_dispersive_bundle: Path, tmp_path: Path,
):
    output = tmp_path / "dispersive-output"
    main([
        "--neutral-bundle", str(synthetic_dispersive_bundle),
        "--output", str(output),
    ])
    assert {path.name for path in output.iterdir()} == {
        "README.md",
        "dispersive_alkane_lines.pdf",
        "dispersive_alkane_lines.png",
        "dispersive_injections.csv",
        "dispersive_interpolated_vn.csv",
        "dispersive_profile.csv",
        "dispersive_profile.pdf",
        "dispersive_profile.png",
        "dispersive_run.json",
    }
    record_text = (output / "dispersive_run.json").read_text(encoding="utf-8")
    record = json.loads(record_text)
    assert record["settings"]["package_version"] == __version__
    assert record["settings"]["primary_retention_mode"] == "cofm"
    assert record["settings"]["homologous_probe_selection"] == "all_carbon_numbered"
    assert record["result"]["reportable"] is True
    assert record["input"]["surface_area_source"] == "synthetic-declared-ssa-v1"
    assert record["input"]["pressure_roles"] == ["measured"]
    assert str(synthetic_dispersive_bundle) not in record_text
    assert stat.S_IMODE(output.stat().st_mode) == 0o700
    with pytest.raises(SystemExit, match="output already exists"):
        main([
            "--neutral-bundle", str(synthetic_dispersive_bundle),
            "--output", str(output),
        ])


def test_cli_records_explicit_homologous_probe_selection(
    synthetic_dispersive_bundle: Path, tmp_path: Path,
):
    output = tmp_path / "explicit-selection-output"
    main([
        "--neutral-bundle", str(synthetic_dispersive_bundle),
        "--homologous-probe-id", "probe-homolog-08",
        "--homologous-probe-id", "probe-homolog-09",
        "--homologous-probe-id", "probe-homolog-10",
        "--output", str(output),
    ])
    record = json.loads(
        (output / "dispersive_run.json").read_text(encoding="utf-8")
    )
    assert record["settings"]["homologous_probe_selection"] == "explicit"


def test_cli_error_does_not_create_output(
    synthetic_dispersive_bundle: Path, tmp_path: Path,
):
    copied = _copy_bundle(synthetic_dispersive_bundle, tmp_path)
    experiment = pd.read_csv(copied / "experiment.csv", keep_default_na=False)
    experiment.loc[0, "surface_area_source"] = ""
    _write_table_and_refresh_manifest(copied, "experiment.csv", experiment)
    output = tmp_path / "should-not-exist"
    with pytest.raises(SystemExit, match="SSA and surface_area_source must occur together"):
        main(["--neutral-bundle", str(copied), "--output", str(output)])
    assert not output.exists()


def test_cli_runs_packaged_closed_form_example(tmp_path: Path):
    output = tmp_path / "packaged-example"
    main(["--synthetic-example", "--output", str(output)])
    record = json.loads(
        (output / "dispersive_run.json").read_text(encoding="utf-8")
    )
    assert record["result"]["reportable"] is True
    assert record["result"]["gamma_d_min_mJm2"] == pytest.approx(40.0, abs=0.20)


def test_cli_handles_coverage_with_no_surviving_fit(
    synthetic_dispersive_bundle: Path, tmp_path: Path,
):
    copied = _copy_bundle(synthetic_dispersive_bundle, tmp_path)
    injections = pd.read_csv(copied / "injections.csv", keep_default_na=False)
    injections.loc[
        injections["injection_id"] == "injection-08-001",
        "target_coverage_fraction",
    ] = "0.001"
    _write_table_and_refresh_manifest(copied, "injections.csv", injections)
    output = tmp_path / "no-fit-output"
    main([
        "--neutral-bundle", str(copied),
        "--no-extrapolation",
        "--output", str(output),
    ])
    record = json.loads(
        (output / "dispersive_run.json").read_text(encoding="utf-8")
    )
    profile = pd.read_csv(output / "dispersive_profile.csv")
    assert record["result"]["reportable"] is False
    assert (profile["n_alkanes"] == 0).any()


def test_committed_fixture_is_reproducible(
    synthetic_dispersive_bundle: Path, tmp_path: Path,
):
    generated = tmp_path / "regenerated"
    subprocess.run(
        [
            sys.executable,
            "scripts/generate_synthetic_dispersive_fixture.py",
            "--output",
            str(generated),
        ],
        check=True,
    )
    generated_bundle = read_neutral_bundle(generated)
    committed_bundle = read_neutral_bundle(synthetic_dispersive_bundle)
    assert generated_bundle.dataset_id == committed_bundle.dataset_id
    assert generated_bundle.tables.keys() == committed_bundle.tables.keys()
    for filename in generated_bundle.tables:
        pd.testing.assert_frame_equal(
            generated_bundle.table(filename),
            committed_bundle.table(filename),
            check_dtype=False,
            check_exact=False,
            rtol=1e-11,
            atol=1e-14,
        )
