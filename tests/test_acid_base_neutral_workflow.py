"""End-to-end tests for the source-neutral acid/base workflow."""

from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import igc_analysis.analysis.acid_base as legacy_acid_base
from igc_analysis.analysis.acid_base_workflow import run_acid_base_from_neutral
from igc_analysis.cli.acid_base import main
from igc_analysis.io.neutral_data import bundled_contract_path


HOMOLOGS = ("probe-homolog-08", "probe-homolog-09", "probe-homolog-10")
POLAR = ("probe-polar-01", "probe-polar-02", "probe-polar-03")


@pytest.fixture()
def synthetic_bundle() -> Path:
    return bundled_contract_path() / "examples" / "synthetic_dispersive_profile"


def _refresh(bundle: Path, filename: str, table: pd.DataFrame) -> None:
    path = bundle / filename
    table.to_csv(path, index=False)
    manifest_path = bundle / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["files"][filename]["row_count"] = len(table)
    manifest["files"][filename]["sha256"] = hashlib.sha256(
        path.read_bytes()
    ).hexdigest()
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


def test_closed_form_profile_recovers_known_ka_kb(synthetic_bundle: Path):
    result = run_acid_base_from_neutral(
        synthetic_bundle, homologous_probe_ids=HOMOLOGS, polar_probe_ids=POLAR
    )
    primary = result.profile[result.profile["retention_mode"] == "cofm"]
    assert result.reportable is True
    assert result.qc["pass"] is True
    assert len(primary) == 4
    assert (primary["fit_method"] == "regression").all()
    assert (primary["n_homologs"] == 3).all()
    assert (primary["n_polar_probes"] == 3).all()
    assert primary["Ka"].tolist() == pytest.approx([0.03] * 4, abs=0.001)
    assert primary["Kb"].tolist() == pytest.approx([0.05] * 4, abs=0.001)
    assert primary["r_squared"].min() > 0.999


def test_peak_max_is_retained_only_as_sensitivity(synthetic_bundle: Path):
    result = run_acid_base_from_neutral(
        synthetic_bundle, homologous_probe_ids=HOMOLOGS, polar_probe_ids=POLAR
    )
    assert set(result.profile["retention_mode"]) == {"cofm", "peak_max"}
    cofm = result.profile[result.profile["retention_mode"] == "cofm"]["Ka"].to_numpy()
    peak = result.profile[result.profile["retention_mode"] == "peak_max"]["Ka"].to_numpy()
    assert abs(cofm - peak).max() > 0.0005


def test_workflow_does_not_use_name_based_legacy_properties(
    synthetic_bundle: Path, monkeypatch: pytest.MonkeyPatch,
):
    def forbidden_lookup(*args, **kwargs):
        raise AssertionError("legacy property lookup was called")

    monkeypatch.setattr(legacy_acid_base, "get_probe", forbidden_lookup)
    result = run_acid_base_from_neutral(
        synthetic_bundle, homologous_probe_ids=HOMOLOGS, polar_probe_ids=POLAR
    )
    assert result.reportable is True


def test_three_explicit_polar_probes_are_required(synthetic_bundle: Path):
    with pytest.raises(ValueError, match="at least three unique polar"):
        run_acid_base_from_neutral(
            synthetic_bundle,
            homologous_probe_ids=HOMOLOGS,
            polar_probe_ids=POLAR[:2],
        )


def test_missing_polar_property_or_source_is_rejected(
    synthetic_bundle: Path, tmp_path: Path,
):
    copied = tmp_path / "bundle"
    shutil.copytree(synthetic_bundle, copied)
    properties = pd.read_csv(copied / "probe_properties.csv", keep_default_na=False)
    properties.loc[
        properties["probe_id"] == POLAR[0], "donor_number_kJ_mol"
    ] = ""
    _refresh(copied, "probe_properties.csv", properties)
    with pytest.raises(ValueError, match="donor_number_kJ_mol is required"):
        run_acid_base_from_neutral(
            copied, homologous_probe_ids=HOMOLOGS, polar_probe_ids=POLAR
        )


def test_mislabeled_j_per_m2_liquid_tension_is_rejected(
    synthetic_bundle: Path, tmp_path: Path,
):
    copied = tmp_path / "bundle"
    shutil.copytree(synthetic_bundle, copied)
    properties = pd.read_csv(copied / "probe_properties.csv", keep_default_na=False)
    properties.loc[
        properties["probe_id"] == POLAR[0], "gamma_l_d_mJ_m2"
    ] = "0.018"
    _refresh(copied, "probe_properties.csv", properties)
    with pytest.raises(ValueError, match="unit/plausibility gate"):
        run_acid_base_from_neutral(
            copied, homologous_probe_ids=HOMOLOGS, polar_probe_ids=POLAR
        )


def test_nonpositive_polar_retention_is_retained_and_nonreportable(
    synthetic_bundle: Path, tmp_path: Path,
):
    copied = tmp_path / "bundle"
    shutil.copytree(synthetic_bundle, copied)
    traces = pd.read_csv(copied / "traces.csv", keep_default_na=False)
    injections = pd.read_csv(copied / "injections.csv", keep_default_na=False)
    affected = injections.loc[
        injections["injection_id"].str.startswith("injection-polar-"),
        "injection_id",
    ].iloc[:4]
    for injection_id in affected:
        mask = traces["injection_id"] == injection_id
        time_s = traces.loc[mask, "time_s"].to_numpy(dtype=float)
        traces.loc[mask, "signal_raw"] = np.exp(
            -0.5 * ((time_s - 6.0) / 1.0) ** 2
        )
    _refresh(copied, "traces.csv", traces)
    result = run_acid_base_from_neutral(
        copied, homologous_probe_ids=HOMOLOGS, polar_probe_ids=POLAR
    )
    assert result.reportable is False
    assert any(
        flag["check"] == "positive_net_retention"
        and flag["severity"] == "critical"
        for flag in result.qc["flags"]
    )


def test_all_nonpositive_polar_retention_returns_structured_failure(
    synthetic_bundle: Path, tmp_path: Path,
):
    copied = tmp_path / "bundle"
    shutil.copytree(synthetic_bundle, copied)
    traces = pd.read_csv(copied / "traces.csv", keep_default_na=False)
    injections = pd.read_csv(copied / "injections.csv", keep_default_na=False)
    affected = injections.loc[
        injections["injection_id"].str.startswith("injection-polar-"),
        "injection_id",
    ]
    for injection_id in affected:
        mask = traces["injection_id"] == injection_id
        time_s = traces.loc[mask, "time_s"].to_numpy(dtype=float)
        traces.loc[mask, "signal_raw"] = np.exp(
            -0.5 * ((time_s - 6.0) / 1.0) ** 2
        )
    _refresh(copied, "traces.csv", traces)
    result = run_acid_base_from_neutral(
        copied, homologous_probe_ids=HOMOLOGS, polar_probe_ids=POLAR
    )
    assert result.delta_g_sp.empty
    assert tuple(result.delta_g_sp.columns)
    assert result.reportable is False
    assert result.qc["pass"] is False


def test_coverage_extrapolation_is_nonreportable(
    synthetic_bundle: Path, tmp_path: Path,
):
    copied = tmp_path / "bundle"
    shutil.copytree(synthetic_bundle, copied)
    injections = pd.read_csv(copied / "injections.csv", keep_default_na=False)
    injections.loc[
        injections["target_coverage_fraction"] == "0.04",
        "target_coverage_fraction",
    ] = "0.05"
    _refresh(copied, "injections.csv", injections)
    result = run_acid_base_from_neutral(
        copied, homologous_probe_ids=HOMOLOGS, polar_probe_ids=POLAR
    )
    assert result.reportable is False
    assert any(flag["check"] == "coverage_mapping" for flag in result.qc["flags"])


def test_detector_gain_variation_is_nonreportable(
    synthetic_bundle: Path, tmp_path: Path,
):
    copied = tmp_path / "bundle"
    shutil.copytree(synthetic_bundle, copied)
    injections = pd.read_csv(copied / "injections.csv", keep_default_na=False)
    injections.loc[injections["role"] == "probe", "detector_gain"] = 2.0
    injections.loc[injections.index[-2], "detector_gain"] = 1.0
    _refresh(copied, "injections.csv", injections)
    result = run_acid_base_from_neutral(
        copied, homologous_probe_ids=HOMOLOGS, polar_probe_ids=POLAR
    )
    assert result.reportable is False
    assert any(
        flag["check"] == "detector_gain_variation"
        for flag in result.qc["flags"]
    )


def test_negative_delta_g_sp_is_visible_as_warning(
    synthetic_bundle: Path, tmp_path: Path,
):
    copied = tmp_path / "bundle"
    shutil.copytree(synthetic_bundle, copied)
    properties = pd.read_csv(copied / "probe_properties.csv", keep_default_na=False)
    properties.loc[
        properties["probe_id"] == POLAR[0], "gamma_l_d_mJ_m2"
    ] = "100.0"
    _refresh(copied, "probe_properties.csv", properties)
    result = run_acid_base_from_neutral(
        copied, homologous_probe_ids=HOMOLOGS, polar_probe_ids=POLAR
    )
    assert any(
        flag["check"] == "negative_delta_g_sp"
        for flag in result.qc["flags"]
    )


def test_poor_schultz_reference_line_is_not_reportable(
    synthetic_bundle: Path, tmp_path: Path,
):
    copied = tmp_path / "bundle"
    shutil.copytree(synthetic_bundle, copied)
    properties = pd.read_csv(copied / "probe_properties.csv", keep_default_na=False)
    properties.loc[
        properties["probe_id"] == HOMOLOGS[1], "gamma_l_d_mJ_m2"
    ] = "100.0"
    _refresh(copied, "probe_properties.csv", properties)
    result = run_acid_base_from_neutral(
        copied, homologous_probe_ids=HOMOLOGS, polar_probe_ids=POLAR
    )
    assert result.reportable is False
    assert any(
        flag["check"] in {"schultz_r_squared", "schultz_gamma_d_bounds"}
        and flag["severity"] == "critical"
        for flag in result.qc["flags"]
    )


def test_cli_writes_auditable_atomic_outputs(
    synthetic_bundle: Path, tmp_path: Path,
):
    output = tmp_path / "acid-base-output"
    main(["--synthetic-example", "--output", str(output)])
    expected = {
        "README.md", "acid_base_injections.csv", "acid_base_interpolated_vn.csv",
        "acid_base_schultz_lines.csv", "acid_base_delta_g_sp.csv",
        "acid_base_profile.csv", "acid_base_profile.pdf", "acid_base_profile.png",
        "acid_base_gutmann_fits.pdf", "acid_base_gutmann_fits.png",
        "acid_base_run.json",
    }
    assert {path.name for path in output.iterdir()} == expected
    record_text = (output / "acid_base_run.json").read_text(encoding="utf-8")
    record = json.loads(record_text)
    assert record["result"]["reportable"] is True
    assert record["settings"]["polar_probe_policy"] == "explicit_inclusion"
    assert record["settings"]["van_oss"].startswith("not_available")
    assert str(synthetic_bundle) not in record_text


def test_cli_refuses_existing_output(tmp_path: Path):
    output = tmp_path / "existing"
    output.mkdir()
    with pytest.raises(SystemExit, match="output already exists"):
        main(["--synthetic-example", "--output", str(output)])
