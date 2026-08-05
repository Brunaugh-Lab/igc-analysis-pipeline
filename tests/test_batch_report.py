"""Tests for explicit source-neutral batch orchestration."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import matplotlib
import pandas as pd
import pytest

from igc_analysis.cli import bet, report
from igc_analysis.io.neutral_data import bundled_contract_path


def _write_manifest(path: Path, data: dict) -> Path:
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    return path


def _one_bundle_manifest(bundle_path: str, *, settings: dict | None = None) -> dict:
    job = {
        "job_id": "bet-001",
        "analysis": "bet",
        "bundle_ids": ["bundle-001"],
    }
    if settings is not None:
        job["settings"] = settings
    return {
        "schema_version": report.BATCH_SCHEMA_VERSION,
        "batch_id": "batch-001",
        "bundles": [{"bundle_id": "bundle-001", "path": bundle_path}],
        "jobs": [job],
    }


def test_packaged_schema_matches_supported_analyses_and_settings():
    schema_path = report.bundled_batch_contract_path() / "schema.json"
    schema = json.loads(schema_path.read_text(encoding="utf-8"))
    job = schema["properties"]["jobs"]["items"]

    assert set(job["properties"]["analysis"]["enum"]) == report.SUPPORTED_ANALYSES
    assert job["properties"]["bundle_ids"]["maxItems"] == 1
    for analysis in report.SUPPORTED_ANALYSES:
        settings_schema = schema["$defs"][f"{analysis}_settings"]
        if "$ref" in settings_schema:
            settings_schema = schema["$defs"][
                settings_schema["$ref"].rsplit("/", 1)[-1]
            ]
        assert set(settings_schema["properties"]) == report.SETTING_KEYS[analysis]


def test_synthetic_batch_runs_every_neutral_workflow_transactionally(tmp_path):
    output = tmp_path / "results"
    font_settings = {
        key: matplotlib.rcParams[key] for key in ("pdf.fonttype", "ps.fonttype")
    }

    report.main(["--synthetic-example", "--output", str(output)])

    assert {
        key: matplotlib.rcParams[key] for key in font_settings
    } == font_settings
    summary = pd.read_csv(output / "batch_summary.csv", keep_default_na=False)
    assert list(summary.columns) == [
        "job_id", "analysis", "bundle_ids", "dataset_ids", "result_directory",
        "record_file", "reportability_scope", "reportable", "qc_status",
        "qc_summary", "selected_model",
    ]
    assert list(summary["analysis"]) == [
        "bet", "dispersive", "acid_base", "full_peak",
    ]
    assert list(summary["reportability_scope"]) == [
        "bet_ssa", "dispersive_profile", "acid_base_profile",
        "full_peak_recovered_ssa",
    ]
    assert summary.loc[0:2, "reportable"].tolist() == [True, True, True]
    assert summary.loc[3, "selected_model"] == "henry"
    assert summary.loc[3, "reportable"] in (False, 0)
    assert summary.loc[3, "qc_status"] == "NOT_COMBINED"
    assert "minimum mass balance" in summary.loc[3, "qc_summary"]
    record = json.loads((output / "batch_run.json").read_text(encoding="utf-8"))
    assert record["settings"]["batch_schema_version"] == report.BATCH_SCHEMA_VERSION
    assert record["input"]["batch_id"] == "batch-synthetic-001"
    assert len(record["jobs"]) == 4
    assert all("path" not in bundle for bundle in record["input"]["bundles"])
    assert str(tmp_path) not in json.dumps(record)
    for job_id in ("bet-001", "dispersive-001", "acid-base-001", "full-peak-001"):
        assert (output / job_id / "README.md").is_file()


def test_every_full_peak_job_uses_and_contains_its_pdf_font_settings(
    tmp_path, monkeypatch,
):
    from igc_analysis.cli import full_peak

    observed: list[tuple[int, int]] = []

    def fake_main(argv):
        observed.append((
            matplotlib.rcParams["pdf.fonttype"],
            matplotlib.rcParams["ps.fonttype"],
        ))
        output = Path(argv[argv.index("--output") + 1])
        output.mkdir()
        (output / "full_peak_run.json").write_text(json.dumps({
            "selected_model": "henry",
            "qc": {
                "mass_balance_min": 0.9,
                "all_params_identifiable": True,
                "cooperative": False,
            },
            "ssa": {"reportable": False},
        }), encoding="utf-8")

    monkeypatch.setattr(full_peak, "main", fake_main)
    fixture = bundled_contract_path() / "examples" / "synthetic_peak_shape"
    job = {
        "analysis": "full_peak", "bundle_ids": ["bundle"], "settings": {},
    }
    with matplotlib.rc_context({"pdf.fonttype": 3, "ps.fonttype": 3}):
        report._run_job(job, {"bundle": fixture}, tmp_path / "first")
        report._run_job(job, {"bundle": fixture}, tmp_path / "second")
        assert matplotlib.rcParams["pdf.fonttype"] == 3
        assert matplotlib.rcParams["ps.fonttype"] == 3

    assert observed == [(42, 42), (42, 42)]


def test_relative_bundle_path_is_resolved_from_manifest_directory(tmp_path):
    fixture = bundled_contract_path() / "examples" / "synthetic_bet_isotherm"
    shutil.copytree(fixture, tmp_path / "neutral" / "bundle")
    manifest = _write_manifest(
        tmp_path / "batch.json",
        _one_bundle_manifest("neutral/bundle"),
    )

    report.main(["--manifest", str(manifest), "--output", str(tmp_path / "out")])

    record = json.loads((tmp_path / "out" / "batch_run.json").read_text())
    assert record["input"]["bundles"][0]["dataset_id"] == "synthetic-bet-isotherm-001"


def test_batch_orchestration_preserves_child_scientific_tables(tmp_path):
    fixture = bundled_contract_path() / "examples" / "synthetic_bet_isotherm"
    manifest = _write_manifest(
        tmp_path / "batch.json", _one_bundle_manifest(str(fixture))
    )
    batch_output = tmp_path / "batch"
    standalone_output = tmp_path / "standalone"

    report.main(["--manifest", str(manifest), "--output", str(batch_output)])
    bet.main(["--neutral-bundle", str(fixture), "--output", str(standalone_output)])

    for name in ("bet_injections.csv", "bet_isotherm.csv", "bet_linearization.csv"):
        assert (batch_output / "bet-001" / name).read_bytes() == (
            standalone_output / name
        ).read_bytes()


@pytest.mark.parametrize(
    "mutate, message",
    [
        (lambda data: data.update(extra=True), "unknown keys"),
        (lambda data: data.update(schema_version="0.0.0"), "unsupported batch schema"),
        (
            lambda data: data["jobs"][0].update(bundle_ids=["missing"]),
            "unknown bundles",
        ),
        (
            lambda data: data["jobs"][0].update(analysis="humidity"),
            "unsupported analysis",
        ),
        (
            lambda data: data["jobs"][0].update(settings={"mystery": 1}),
            "unknown keys",
        ),
    ],
)
def test_manifest_ambiguity_fails_without_output(tmp_path, mutate, message):
    fixture = bundled_contract_path() / "examples" / "synthetic_bet_isotherm"
    data = _one_bundle_manifest(str(fixture))
    mutate(data)
    manifest = _write_manifest(tmp_path / "batch.json", data)
    output = tmp_path / "out"

    with pytest.raises(SystemExit, match=message):
        report.main(["--manifest", str(manifest), "--output", str(output)])

    assert not output.exists()


def test_duplicate_dataset_aliases_are_rejected(tmp_path):
    fixture = bundled_contract_path() / "examples" / "synthetic_bet_isotherm"
    data = _one_bundle_manifest(str(fixture))
    data["bundles"].append({"bundle_id": "bundle-002", "path": str(fixture)})
    manifest = _write_manifest(tmp_path / "batch.json", data)

    with pytest.raises(SystemExit, match="declared by more than one bundle"):
        report.main(["--manifest", str(manifest), "--output", str(tmp_path / "out")])


def test_job_failure_removes_all_staged_outputs(tmp_path):
    fixture = bundled_contract_path() / "examples" / "synthetic_bet_isotherm"
    data = _one_bundle_manifest(str(fixture), settings={"p0_min": "not-a-number"})
    manifest = _write_manifest(tmp_path / "batch.json", data)
    output = tmp_path / "out"

    with pytest.raises(SystemExit, match="job 'bet-001'.*p0_min"):
        report.main(["--manifest", str(manifest), "--output", str(output)])

    assert not output.exists()


def test_child_partial_output_is_removed_transactionally(tmp_path, monkeypatch):
    fixture = bundled_contract_path() / "examples" / "synthetic_bet_isotherm"
    manifest = _write_manifest(
        tmp_path / "batch.json", _one_bundle_manifest(str(fixture))
    )
    output = tmp_path / "out"

    def fail_after_partial_output(argv):
        staged_job = Path(argv[argv.index("--output") + 1])
        staged_job.mkdir()
        (staged_job / "partial.txt").write_text("partial", encoding="utf-8")
        raise ValueError("synthetic child failure")

    monkeypatch.setattr(report.bet, "main", fail_after_partial_output)

    with pytest.raises(SystemExit, match="synthetic child failure"):
        report.main(["--manifest", str(manifest), "--output", str(output)])

    assert not output.exists()
    assert not list(tmp_path.glob(".igc-report-*"))


def test_later_job_failure_removes_earlier_completed_job(tmp_path):
    fixture = bundled_contract_path() / "examples" / "synthetic_bet_isotherm"
    data = _one_bundle_manifest(str(fixture))
    data["jobs"].append({
        "job_id": "bet-invalid",
        "analysis": "bet",
        "bundle_ids": ["bundle-001"],
        "settings": {"p0_min": 0.8, "p0_max": 0.2},
    })
    manifest = _write_manifest(tmp_path / "batch.json", data)
    output = tmp_path / "out"

    with pytest.raises(SystemExit, match="job 'bet-invalid'.*BET bounds"):
        report.main(["--manifest", str(manifest), "--output", str(output)])

    assert not output.exists()


def test_cross_bundle_full_peak_fit_is_outside_batch_scope(tmp_path):
    fixture = bundled_contract_path() / "examples" / "synthetic_peak_shape"
    second = bundled_contract_path() / "examples" / "synthetic_dispersive_profile"
    data = {
        "schema_version": report.BATCH_SCHEMA_VERSION,
        "batch_id": "batch-001",
        "bundles": [
            {"bundle_id": "first", "path": str(fixture)},
            {"bundle_id": "second", "path": str(second)},
        ],
        "jobs": [{
            "job_id": "full-peak-001", "analysis": "full_peak",
            "bundle_ids": ["first", "second"],
        }],
    }
    manifest = _write_manifest(tmp_path / "batch.json", data)

    with pytest.raises(SystemExit, match="requires exactly one bundle"):
        report.main(["--manifest", str(manifest), "--output", str(tmp_path / "out")])


@pytest.mark.parametrize(
    "analysis, settings, message",
    [
        ("bet", {"p0_min": 0.8, "p0_max": 0.2}, "BET bounds"),
        ("bet", {"ambient_pressure_pa": 0}, "must be positive"),
        ("dispersive", {"max_temperature_span_k": -1}, "cannot be negative"),
        ("full_peak", {"n_cells": 1}, "must be at least 2"),
        ("full_peak", {"n_starts": 0}, "must be at least 1"),
        ("full_peak", {"cross_section_m2": 0}, "must be positive"),
        ("full_peak", {"models": ["henry,bad"]}, "opaque identifiers"),
    ],
)
def test_invalid_setting_domains_fail_closed(tmp_path, analysis, settings, message):
    fixture_name = (
        "synthetic_peak_shape" if analysis == "full_peak"
        else "synthetic_dispersive_profile" if analysis == "dispersive"
        else "synthetic_bet_isotherm"
    )
    fixture = bundled_contract_path() / "examples" / fixture_name
    data = _one_bundle_manifest(str(fixture), settings=settings)
    data["jobs"][0]["analysis"] = analysis
    manifest = _write_manifest(tmp_path / "batch.json", data)

    with pytest.raises(SystemExit, match=message):
        report.main(["--manifest", str(manifest), "--output", str(tmp_path / "out")])


def test_existing_output_is_never_overwritten(tmp_path):
    output = tmp_path / "out"
    output.mkdir()
    marker = output / "owner-data.txt"
    marker.write_text("preserve", encoding="utf-8")

    with pytest.raises(SystemExit, match="output already exists"):
        report.main(["--synthetic-example", "--output", str(output)])

    assert marker.read_text(encoding="utf-8") == "preserve"


def test_output_claim_preserves_unexpected_concurrent_content(tmp_path, monkeypatch):
    output = tmp_path / "out"
    original = report._write_batch_outputs

    def add_concurrent_content(staged, manifest, bundle_records, rows):
        original(staged, manifest, bundle_records, rows)
        (output / "owner-data.txt").write_text("preserve", encoding="utf-8")

    monkeypatch.setattr(report, "_write_batch_outputs", add_concurrent_content)

    with pytest.raises(SystemExit, match="changed during execution"):
        report.main(["--synthetic-example", "--output", str(output)])

    assert (output / "owner-data.txt").read_text(encoding="utf-8") == "preserve"


def test_output_cannot_be_created_inside_a_neutral_bundle(tmp_path):
    fixture = bundled_contract_path() / "examples" / "synthetic_bet_isotherm"
    bundle = tmp_path / "neutral"
    shutil.copytree(fixture, bundle)
    manifest = _write_manifest(
        tmp_path / "batch.json", _one_bundle_manifest(str(bundle))
    )
    output = bundle / "derived-output"

    with pytest.raises(SystemExit, match="must not be inside neutral bundle"):
        report.main(["--manifest", str(manifest), "--output", str(output)])

    assert not output.exists()
