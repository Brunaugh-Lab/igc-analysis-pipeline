"""Tests for explicit source-neutral batch orchestration."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

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


def test_synthetic_batch_runs_every_neutral_workflow_transactionally(tmp_path):
    output = tmp_path / "results"

    report.main(["--synthetic-example", "--output", str(output)])

    summary = pd.read_csv(output / "batch_summary.csv", keep_default_na=False)
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
    record = json.loads((output / "batch_run.json").read_text(encoding="utf-8"))
    assert record["settings"]["batch_schema_version"] == report.BATCH_SCHEMA_VERSION
    assert record["input"]["batch_id"] == "batch-synthetic-001"
    assert len(record["jobs"]) == 4
    assert all("path" not in bundle for bundle in record["input"]["bundles"])
    assert str(tmp_path) not in json.dumps(record)
    for job_id in ("bet-001", "dispersive-001", "acid-base-001", "full-peak-001"):
        assert (output / job_id / "README.md").is_file()


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


def test_existing_output_is_never_overwritten(tmp_path):
    output = tmp_path / "out"
    output.mkdir()
    marker = output / "owner-data.txt"
    marker.write_text("preserve", encoding="utf-8")

    with pytest.raises(SystemExit, match="output already exists"):
        report.main(["--synthetic-example", "--output", str(output)])

    assert marker.read_text(encoding="utf-8") == "preserve"


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
