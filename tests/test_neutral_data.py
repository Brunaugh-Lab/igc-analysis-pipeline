"""Source-neutral contract ingestion and peak-shape integration tests."""

from __future__ import annotations

import json
import shutil
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from igc_analysis.analysis.full_peak import (
    build_trace_dataset_from_neutral,
    traces_to_dataframe,
)
from igc_analysis.io.neutral_data import (
    NeutralContractError,
    bundled_contract_path,
    read_neutral_bundle,
)
from igc_analysis.cli.full_peak import _json_safe, _parse_neutral_bundles


@pytest.fixture()
def synthetic_bundle() -> Path:
    return bundled_contract_path() / "examples" / "synthetic_peak_shape"


def test_valid_020_bundle_is_loaded_after_structural_validation(synthetic_bundle: Path):
    bundle = read_neutral_bundle(synthetic_bundle)
    assert bundle.contract_version == "0.2.0"
    assert bundle.dataset_id == "synthetic-peak-shape-001"
    assert set(bundle.tables) == set(bundle.manifest["files"])
    assert len(bundle.table("injections.csv")) == 4
    assert len(bundle.table("traces.csv")) == 60


def test_manifest_hash_tampering_is_rejected(synthetic_bundle: Path, tmp_path: Path):
    copied = tmp_path / "bundle"
    shutil.copytree(synthetic_bundle, copied)
    traces = copied / "traces.csv"
    traces.write_text(
        traces.read_text(encoding="utf-8").replace(",0.10,arbitrary_unit", ",9.99,arbitrary_unit", 1),
        encoding="utf-8",
    )
    with pytest.raises(NeutralContractError, match="SHA-256 digest does not match"):
        read_neutral_bundle(copied)


def test_wrong_contract_version_is_rejected(synthetic_bundle: Path, tmp_path: Path):
    copied = tmp_path / "bundle"
    shutil.copytree(synthetic_bundle, copied)
    manifest_path = copied / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["contract_version"] = "0.1.0"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    with pytest.raises(NeutralContractError, match="contract_version"):
        read_neutral_bundle(copied)


def test_neutral_input_labels_must_be_opaque():
    assert _parse_neutral_bundles(
        SimpleNamespace(neutral_bundle=["block-001=/governed/location"])
    ) == {"block-001": "/governed/location"}
    with pytest.raises(SystemExit, match="neutral opaque identifier"):
        _parse_neutral_bundles(
            SimpleNamespace(neutral_bundle=["non neutral label=/governed/location"])
        )


def test_run_record_values_are_strict_json_safe():
    assert _json_safe({"nan": np.nan, "inf": np.inf, "value": np.float64(2.0)}) == {
        "nan": None,
        "inf": None,
        "value": 2.0,
    }


def test_neutral_bundle_drives_full_peak_trace_construction(synthetic_bundle: Path):
    blocks = build_trace_dataset_from_neutral(
        {"synthetic-block": synthetic_bundle},
        n_cells=60,
        transport_mode="bracket_interpolated",
        verbose=False,
    )
    assert len(blocks) == 1
    block = blocks[0]
    assert block.block == "synthetic-block"
    assert block.schedule == [
        "injection-000",
        "injection-001",
        "injection-002",
        "injection-003",
    ]
    assert block.methane_names == ["injection-000", "injection-003"]
    assert [injection.name for injection in block.injections] == [
        "injection-001",
        "injection-002",
    ]
    assert all(injection.probe == "SYNTHETIC PROBE ALPHA" for injection in block.injections)
    assert all(injection.n_injected_mol > 0 for injection in block.injections)
    assert all(np.all(injection.c_out_mol_m3 >= 0) for injection in block.injections)
    assert all(injection.t0_assignment_basis.startswith("neutral") for injection in block.injections)

    traces = traces_to_dataframe(blocks)
    assert len(traces) == 30
    assert set(traces["injection_name"]) == {"injection-001", "injection-002"}
