#!/usr/bin/env python3
"""Validate an IGC neutral-data 0.2.0 bundle using the standard library."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any


SPEC_DIR = Path(__file__).resolve().parent
SCHEMA_PATH = SPEC_DIR / "schema.json"
ID_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$")
SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")
WINDOWS_ABSOLUTE_PATH = re.compile(r"^[A-Za-z]:[\\/]")
SENTINELS = {"na", "n/a", "null", "none", "-999", "."}


class ValidationError(Exception):
    """Raised when a bundle does not conform to the contract."""


def _fail(message: str) -> None:
    raise ValidationError(message)


def _read_json(path: Path) -> dict[str, Any]:
    try:
        with path.open("r", encoding="utf-8-sig") as handle:
            value = json.load(handle)
    except (OSError, json.JSONDecodeError) as exc:
        _fail(f"cannot read valid JSON from {path}: {exc}")
    if not isinstance(value, dict):
        _fail(f"expected a JSON object in {path}")
    return value


def _parse_datetime(value: str, location: str) -> str:
    normalized = value[:-1] + "+00:00" if value.endswith("Z") else value
    try:
        parsed = datetime.fromisoformat(normalized)
    except ValueError:
        _fail(f"{location}: invalid ISO 8601 datetime {value!r}")
    if parsed.tzinfo is None or parsed.utcoffset() is None:
        _fail(f"{location}: datetime must include a timezone")
    return value


def _parse_value(value: str | None, definition: dict[str, Any], location: str) -> Any:
    if value is None:
        _fail(f"{location}: row has fewer fields than the header")
    if value == "":
        if definition.get("required", False):
            _fail(f"{location}: required value is empty")
        return None
    stripped = value.strip()
    if stripped == "":
        _fail(f"{location}: whitespace-only values are not allowed")
    if stripped != value:
        _fail(f"{location}: leading or trailing whitespace is not allowed")
    if stripped.lower() in SENTINELS:
        _fail(f"{location}: sentinel missing value {value!r} is not allowed")

    value_type = definition["type"]
    try:
        if value_type == "string":
            parsed: Any = value
        elif value_type == "integer":
            if not re.fullmatch(r"[+-]?\d+", value):
                raise ValueError
            parsed = int(value)
        elif value_type == "number":
            parsed = float(value)
            if not math.isfinite(parsed):
                raise ValueError
        elif value_type == "boolean":
            if value.lower() not in {"true", "false"}:
                raise ValueError
            parsed = value.lower() == "true"
        elif value_type == "datetime":
            parsed = _parse_datetime(value, location)
        else:
            _fail(f"{location}: unsupported schema type {value_type!r}")
    except ValueError:
        _fail(f"{location}: {value!r} is not a valid {value_type}")

    if "enum" in definition and parsed not in definition["enum"]:
        _fail(f"{location}: {parsed!r} is not an allowed value")
    if definition.get("format") == "identifier" and not ID_PATTERN.fullmatch(parsed):
        _fail(f"{location}: value must be a neutral opaque identifier")
    if value_type == "string" and (
        parsed.startswith(("/", "~/", "file://")) or WINDOWS_ABSOLUTE_PATH.match(parsed)
    ):
        _fail(f"{location}: absolute or user-relative paths are not allowed")
    if parsed is not None and value_type in {"integer", "number"}:
        if "minimum" in definition and parsed < definition["minimum"]:
            _fail(f"{location}: value must be >= {definition['minimum']}")
        if "maximum" in definition and parsed > definition["maximum"]:
            _fail(f"{location}: value must be <= {definition['maximum']}")
        if "exclusive_minimum" in definition and parsed <= definition["exclusive_minimum"]:
            _fail(f"{location}: value must be > {definition['exclusive_minimum']}")
    return parsed


def _read_table(path: Path, definition: dict[str, Any]) -> list[dict[str, Any]]:
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, restkey="__extra__", restval=None)
            expected = list(definition["columns"])
            if reader.fieldnames != expected:
                _fail(f"{path.name}: header must be exactly {expected}; found {reader.fieldnames}")
            rows: list[dict[str, Any]] = []
            for line_number, raw_row in enumerate(reader, start=2):
                if raw_row.get("__extra__"):
                    _fail(f"{path.name}:{line_number}: row has more fields than the header")
                parsed_row: dict[str, Any] = {}
                for column, column_definition in definition["columns"].items():
                    parsed_row[column] = _parse_value(
                        raw_row[column], column_definition, f"{path.name}:{line_number}:{column}"
                    )
                rows.append(parsed_row)
    except UnicodeDecodeError as exc:
        _fail(f"{path.name}: file must be UTF-8: {exc}")
    except OSError as exc:
        _fail(f"cannot read {path}: {exc}")
    return rows


def _key(row: dict[str, Any], columns: list[str]) -> tuple[Any, ...]:
    return tuple(row[column] for column in columns)


def _check_unique(table_name: str, rows: list[dict[str, Any]], columns: list[str]) -> None:
    seen: set[tuple[Any, ...]] = set()
    for row_number, row in enumerate(rows, start=2):
        key = _key(row, columns)
        if key in seen:
            _fail(f"{table_name}:{row_number}: duplicate key {columns}={key}")
        seen.add(key)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _validate_manifest(
    bundle: Path, schema: dict[str, Any]
) -> tuple[dict[str, Any], list[str]]:
    manifest = _read_json(bundle / "manifest.json")
    required_manifest_fields = set(schema["manifest"]["required"])
    allowed_manifest_fields = required_manifest_fields | set(schema["manifest"].get("optional", []))
    for field in required_manifest_fields:
        if field not in manifest:
            _fail(f"manifest.json: missing required field {field!r}")
    unknown_manifest_fields = sorted(set(manifest) - allowed_manifest_fields)
    if unknown_manifest_fields:
        _fail(f"manifest.json: undeclared fields are not allowed: {unknown_manifest_fields}")
    if manifest["contract_name"] != schema["contract_name"]:
        _fail("manifest.json: contract_name does not match this specification")
    if manifest["contract_version"] != schema["contract_version"]:
        _fail("manifest.json: contract_version does not match this specification")
    profile = manifest["profile"]
    if profile not in schema["profiles"]:
        _fail(f"manifest.json: unknown profile {profile!r}")
    if not isinstance(manifest["dataset_id"], str) or not ID_PATTERN.fullmatch(manifest["dataset_id"]):
        _fail("manifest.json: dataset_id must be a neutral opaque identifier")
    _parse_datetime(str(manifest["created_at"]), "manifest.json:created_at")
    if not isinstance(manifest["adapter_version"], str) or not manifest["adapter_version"]:
        _fail("manifest.json: adapter_version must be a nonempty string")
    if not ID_PATTERN.fullmatch(manifest["adapter_version"]):
        _fail("manifest.json: adapter_version must be a neutral opaque version identifier")
    source_fingerprint = manifest.get("source_fingerprint")
    if source_fingerprint is not None and (
        not isinstance(source_fingerprint, str) or not SHA256_PATTERN.fullmatch(source_fingerprint)
    ):
        _fail("manifest.json: source_fingerprint must be a lowercase SHA-256 digest")

    required_tables = schema["profiles"][profile]["required_tables"]
    expected_entries = {"manifest.json", *required_tables}
    actual_entries = {path.name for path in bundle.iterdir()}
    unknown_entries = sorted(actual_entries - expected_entries)
    if unknown_entries:
        _fail(f"bundle contains undeclared entries: {unknown_entries}")
    file_metadata = manifest["files"]
    if not isinstance(file_metadata, dict) or set(file_metadata) != set(required_tables):
        _fail("manifest.json: files must contain exactly the tables required by the profile")
    for filename in required_tables:
        path = bundle / filename
        if not path.is_file():
            _fail(f"bundle is missing required file {filename}")
        metadata = file_metadata[filename]
        if not isinstance(metadata, dict) or set(metadata) != {"row_count", "sha256"}:
            _fail(f"manifest.json: {filename} metadata must contain row_count and sha256")
        if not isinstance(metadata["row_count"], int) or metadata["row_count"] < 0:
            _fail(f"manifest.json: {filename} row_count must be a nonnegative integer")
        digest = metadata["sha256"]
        if not isinstance(digest, str) or not SHA256_PATTERN.fullmatch(digest):
            _fail(f"manifest.json: {filename} sha256 must be a lowercase SHA-256 digest")
        if _sha256(path) != digest:
            _fail(f"manifest.json: {filename} SHA-256 digest does not match")
    return manifest, required_tables


def _validate_cross_table(
    manifest: dict[str, Any], schema: dict[str, Any], tables: dict[str, list[dict[str, Any]]]
) -> None:
    for table_name, rows in tables.items():
        definition = schema["tables"][table_name]
        _check_unique(table_name, rows, definition["primary_key"])
        for unique_columns in definition.get("unique", []):
            _check_unique(table_name, rows, unique_columns)

    for table_name, rows in tables.items():
        for foreign_key in schema["tables"][table_name].get("foreign_keys", []):
            reference_rows = tables[foreign_key["reference_table"]]
            valid_keys = {_key(row, foreign_key["reference_columns"]) for row in reference_rows}
            for row_number, row in enumerate(rows, start=2):
                candidate = _key(row, foreign_key["columns"])
                if foreign_key.get("allow_empty") and all(value is None for value in candidate):
                    continue
                if candidate not in valid_keys:
                    _fail(
                        f"{table_name}:{row_number}: unresolved foreign key "
                        f"{foreign_key['columns']}={candidate}"
                    )

    experiments = tables["experiment.csv"]
    if len(experiments) != 1:
        _fail("experiment.csv: version 0.2.0 requires exactly one experiment row")
    experiment = experiments[0]
    if experiment["dataset_id"] != manifest["dataset_id"]:
        _fail("experiment.csv: dataset_id must match manifest.json")
    if (experiment["specific_surface_area_m2_g"] is None) != (experiment["surface_area_source"] is None):
        _fail("experiment.csv: SSA and surface_area_source must occur together")

    for row_number, column in enumerate(tables["columns.csv"], start=2):
        if (column["density_kg_m3"] is None) != (column["density_basis"] is None):
            _fail(f"columns.csv:{row_number}: density and density_basis must occur together")

    injections = tables["injections.csv"]
    sequence = [row["sequence_index"] for row in injections]
    if sequence != list(range(len(sequence))):
        _fail("injections.csv: rows must be in contiguous zero-based sequence_index order")

    components_by_injection: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row_number, component in enumerate(tables["injection_components.csv"], start=2):
        components_by_injection[component["injection_id"]].append(component)
        vapor_fields = [
            component["saturation_vapor_pressure_Pa"],
            component["vapor_pressure_source"],
            component["vapor_pressure_model_id"],
        ]
        if any(value is not None for value in vapor_fields) and not all(value is not None for value in vapor_fields):
            _fail(f"injection_components.csv:{row_number}: vapor-pressure fields must occur together")

    injection_by_id = {row["injection_id"]: row for row in injections}
    traces_by_key: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row_number, trace in enumerate(tables["traces.csv"], start=2):
        declared_channel = injection_by_id[trace["injection_id"]]["detector_channel"]
        if trace["detector_channel"] != declared_channel:
            _fail(f"traces.csv:{row_number}: detector channel does not match the injection")
        traces_by_key[(trace["injection_id"], trace["detector_channel"])].append(trace)
        corrected_fields = [
            trace["signal_corrected"],
            trace["preprocessing_method"],
            trace["preprocessing_version"],
        ]
        if any(value is not None for value in corrected_fields) and not all(
            value is not None for value in corrected_fields
        ):
            _fail(f"traces.csv:{row_number}: corrected signal and both preprocessing fields must occur together")

    for injection in injections:
        injection_id = injection["injection_id"]
        if injection_id not in components_by_injection:
            _fail(f"injections.csv: {injection_id} has no injection component")
        trace_key = (injection_id, injection["detector_channel"])
        if trace_key not in traces_by_key:
            _fail(f"injections.csv: {injection_id} has no trace for its declared detector channel")
        roles = {component["component_role"] for component in components_by_injection[injection_id]}
        if injection["role"] == "dead_time" and "dead_time_marker" not in roles:
            _fail(f"injections.csv: dead-time injection {injection_id} lacks a dead_time_marker component")
        if injection["role"] == "probe" and "analyte" not in roles:
            _fail(f"injections.csv: probe injection {injection_id} lacks an analyte component")

    for (injection_id, channel), points in traces_by_key.items():
        indices = [point["point_index"] for point in points]
        if indices != list(range(len(indices))):
            _fail(f"traces.csv: {injection_id}/{channel} point_index must be contiguous and ordered from zero")
        times = [point["time_s"] for point in points]
        if any(current <= previous for previous, current in zip(times, times[1:])):
            _fail(f"traces.csv: {injection_id}/{channel} time_s must be strictly increasing")
        signal_units = {point["signal_unit"] for point in points}
        if len(signal_units) != 1:
            _fail(f"traces.csv: {injection_id}/{channel} must use one signal_unit")

    allowed_units = schema["tables"]["conditions.csv"]["unit_by_quantity"]
    for row_number, condition in enumerate(tables["conditions.csv"], start=2):
        if condition["unit"] != allowed_units[condition["quantity"]]:
            _fail(f"conditions.csv:{row_number}: unit is incompatible with quantity")
        if condition["quantity"] == "relative_humidity" and not 0 <= condition["value"] <= 1:
            _fail(f"conditions.csv:{row_number}: relative_humidity must be between 0 and 1")
        if condition["quantity"] == "column_temperature" and condition["value"] <= 0:
            _fail(f"conditions.csv:{row_number}: column temperature must be positive")
        if condition["quantity"].startswith("flow_") and condition["value"] <= 0:
            _fail(f"conditions.csv:{row_number}: flow must be positive")
        if (
            condition["quantity"].startswith("flow_")
            and condition["value_role"] == "measured"
            and condition["source_channel"] is None
        ):
            _fail(f"conditions.csv:{row_number}: measured flow requires source_channel")
        if condition["quantity"] in {"pressure_inlet", "pressure_outlet"} and condition["value"] <= 0:
            _fail(f"conditions.csv:{row_number}: absolute pressure must be positive")
        if condition["quantity"] == "pressure_drop" and condition["value"] < 0:
            _fail(f"conditions.csv:{row_number}: pressure drop must be nonnegative")

    calibration_by_id = {
        row["calibration_id"]: row for row in tables["calibration.csv"]
    }
    for row_number, component in enumerate(tables["injection_components.csv"], start=2):
        calibration_id = component["calibration_id"]
        if calibration_id is not None and calibration_by_id[calibration_id]["probe_id"] != component["probe_id"]:
            _fail(f"injection_components.csv:{row_number}: calibration probe does not match component probe")

    for row_number, calibration in enumerate(tables["calibration.csv"], start=2):
        quadratic = calibration["calibration_model"] == "quadratic"
        if quadratic != (calibration["parameter_2"] is not None):
            _fail(f"calibration.csv:{row_number}: parameter_2 must occur only for a quadratic model")


def validate_bundle(bundle: Path) -> dict[str, int]:
    bundle = bundle.expanduser().resolve()
    if not bundle.is_dir():
        _fail(f"bundle directory does not exist: {bundle}")
    schema = _read_json(SCHEMA_PATH)
    manifest, required_tables = _validate_manifest(bundle, schema)
    tables: dict[str, list[dict[str, Any]]] = {}
    for table_name in required_tables:
        rows = _read_table(bundle / table_name, schema["tables"][table_name])
        expected_count = manifest["files"][table_name]["row_count"]
        if len(rows) != expected_count:
            _fail(f"{table_name}: row count {len(rows)} does not match manifest value {expected_count}")
        tables[table_name] = rows
    _validate_cross_table(manifest, schema, tables)
    return {table_name: len(rows) for table_name, rows in tables.items()}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bundle", type=Path, help="normalized bundle directory")
    args = parser.parse_args()
    try:
        counts = validate_bundle(args.bundle)
    except ValidationError as exc:
        print(f"INVALID: {exc}", file=sys.stderr)
        return 1
    total_rows = sum(counts.values())
    print(f"VALID igc-neutral-data 0.2.0 trace-core bundle: {total_rows} rows")
    for table_name, count in counts.items():
        print(f"  {table_name}: {count}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
