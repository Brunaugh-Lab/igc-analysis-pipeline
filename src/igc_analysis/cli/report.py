"""Manifest-driven orchestration of source-neutral IGC analyses."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import re
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from igc_analysis import __version__
from igc_analysis.cli import acid_base, bet, dispersive
from igc_analysis.io.neutral_data import (
    NeutralContractError,
    bundled_contract_path,
    read_neutral_bundle,
)


BATCH_SCHEMA_VERSION = "igc-analysis-batch/0.1.0"
BATCH_CONTRACT_PATH = (
    Path(__file__).resolve().parents[1] / "contracts" / "igc-analysis-batch" / "0.1.0"
)
ID_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
SUPPORTED_ANALYSES = {"bet", "dispersive", "acid_base", "full_peak"}
TOP_LEVEL_KEYS = {"schema_version", "batch_id", "bundles", "jobs"}
BUNDLE_KEYS = {"bundle_id", "path"}
JOB_KEYS = {"job_id", "analysis", "bundle_ids", "settings"}

SETTING_KEYS = {
    "bet": {
        "probe", "retention", "concentration", "origin", "p0_min", "p0_max",
        "ambient_pressure_pa", "pressure_correction", "sensitivity",
    },
    "dispersive": {
        "homologous_probe_ids", "ambient_pressure_pa", "pressure_correction",
        "extrapolate", "max_temperature_span_k", "max_flow_relative_span",
    },
    "acid_base": {
        "homologous_probe_ids", "polar_probe_ids", "ambient_pressure_pa",
        "pressure_correction", "extrapolate", "max_temperature_span_k",
        "max_flow_relative_span",
    },
    "full_peak": {
        "probe", "transport_mode", "models", "n_cells", "n_starts", "lodo",
        "lodo_models", "cross_section_m2",
    },
}


class BatchManifestError(ValueError):
    """Raised when a batch manifest is ambiguous or invalid."""


def bundled_batch_contract_path() -> Path:
    """Return the packaged batch-contract directory."""

    return BATCH_CONTRACT_PATH


def _canonical_digest(value: object) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def _bundle_manifest_digest(manifest: dict) -> str:
    entries = [
        (filename, metadata["sha256"])
        for filename, metadata in sorted(manifest["files"].items())
    ]
    return _canonical_digest(entries)


def _require_identifier(value: object, field: str) -> str:
    if not isinstance(value, str) or not ID_PATTERN.fullmatch(value):
        raise BatchManifestError(
            f"{field} must be a nonempty opaque identifier containing only "
            "letters, numbers, periods, underscores, or hyphens"
        )
    return value


def _reject_unknown(mapping: dict, allowed: set[str], field: str) -> None:
    unknown = sorted(set(mapping) - allowed)
    if unknown:
        raise BatchManifestError(f"{field} contains unknown keys: {', '.join(unknown)}")


def _require_bool(settings: dict, key: str, default: bool) -> bool:
    value = settings.get(key, default)
    if not isinstance(value, bool):
        raise BatchManifestError(f"setting {key!r} must be true or false")
    return value


def _require_number(
    settings: dict, key: str, default: float | int, *, integer: bool = False,
) -> float | int:
    value = settings.get(key, default)
    valid = isinstance(value, int) if integer else isinstance(value, (int, float))
    if isinstance(value, bool) or not valid:
        kind = "integer" if integer else "number"
        raise BatchManifestError(f"setting {key!r} must be a {kind}")
    if not math.isfinite(float(value)):
        raise BatchManifestError(f"setting {key!r} must be finite")
    return int(value) if integer else float(value)


def _require_string_list(
    settings: dict, key: str, *, minimum: int = 0, optional: bool = False,
) -> list[str] | None:
    if key not in settings and optional:
        return None
    value = settings.get(key)
    if not isinstance(value, list) or any(not isinstance(item, str) for item in value):
        raise BatchManifestError(f"setting {key!r} must be a list of strings")
    if any(not ID_PATTERN.fullmatch(item) for item in value):
        raise BatchManifestError(
            f"setting {key!r} must contain only nonempty opaque identifiers"
        )
    if len(value) < minimum:
        raise BatchManifestError(
            f"setting {key!r} requires at least {minimum} explicit identifiers"
        )
    if len(set(value)) != len(value):
        raise BatchManifestError(f"setting {key!r} contains duplicate identifiers")
    return value


def _load_manifest(path: Path) -> dict:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise BatchManifestError(f"cannot read batch manifest: {exc}") from exc
    if not isinstance(value, dict):
        raise BatchManifestError("batch manifest must be a JSON object")
    return value


def _validate_manifest(data: dict, base_dir: Path) -> tuple[dict, dict, dict]:
    _reject_unknown(data, TOP_LEVEL_KEYS, "batch manifest")
    if data.get("schema_version") != BATCH_SCHEMA_VERSION:
        raise BatchManifestError(
            f"unsupported batch schema {data.get('schema_version')!r}; "
            f"expected {BATCH_SCHEMA_VERSION!r}"
        )
    batch_id = _require_identifier(data.get("batch_id"), "batch_id")
    bundles_value = data.get("bundles")
    jobs_value = data.get("jobs")
    if not isinstance(bundles_value, list) or not bundles_value:
        raise BatchManifestError("bundles must be a nonempty list")
    if not isinstance(jobs_value, list) or not jobs_value:
        raise BatchManifestError("jobs must be a nonempty list")

    bundle_paths: dict[str, Path] = {}
    bundle_records: dict[str, dict] = {}
    dataset_ids: set[str] = set()
    for index, item in enumerate(bundles_value):
        if not isinstance(item, dict):
            raise BatchManifestError(f"bundles[{index}] must be an object")
        _reject_unknown(item, BUNDLE_KEYS, f"bundles[{index}]")
        bundle_id = _require_identifier(item.get("bundle_id"), f"bundles[{index}].bundle_id")
        if bundle_id in bundle_paths:
            raise BatchManifestError(f"duplicate bundle_id {bundle_id!r}")
        raw_path = item.get("path")
        if not isinstance(raw_path, str) or not raw_path.strip():
            raise BatchManifestError(f"bundles[{index}].path must be a nonempty string")
        path = Path(raw_path).expanduser()
        if not path.is_absolute():
            path = base_dir / path
        try:
            neutral = read_neutral_bundle(path)
        except NeutralContractError as exc:
            raise BatchManifestError(f"bundle {bundle_id!r} is invalid: {exc}") from exc
        if neutral.dataset_id in dataset_ids:
            raise BatchManifestError(
                f"dataset_id {neutral.dataset_id!r} is declared by more than one bundle"
            )
        dataset_ids.add(neutral.dataset_id)
        bundle_paths[bundle_id] = neutral.path
        bundle_records[bundle_id] = {
            "dataset_id": neutral.dataset_id,
            "contract_version": neutral.contract_version,
            "manifest_digest": _bundle_manifest_digest(neutral.manifest),
            "source_fingerprint": neutral.manifest.get("source_fingerprint"),
        }

    jobs: list[dict] = []
    job_ids: set[str] = set()
    for index, item in enumerate(jobs_value):
        if not isinstance(item, dict):
            raise BatchManifestError(f"jobs[{index}] must be an object")
        _reject_unknown(item, JOB_KEYS, f"jobs[{index}]")
        job_id = _require_identifier(item.get("job_id"), f"jobs[{index}].job_id")
        if job_id in job_ids:
            raise BatchManifestError(f"duplicate job_id {job_id!r}")
        job_ids.add(job_id)
        analysis = item.get("analysis")
        if analysis not in SUPPORTED_ANALYSES:
            raise BatchManifestError(
                f"job {job_id!r} has unsupported analysis {analysis!r}; "
                f"choose from {', '.join(sorted(SUPPORTED_ANALYSES))}"
            )
        refs = item.get("bundle_ids")
        if not isinstance(refs, list) or not refs or any(not isinstance(ref, str) for ref in refs):
            raise BatchManifestError(f"job {job_id!r} bundle_ids must be a nonempty list")
        if len(set(refs)) != len(refs):
            raise BatchManifestError(f"job {job_id!r} contains duplicate bundle_ids")
        missing = [ref for ref in refs if ref not in bundle_paths]
        if missing:
            raise BatchManifestError(
                f"job {job_id!r} references unknown bundles: {', '.join(missing)}"
            )
        if len(refs) != 1:
            raise BatchManifestError(
                f"job {job_id!r} analysis {analysis!r} requires exactly one bundle"
            )
        settings = item.get("settings", {})
        if not isinstance(settings, dict):
            raise BatchManifestError(f"job {job_id!r} settings must be an object")
        _reject_unknown(settings, SETTING_KEYS[analysis], f"job {job_id!r} settings")
        jobs.append({
            "job_id": job_id,
            "analysis": analysis,
            "bundle_ids": refs,
            "settings": settings,
        })

    normalized = {
        "schema_version": BATCH_SCHEMA_VERSION,
        "batch_id": batch_id,
        "bundles": [
            {"bundle_id": bundle_id, **bundle_records[bundle_id]}
            for bundle_id in bundle_paths
        ],
        "jobs": jobs,
    }
    return normalized, bundle_paths, bundle_records


def _append_string(argv: list[str], flag: str, value: object, field: str) -> None:
    if not isinstance(value, str) or not value:
        raise BatchManifestError(f"setting {field!r} must be a nonempty string")
    argv.extend([flag, value])


def _bet_argv(bundle: Path, output: Path, settings: dict) -> list[str]:
    argv = ["--neutral-bundle", str(bundle), "--output", str(output)]
    for key, flag in (
        ("probe", "--probe"), ("retention", "--retention"),
        ("concentration", "--concentration"), ("origin", "--origin"),
    ):
        if key in settings:
            _append_string(argv, flag, settings[key], key)
    p0_min = _require_number(settings, "p0_min", 0.05)
    p0_max = _require_number(settings, "p0_max", 0.35)
    if not (0 <= p0_min < p0_max < 1):
        raise BatchManifestError("BET bounds must satisfy 0 <= p0_min < p0_max < 1")
    ambient = _require_number(settings, "ambient_pressure_pa", 101325.0)
    if ambient <= 0:
        raise BatchManifestError("setting 'ambient_pressure_pa' must be positive")
    if "p0_min" in settings:
        argv.extend(["--p0-min", str(p0_min)])
    if "p0_max" in settings:
        argv.extend(["--p0-max", str(p0_max)])
    if "ambient_pressure_pa" in settings:
        argv.extend(["--ambient-pressure-pa", str(ambient)])
    if not _require_bool(settings, "pressure_correction", True):
        argv.append("--no-pressure-correction")
    if not _require_bool(settings, "sensitivity", True):
        argv.append("--no-sensitivity")
    return argv


def _surface_argv(
    analysis: str, bundle: Path, output: Path, settings: dict,
) -> list[str]:
    argv = ["--neutral-bundle", str(bundle), "--output", str(output)]
    homologs = _require_string_list(
        settings, "homologous_probe_ids", minimum=3 if analysis == "acid_base" else 0,
        optional=analysis == "dispersive",
    )
    for probe_id in homologs or []:
        argv.extend(["--homologous-probe-id", probe_id])
    if analysis == "acid_base":
        polar = _require_string_list(settings, "polar_probe_ids", minimum=3)
        for probe_id in polar or []:
            argv.extend(["--polar-probe-id", probe_id])
    for key, flag, default in (
        ("ambient_pressure_pa", "--ambient-pressure-pa", 101325.0),
        ("max_temperature_span_k", "--max-temperature-span-k", 1.0),
        ("max_flow_relative_span", "--max-flow-relative-span", 0.05),
    ):
        if key in settings:
            value = _require_number(settings, key, default)
            if key == "ambient_pressure_pa" and value <= 0:
                raise BatchManifestError(f"setting {key!r} must be positive")
            if key != "ambient_pressure_pa" and value < 0:
                raise BatchManifestError(f"setting {key!r} cannot be negative")
            argv.extend([flag, str(value)])
    if not _require_bool(settings, "pressure_correction", True):
        argv.append("--no-pressure-correction")
    if not _require_bool(settings, "extrapolate", True):
        argv.append("--no-extrapolation")
    return argv


def _full_peak_argv(
    bundle_ids: list[str], bundle_paths: dict[str, Path], output: Path, settings: dict,
) -> list[str]:
    argv: list[str] = []
    for bundle_id in bundle_ids:
        argv.extend(["--neutral-bundle", f"{bundle_id}={bundle_paths[bundle_id]}"])
    argv.extend(["--output", str(output)])
    for key, flag in (("probe", "--probe"), ("transport_mode", "--transport-mode")):
        if key in settings:
            _append_string(argv, flag, settings[key], key)
    for key, flag, default in (
        ("n_cells", "--n-cells", 50), ("n_starts", "--n-starts", 4),
    ):
        if key in settings:
            value = _require_number(settings, key, default, integer=True)
            minimum = 2 if key == "n_cells" else 1
            if value < minimum:
                raise BatchManifestError(
                    f"setting {key!r} must be at least {minimum}"
                )
            argv.extend([flag, str(value)])
    for key, flag in (("models", "--models"), ("lodo_models", "--lodo-models")):
        if key in settings:
            values = _require_string_list(settings, key, minimum=1)
            argv.extend([flag, ",".join(values or [])])
    if not _require_bool(settings, "lodo", True):
        argv.append("--no-lodo")
    if "cross_section_m2" in settings:
        cross_section = _require_number(settings, "cross_section_m2", 0.0)
        if cross_section <= 0:
            raise BatchManifestError("setting 'cross_section_m2' must be positive")
        argv.extend([
            "--cross-section",
            str(cross_section),
        ])
    return argv


def _run_job(job: dict, bundle_paths: dict[str, Path], output: Path) -> dict:
    analysis = job["analysis"]
    settings = job["settings"]
    refs = job["bundle_ids"]
    if analysis == "bet":
        bet.main(_bet_argv(bundle_paths[refs[0]], output, settings))
        record_file = "bet_run.json"
    elif analysis == "dispersive":
        dispersive.main(_surface_argv(analysis, bundle_paths[refs[0]], output, settings))
        record_file = "dispersive_run.json"
    elif analysis == "acid_base":
        acid_base.main(_surface_argv(analysis, bundle_paths[refs[0]], output, settings))
        record_file = "acid_base_run.json"
    else:
        # full_peak sets Matplotlib font defaults when imported. Contain that
        # process-global state so job order cannot alter another CLI's figures.
        import matplotlib

        with matplotlib.rc_context({"pdf.fonttype": 42, "ps.fonttype": 42}):
            from igc_analysis.cli import full_peak

            full_peak.main(_full_peak_argv(refs, bundle_paths, output, settings))
        record_file = "full_peak_run.json"
    try:
        record = json.loads((output / record_file).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise BatchManifestError(
            f"{analysis} did not produce a readable {record_file}"
        ) from exc
    if not isinstance(record, dict):
        raise BatchManifestError(f"{analysis} produced an incompatible {record_file}")
    try:
        return _summarize_job_record(analysis, record_file, record)
    except (AttributeError, KeyError, TypeError, ValueError) as exc:
        raise BatchManifestError(
            f"{analysis} produced an incompatible {record_file}"
        ) from exc


def _summarize_job_record(analysis: str, record_file: str, record: dict) -> dict:
    """Extract stable batch fields from one child command's run record."""

    if analysis == "bet":
        reportable = bool(
            record.get("reportability")
            and record["reportability"].get("bet_applicable")
        )
        scope = "bet_ssa"
        qc_status = "PASS" if record["qc"]["passed"] else "REVIEW"
        qc_summary = "PASS" if record["qc"]["passed"] else "REVIEW"
        selected_model = None
    elif analysis in {"dispersive", "acid_base"}:
        reportable = bool(record["result"]["reportable"])
        scope = f"{analysis}_profile"
        qc_status = "PASS" if record["result"]["qc"]["pass"] else "REVIEW"
        qc_summary = record["result"]["qc"]["summary"]
        selected_model = None
    else:
        reportable = bool(record["ssa"]["reportable"])
        scope = "full_peak_recovered_ssa"
        qc_status = "NOT_COMBINED"
        qc_summary = (
            "all parameters identifiable: "
            f"{'yes' if record['qc']['all_params_identifiable'] else 'no'}; "
            f"minimum mass balance: {record['qc']['mass_balance_min']:.6g}; "
            f"cooperative: {'yes' if record['qc']['cooperative'] else 'no'}"
        )
        selected_model = record["selected_model"]
    return {
        "record_file": record_file,
        "reportability_scope": scope,
        "reportable": reportable,
        "qc_status": qc_status,
        "qc_summary": qc_summary,
        "selected_model": selected_model,
    }


def _write_batch_outputs(
    staged: Path, manifest: dict, bundle_records: dict[str, dict], rows: list[dict],
) -> None:
    csv_rows = [
        {
            **row,
            "bundle_ids": ";".join(row["bundle_ids"]),
            "dataset_ids": ";".join(row["dataset_ids"]),
        }
        for row in rows
    ]
    pd.DataFrame(csv_rows).to_csv(staged / "batch_summary.csv", index=False)
    record = {
        "settings": {
            "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "package_version": __version__,
            "python": platform.python_version(),
            "command_name": "igc-report",
            "batch_schema_version": BATCH_SCHEMA_VERSION,
            "orchestration_rule": (
                "explicit single-bundle jobs only; no filename-derived grouping, "
                "cross-bundle fitting, pooling, or "
                "cross-job scientific aggregation"
            ),
        },
        "input": {
            "batch_id": manifest["batch_id"],
            "manifest_digest": _canonical_digest(manifest),
            "bundles": [
                {"bundle_id": bundle_id, **bundle_records[bundle_id]}
                for bundle_id in bundle_records
            ],
        },
        "jobs": rows,
    }
    (staged / "batch_run.json").write_text(
        json.dumps(record, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    lines = [
        "# IGC batch result", "", f"- Batch: `{manifest['batch_id']}`",
        f"- Jobs completed: {len(rows)}", "",
        "Each job retains its own scientific reportability scope and output directory.",
        "This orchestrator accepts one bundle per job. It does not pool bundles,",
        "infer replicate groups from names,",
        "or turn several job verdicts into one batch-level scientific verdict.", "",
        "Review `batch_summary.csv`, `batch_run.json`, and every job's README, run",
        "record, QC flags, figures, and reportability verdict before reporting results.",
        "",
    ]
    (staged / "README.md").write_text("\n".join(lines), encoding="utf-8")


def _synthetic_manifest() -> dict:
    examples = bundled_contract_path() / "examples"
    return {
        "schema_version": BATCH_SCHEMA_VERSION,
        "batch_id": "batch-synthetic-001",
        "bundles": [
            {"bundle_id": "bundle-bet", "path": str(examples / "synthetic_bet_isotherm")},
            {"bundle_id": "bundle-surface", "path": str(examples / "synthetic_dispersive_profile")},
            {"bundle_id": "bundle-peak", "path": str(examples / "synthetic_peak_shape")},
        ],
        "jobs": [
            {"job_id": "bet-001", "analysis": "bet", "bundle_ids": ["bundle-bet"]},
            {
                "job_id": "dispersive-001", "analysis": "dispersive",
                "bundle_ids": ["bundle-surface"],
                "settings": {"homologous_probe_ids": [
                    "probe-homolog-08", "probe-homolog-09", "probe-homolog-10",
                ]},
            },
            {
                "job_id": "acid-base-001", "analysis": "acid_base",
                "bundle_ids": ["bundle-surface"],
                "settings": {
                    "homologous_probe_ids": [
                        "probe-homolog-08", "probe-homolog-09", "probe-homolog-10",
                    ],
                    "polar_probe_ids": [
                        "probe-polar-01", "probe-polar-02", "probe-polar-03",
                    ],
                },
            },
            {
                "job_id": "full-peak-001", "analysis": "full_peak",
                "bundle_ids": ["bundle-peak"],
                "settings": {
                    "models": ["none", "henry"], "n_cells": 12,
                    "n_starts": 1, "lodo": False,
                },
            },
        ],
    }


def _reject_output_inside_input(output: Path, bundle_paths: dict[str, Path]) -> None:
    resolved_output = output.expanduser().resolve()
    for bundle_id, bundle_path in bundle_paths.items():
        if resolved_output.is_relative_to(bundle_path.resolve()):
            raise BatchManifestError(
                f"output must not be inside neutral bundle {bundle_id!r}"
            )


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Run explicit source-neutral IGC analysis jobs as one batch"
    )
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument("--manifest", help="igc-analysis-batch/0.1.0 JSON manifest")
    inputs.add_argument(
        "--synthetic-example", action="store_true",
        help="run packaged BET, dispersive, acid/base, and full-peak fixtures",
    )
    parser.add_argument(
        "--output", "-o", default="igc_batch_results",
        help="new output directory; the command refuses an existing path",
    )
    args = parser.parse_args(argv)
    output = Path(args.output)
    if output.exists():
        raise SystemExit(f"igc-report: output already exists: {output}")
    if args.synthetic_example:
        raw_manifest = _synthetic_manifest()
        base_dir = Path.cwd()
    else:
        manifest_path = Path(args.manifest).expanduser().resolve()
        try:
            raw_manifest = _load_manifest(manifest_path)
        except BatchManifestError as exc:
            raise SystemExit(f"igc-report: {exc}") from exc
        base_dir = manifest_path.parent
    try:
        manifest, bundle_paths, bundle_records = _validate_manifest(
            raw_manifest, base_dir
        )
        _reject_output_inside_input(output, bundle_paths)
    except BatchManifestError as exc:
        raise SystemExit(f"igc-report: {exc}") from exc

    output.parent.mkdir(parents=True, exist_ok=True)
    try:
        output.mkdir()
    except FileExistsError as exc:
        raise SystemExit(f"igc-report: output already exists: {output}") from exc
    claimed_output = True
    claimed_mode = output.stat().st_mode & 0o777
    try:
        with tempfile.TemporaryDirectory(prefix=".igc-report-", dir=output.parent) as temporary:
            staged = Path(temporary)
            rows: list[dict] = []
            for job in manifest["jobs"]:
                job_output = staged / job["job_id"]
                try:
                    result = _run_job(job, bundle_paths, job_output)
                except (
                    BatchManifestError, NeutralContractError, ValueError, SystemExit,
                ) as exc:
                    raise BatchManifestError(
                        f"job {job['job_id']!r} ({job['analysis']}) failed: {exc}"
                    ) from exc
                rows.append({
                    "job_id": job["job_id"],
                    "analysis": job["analysis"],
                    "bundle_ids": job["bundle_ids"],
                    "dataset_ids": [
                        bundle_records[ref]["dataset_id"] for ref in job["bundle_ids"]
                    ],
                    "result_directory": job["job_id"],
                    **result,
                })
            _write_batch_outputs(staged, manifest, bundle_records, rows)
            staged.chmod(claimed_mode)
            try:
                output.rmdir()
            except OSError as exc:
                raise BatchManifestError(
                    "claimed output directory changed during execution; refusing "
                    "to replace it"
                ) from exc
            try:
                staged.replace(output)
            except OSError as exc:
                raise BatchManifestError(f"cannot publish batch output: {exc}") from exc
            claimed_output = False
    except BatchManifestError as exc:
        raise SystemExit(f"igc-report: {exc}") from exc
    finally:
        if claimed_output:
            try:
                output.rmdir()
            except OSError:
                # Preserve unexpected content instead of deleting it.
                pass

    print(f"Batch outputs written to {output}")
    print(f"Completed jobs: {len(manifest['jobs'])}")


if __name__ == "__main__":
    main()
