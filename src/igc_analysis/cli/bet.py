"""Command-line BET analysis for a validated source-neutral bundle."""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from igc_analysis.analysis.bet_workflow import run_bet_from_neutral
from igc_analysis import __version__
from igc_analysis.io.neutral_data import (
    NeutralContractError,
    bundled_contract_path,
    read_neutral_bundle,
)


SYNTHETIC_EXAMPLE = bundled_contract_path() / "examples" / "synthetic_bet_isotherm"


def _json_safe(value):
    if isinstance(value, dict):
        return {key: _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def _write_tables(result, output: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    injections = pd.DataFrame([asdict(item) for item in result.injections])
    isotherm = pd.DataFrame([asdict(item) for item in result.isotherm])
    injections.to_csv(output / "bet_injections.csv", index=False)
    isotherm.to_csv(output / "bet_isotherm.csv", index=False)
    pd.DataFrame({"P_over_P0": result.bet_x, "bet_y_g_per_mmol": result.bet_y}).to_csv(
        output / "bet_linearization.csv", index=False
    )
    return injections, isotherm


def _write_figure(result, output: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    if result.isotherm:
        axes[0].plot(
            [point.P_over_P0 for point in result.isotherm],
            [point.n_adsorbed_mmol_g for point in result.isotherm],
            "o-",
            color="tab:blue",
        )
    axes[0].set_xlabel("P/P0")
    axes[0].set_ylabel("adsorbed amount (mmol/g)")
    axes[0].set_title("Adsorption isotherm")

    if len(result.bet_x):
        axes[1].plot(result.bet_x, result.bet_y, "o", color="0.25", label="fit data")
        if np.isfinite(result.slope) and np.isfinite(result.intercept):
            grid = np.linspace(float(result.bet_x.min()), float(result.bet_x.max()), 100)
            axes[1].plot(
                grid,
                result.slope * grid + result.intercept,
                color="tab:red",
                label="linear fit",
            )
            axes[1].legend()
    axes[1].set_xlabel("P/P0")
    axes[1].set_ylabel("P/P0 / [n(1-P/P0)] (g/mmol)")
    axes[1].set_title("BET linearization")
    fig.tight_layout()
    fig.savefig(output / "bet_diagnostics.pdf", bbox_inches="tight")
    fig.savefig(output / "bet_diagnostics.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def _write_readme(result, record: dict, output: Path) -> None:
    classification = result.classification
    reportable = bool(classification and classification.bet_applicable)
    ssa = f"{result.ssa_m2_g:.4g} m2/g" if np.isfinite(result.ssa_m2_g) else "undefined"
    lines = [
        "# BET analysis result",
        "",
        f"- Dataset: `{record['input']['dataset_id']}`",
        f"- Probe: {result.probe}",
        f"- Numerical SSA: {ssa}",
        f"- Reportable BET SSA: {'yes' if reportable else 'no'}",
        f"- QC: {result.qc.flag_string}",
    ]
    if classification is not None:
        lines.extend([
            f"- Isotherm classification: {classification.isotherm_type}",
            f"- Rationale: {classification.rationale}",
            "",
            classification.recommendation,
        ])
    lines.extend([
        "",
        "The numerical fit is not permission to report an SSA. Use the strict",
        "reportability verdict, inspect `bet_diagnostics.pdf`, and review every",
        "QC message in `bet_run.json` before interpretation.",
        "",
        "Conventions used:",
        "",
        f"- Concentration: `{result.concentration_mode}`",
        f"- Retention: `{result.retention_mode}` with its matched dead-time convention",
        f"- Isotherm origin: `{result.origin_strategy}`",
        f"- Pressure correction: "
        f"{'applied' if record['settings']['pressure_correction'] else 'DISABLED'}",
        f"- Pressure basis: {', '.join(result.pressure_basis)}",
        "",
    ])
    if record["settings"]["sensitivity"] and result.diagnostics is not None:
        alt_retention = result.diagnostics.ssa_alt_retention
        alt_origin = result.diagnostics.ssa_alt_origin
        lines.extend([
            "Sensitivity calculations:",
            "",
            f"- `{result.diagnostics.alt_retention_mode}` retention SSA: "
            f"{alt_retention:.4g} m2/g" if np.isfinite(alt_retention)
            else f"- `{result.diagnostics.alt_retention_mode}` retention SSA: undefined",
            f"- `{result.diagnostics.alt_origin_strategy}` origin SSA: "
            f"{alt_origin:.4g} m2/g" if np.isfinite(alt_origin)
            else f"- `{result.diagnostics.alt_origin_strategy}` origin SSA: undefined",
            "",
        ])
    (output / "README.md").write_text("\n".join(lines), encoding="utf-8")


def _manifest_digest(manifest: dict) -> str:
    entries = [
        (filename, metadata["sha256"])
        for filename, metadata in sorted(manifest["files"].items())
    ]
    payload = json.dumps(entries, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="BET surface-area analysis from igc-neutral-data/0.2.0"
    )
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument("--neutral-bundle")
    inputs.add_argument(
        "--synthetic-example",
        action="store_true",
        help="run the packaged closed-form Type II verification fixture",
    )
    parser.add_argument("--output", "-o", default="bet_results")
    parser.add_argument("--probe", default="auto")
    parser.add_argument("--retention", choices=("peak_max", "cofm"), default="peak_max")
    parser.add_argument(
        "--concentration", choices=("eluted", "loop"), default="eluted",
        help="Default eluted is the corrected peak-apex convention; loop is legacy comparison only",
    )
    parser.add_argument(
        "--origin", choices=("legacy", "rectangular", "linear"), default="legacy"
    )
    parser.add_argument("--p0-min", type=float, default=0.05)
    parser.add_argument("--p0-max", type=float, default=0.35)
    parser.add_argument("--ambient-pressure-pa", type=float, default=101325.0)
    parser.add_argument("--no-pressure-correction", action="store_true")
    parser.add_argument("--no-sensitivity", action="store_true")
    args = parser.parse_args(argv)
    bundle_path = SYNTHETIC_EXAMPLE if args.synthetic_example else args.neutral_bundle

    try:
        neutral = read_neutral_bundle(bundle_path)
        result = run_bet_from_neutral(
            neutral,
            p0_min=args.p0_min,
            p0_max=args.p0_max,
            probe=args.probe,
            retention_mode=args.retention,
            concentration_mode=args.concentration,
            origin=args.origin,
            pressure_correction=not args.no_pressure_correction,
            ambient_pressure_pa=args.ambient_pressure_pa,
            sensitivity=not args.no_sensitivity,
        )
    except (NeutralContractError, ValueError) as exc:
        raise SystemExit(f"igc-bet: {exc}") from exc

    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    _write_tables(result, output)
    _write_figure(result, output)

    neutral_injections = neutral.table("injections.csv")
    neutral_components = neutral.table("injection_components.csv")
    neutral_properties = neutral.table("probe_properties.csv")
    neutral_calibrations = neutral.table("calibration.csv")
    classification = asdict(result.classification) if result.classification else None
    record = {
        "settings": {
            "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "package_version": __version__,
            "python": platform.python_version(),
            "command_name": "igc-bet",
            "p0_bounds": [args.p0_min, args.p0_max],
            "retention_mode": args.retention,
            "concentration_mode": args.concentration,
            "origin_strategy": args.origin,
            "pressure_correction": not args.no_pressure_correction,
            "ambient_pressure_pa": args.ambient_pressure_pa,
            "sensitivity": not args.no_sensitivity,
        },
        "input": {
            "dataset_id": neutral.dataset_id,
            "contract_version": neutral.contract_version,
            "source_fingerprint": neutral.manifest.get("source_fingerprint"),
            "manifest_digest": _manifest_digest(neutral.manifest),
            "block_ids": sorted(set(neutral_injections["block_id"].astype(str))),
            "dead_time_injection_ids": list(
                neutral_injections.loc[
                    neutral_injections["role"] == "dead_time", "injection_id"
                ].astype(str)
            ),
        },
        "declared_provenance": {
            "probe_property_sources": sorted(
                set(neutral_properties["properties_source"].astype(str))
            ),
            "vapor_pressure_sources": sorted(
                value for value in set(
                    neutral_components["vapor_pressure_source"].astype(str)
                ) if value
            ),
            "vapor_pressure_model_ids": sorted(
                value for value in set(
                    neutral_components["vapor_pressure_model_id"].astype(str)
                ) if value
            ),
            "calibration_sources": sorted(
                set(neutral_calibrations["calibration_source"].astype(str))
            ),
        },
        "result": {
            "ssa_m2_g": result.ssa_m2_g,
            "n_monolayer_mmol_g": result.n_monolayer_mmol_g,
            "C_bet": result.C_bet,
            "r_squared": result.r_squared,
            "n_points": result.n_points,
            "fit_range": result.p_over_p0_range if result.n_points >= 2 else None,
            "probe": result.probe,
            "sample_id": result.sample_name,
            "cross_section_m2": result.cross_section_m2,
            "mean_p_sat_Pa": result.p_sat_Pa,
            "mean_james_martin_j": result.james_martin_j,
            "mean_standard_flow_sccm": result.flow_sccm,
            "flow_source_channels": result.flow_source_channels,
            "conditions_source": result.conditions_source,
            "pressure_source": result.pressure_source,
            "pressure_basis": result.pressure_basis,
        },
        "reportability": classification,
        "qc": {
            "passed": result.qc.passed,
            "flags": result.qc.flags,
            "messages": result.qc.messages,
        },
        "diagnostics": asdict(result.diagnostics) if result.diagnostics else None,
    }
    safe_record = _json_safe(record)
    (output / "bet_run.json").write_text(
        json.dumps(safe_record, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    _write_readme(result, safe_record, output)

    reportable = bool(result.classification and result.classification.bet_applicable)
    print(f"BET outputs written to {output}")
    print(f"QC: {result.qc.flag_string}")
    print(f"Reportable BET SSA: {'yes' if reportable else 'no'}")


if __name__ == "__main__":
    main()
