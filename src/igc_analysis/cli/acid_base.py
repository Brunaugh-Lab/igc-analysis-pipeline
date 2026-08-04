"""Command-line Schultz--Gutmann acid/base analysis for a neutral bundle."""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from igc_analysis import __version__
from igc_analysis.analysis.acid_base_workflow import run_acid_base_from_neutral
from igc_analysis.io.neutral_data import (
    NeutralContractError,
    bundled_contract_path,
    read_neutral_bundle,
)


SYNTHETIC_EXAMPLE = (
    bundled_contract_path() / "examples" / "synthetic_dispersive_profile"
)
SYNTHETIC_HOMOLOGS = (
    "probe-homolog-08", "probe-homolog-09", "probe-homolog-10",
)
SYNTHETIC_POLAR = ("probe-polar-01", "probe-polar-02", "probe-polar-03")


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


def _manifest_digest(manifest: dict) -> str:
    entries = [
        (filename, metadata["sha256"])
        for filename, metadata in sorted(manifest["files"].items())
    ]
    payload = json.dumps(entries, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def _save(fig, output: Path, stem: str) -> None:
    fig.savefig(output / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(output / f"{stem}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def _write_figures(result, output: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
    for mode, style, label in (
        ("cofm", "o-", "center of mass (primary)"),
        ("peak_max", "s--", "peak maximum (sensitivity)"),
    ):
        data = result.profile[result.profile["retention_mode"] == mode].sort_values(
            "coverage"
        )
        axes[0].plot(data["coverage"], data["Ka"], style, label=label)
        axes[1].plot(data["coverage"], data["Kb"], style, label=label)
    for ax, label in zip(axes, ("Ka", "Kb")):
        ax.set_xlabel("fractional surface coverage")
        ax.set_ylabel(label)
        ax.grid(alpha=0.25)
        ax.legend()
    fig.suptitle("Coverage-resolved Gutmann acid/base parameters")
    fig.tight_layout()
    _save(fig, output, "acid_base_profile")

    primary = result.delta_g_sp[
        result.delta_g_sp["retention_mode"] == "cofm"
    ]
    coverages = sorted(primary["coverage"].unique())
    fig, axes = plt.subplots(
        1, len(coverages), figsize=(4.0 * len(coverages), 3.6), squeeze=False
    )
    for index, coverage in enumerate(coverages):
        ax = axes[0][index]
        data = primary[np.isclose(primary["coverage"], coverage)].copy()
        x = data["dn"] / data["an_star"]
        y = data["delta_g_sp_kJmol"] / data["an_star"]
        ax.plot(x, y, "o")
        profile = result.profile[
            (result.profile["retention_mode"] == "cofm")
            & np.isclose(result.profile["coverage"], coverage)
        ].iloc[0]
        if len(data) >= 2 and np.isfinite(profile["Ka"]):
            grid = np.linspace(float(x.min()), float(x.max()), 100)
            ax.plot(grid, profile["Ka"] * grid + profile["Kb"], "--", color="0.25")
        ax.set_title(f"coverage={coverage:.3g}; R2={profile['r_squared']:.3f}")
        ax.set_xlabel("DN / AN*")
        ax.set_ylabel("delta Gsp / AN*")
        ax.grid(alpha=0.2)
    fig.tight_layout()
    _save(fig, output, "acid_base_gutmann_fits")


def _write_readme(result, output: Path) -> None:
    lines = [
        "# Acid/base surface-characterization result",
        "",
        f"- Dataset: `{result.dataset_id}`",
        f"- Reportable profile: {'yes' if result.reportable else 'no'}",
        f"- QC: {result.qc['summary']}",
        "",
        "Center-of-mass retention is primary. Peak maximum reproduces the",
        "historical retention convention as a sensitivity calculation.",
        "",
        "Ka and Kb are convention-dependent Gutmann descriptors. They are not",
        "a unique molecular mechanism or interchangeable with van Oss surface-",
        "energy components. Review the declared probe-property sources, every",
        "coverage mapping, the per-probe delta Gsp table, and all QC flags before",
        "reporting the profile.",
        "",
        "Van Oss analysis is intentionally unavailable in this source-neutral",
        "workflow until its additional liquid-component properties and sources",
        "are represented explicitly in the neutral contract.",
        "",
    ]
    (output / "README.md").write_text("\n".join(lines), encoding="utf-8")


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Schultz--Gutmann acid/base analysis from igc-neutral-data/0.2.0"
    )
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument("--neutral-bundle")
    inputs.add_argument("--synthetic-example", action="store_true")
    parser.add_argument("--homologous-probe-id", action="append")
    parser.add_argument("--polar-probe-id", action="append")
    parser.add_argument("--output", "-o", default="acid_base_results")
    parser.add_argument("--ambient-pressure-pa", type=float, default=101325.0)
    parser.add_argument("--no-pressure-correction", action="store_true")
    parser.add_argument("--no-extrapolation", action="store_true")
    parser.add_argument("--max-temperature-span-k", type=float, default=1.0)
    parser.add_argument("--max-flow-relative-span", type=float, default=0.05)
    args = parser.parse_args(argv)
    if args.synthetic_example:
        bundle_path = SYNTHETIC_EXAMPLE
        homologs = args.homologous_probe_id or list(SYNTHETIC_HOMOLOGS)
        polar = args.polar_probe_id or list(SYNTHETIC_POLAR)
    else:
        bundle_path = args.neutral_bundle
        homologs = args.homologous_probe_id
        polar = args.polar_probe_id
    if not homologs or not polar:
        parser.error(
            "real bundles require repeated --homologous-probe-id and "
            "--polar-probe-id declarations"
        )
    try:
        neutral = read_neutral_bundle(bundle_path)
        result = run_acid_base_from_neutral(
            neutral,
            homologous_probe_ids=homologs,
            polar_probe_ids=polar,
            pressure_correction=not args.no_pressure_correction,
            ambient_pressure_pa=args.ambient_pressure_pa,
            extrapolate=not args.no_extrapolation,
            max_temperature_span_K=args.max_temperature_span_k,
            max_flow_relative_span=args.max_flow_relative_span,
        )
    except (NeutralContractError, ValueError) as exc:
        raise SystemExit(f"igc-acid-base: {exc}") from exc

    record = _json_safe({
        "settings": {
            "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "package_version": __version__,
            "python": platform.python_version(),
            "command_name": "igc-acid-base",
            "primary_retention_mode": "cofm",
            "sensitivity_retention_mode": "peak_max",
            "homologous_probe_ids": list(homologs),
            "polar_probe_ids": list(polar),
            "polar_probe_policy": "explicit_inclusion",
            "pressure_correction": not args.no_pressure_correction,
            "ambient_pressure_pa": args.ambient_pressure_pa,
            "extrapolate": not args.no_extrapolation,
            "van_oss": "not_available_without_explicit_contract_properties",
            "reportability_rule": (
                "zero critical flags; at least three homologs and three polar "
                "probes at every coverage; finite regression-derived Ka, Kb, "
                "and R2; no selected-probe coverage extrapolation; one detector gain"
            ),
        },
        "input": {
            "dataset_id": neutral.dataset_id,
            "manifest_fingerprint": _manifest_digest(neutral.manifest),
            "contract_name": neutral.manifest["contract_name"],
            "contract_version": neutral.manifest["contract_version"],
            "sample_id": result.sample_id,
            "properties_sources": result.properties_sources,
            "calibration_sources": result.calibration_sources,
            "flow_source_channels": result.flow_source_channels,
            "pressure_basis": result.pressure_basis,
            "pressure_roles": result.pressure_roles,
            "detector_gains": result.detector_gains,
        },
        "result": {
            "reportable": result.reportable,
            "temperature_K": result.temperature_K,
            "coverage_points": int((result.profile["retention_mode"] == "cofm").sum()),
            "qc": result.qc,
        },
    })
    output = Path(args.output)
    if output.exists():
        raise SystemExit(f"igc-acid-base: output already exists: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".igc-acid-base-", dir=output.parent) as temporary:
        staged = Path(temporary)
        result.injections.to_csv(staged / "acid_base_injections.csv", index=False)
        result.interpolated.to_csv(staged / "acid_base_interpolated_vn.csv", index=False)
        result.schultz_lines.to_csv(staged / "acid_base_schultz_lines.csv", index=False)
        result.delta_g_sp.to_csv(staged / "acid_base_delta_g_sp.csv", index=False)
        result.profile.to_csv(staged / "acid_base_profile.csv", index=False)
        _write_figures(result, staged)
        (staged / "acid_base_run.json").write_text(
            json.dumps(record, indent=2, allow_nan=False) + "\n", encoding="utf-8"
        )
        _write_readme(result, staged)
        staged.replace(output)
    print(f"Acid/base outputs written to {output}")
    print(f"Reportable profile: {'yes' if result.reportable else 'no'}")
    print(result.qc["summary"])


if __name__ == "__main__":
    main()
