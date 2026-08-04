"""Command-line dispersive surface-energy analysis for a neutral bundle."""

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
from igc_analysis.analysis.dispersive_workflow import (
    GAMMA_D_CRITICAL_BOUNDS_MJ_M2,
    GAMMA_D_WARNING_BOUNDS_MJ_M2,
    run_dispersive_from_neutral,
)
from igc_analysis.constants import R_GAS
from igc_analysis.io.neutral_data import (
    NeutralContractError,
    bundled_contract_path,
    read_neutral_bundle,
)


SYNTHETIC_EXAMPLE = (
    bundled_contract_path() / "examples" / "synthetic_dispersive_profile"
)


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
    profile = result.gamma_d.sort_values("coverage")
    fig, ax = plt.subplots(figsize=(6.2, 4.4))
    ax.plot(
        profile["coverage"], profile["gamma_d_mJm2"], "o-",
        label="center of mass (primary)", color="tab:blue",
    )
    ax.plot(
        profile["coverage"], profile["gamma_d_pm_mJm2"], "s--",
        label="peak maximum (sensitivity)", color="tab:orange",
    )
    ax.set_xlabel("fractional surface coverage")
    ax.set_ylabel("dispersive surface energy, gamma_d (mJ/m2)")
    ax.set_title("Coverage-resolved dispersive surface energy")
    ax.legend()
    ax.grid(alpha=0.25)
    fig.tight_layout()
    _save(fig, output, "dispersive_profile")

    cofm = result.interpolated[result.interpolated["retention_mode"] == "cofm"]
    coverages = sorted(cofm["target_coverage"].unique())
    ncols = 3
    nrows = int(np.ceil(len(coverages) / ncols))
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(4.2 * ncols, 3.2 * nrows), squeeze=False
    )
    for index, coverage in enumerate(coverages):
        ax = axes[index // ncols][index % ncols]
        data = cofm[cofm["target_coverage"] == coverage].dropna(
            subset=["VN_mL_g"]
        ).sort_values("carbon_number")
        x = data["carbon_number"].to_numpy(dtype=float)
        y = R_GAS * result.temperature_K * np.log(
            data["VN_mL_g"].to_numpy(dtype=float)
        ) / 1000.0
        ax.plot(x, y, "o", color="tab:blue")
        fits = profile[np.isclose(profile["coverage"], coverage)]
        if len(x) >= 2 and not fits.empty and np.isfinite(
            float(fits.iloc[0]["slope_Jmol"])
        ):
            fit = fits.iloc[0]
            grid = np.linspace(x.min(), x.max(), 100)
            ax.plot(
                grid,
                (fit["slope_Jmol"] * grid + fit["intercept"]) / 1000.0,
                "--",
                color="0.25",
            )
            ax.set_title(
                f"coverage={coverage:.3g}; gamma_d={fit['gamma_d_mJm2']:.2f}; "
                f"R2={fit['r_squared']:.5f}"
            )
        else:
            ax.set_title(f"coverage={coverage:.3g}; insufficient homologs")
        ax.set_xlabel("carbon number")
        ax.set_ylabel("RT ln(VN) (kJ/mol)")
        ax.grid(alpha=0.2)
    for index in range(len(coverages), nrows * ncols):
        axes[index // ncols][index % ncols].axis("off")
    fig.tight_layout()
    _save(fig, output, "dispersive_alkane_lines")


def _write_readme(result, record: dict, output: Path) -> None:
    profile = result.gamma_d.sort_values("coverage")
    gd_low = float(profile["gamma_d_mJm2"].min())
    gd_high = float(profile["gamma_d_mJm2"].max())
    lines = [
        "# Dispersive surface-energy result",
        "",
        f"- Dataset: `{result.dataset_id}`",
        f"- Numerical gamma_d range: {gd_low:.4g} to {gd_high:.4g} mJ/m2",
        f"- Reportable profile: {'yes' if result.reportable else 'no'}",
        f"- QC: {result.qc['summary']}",
        f"- Profile-shape label: `{result.qc['profile_shape']}`",
        f"- Supplied SSA: {result.specific_surface_area_m2_g:.6g} m2/g",
        f"- SSA provenance: {result.surface_area_source}",
        "",
        "Center-of-mass retention is the primary result. Peak maximum is",
        "reported as a sensitivity calculation because peak tailing can change",
        "the inferred profile. Review both figures and every QC flag in",
        "`dispersive_run.json` before reporting a value.",
        "",
        "The profile-shape label is descriptive, not a unique mechanistic",
        "assignment. Coverage dependence can reflect surface energetics,",
        "transport, packing, peak shape, calibration, or extrapolation.",
        "",
        "Conventions used:",
        "",
        "- actual coverage from calibrated amount and the declared SSA",
        "- piecewise-linear mapping from actual to target coverage",
        f"- extrapolation: {'enabled and flagged' if record['settings']['extrapolate'] else 'disabled'}",
        f"- pressure correction: {'applied' if record['settings']['pressure_correction'] else 'disabled'}",
        "- work of cohesion: two times gamma_d",
        "",
    ]
    (output / "README.md").write_text("\n".join(lines), encoding="utf-8")


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Dorris-Gray dispersive surface energy from "
            "igc-neutral-data/0.2.0"
        )
    )
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument("--neutral-bundle")
    inputs.add_argument(
        "--synthetic-example",
        action="store_true",
        help="run the packaged closed-form Dorris-Gray verification fixture",
    )
    parser.add_argument(
        "--output", "-o", default="dispersive_results",
        help="new output directory; the command refuses an existing path",
    )
    parser.add_argument("--ambient-pressure-pa", type=float, default=101325.0)
    parser.add_argument(
        "--homologous-probe-id",
        action="append",
        help=(
            "opaque probe ID in the homologous series; repeat at least three "
            "times when the bundle also contains other analytes"
        ),
    )
    parser.add_argument("--no-pressure-correction", action="store_true")
    parser.add_argument(
        "--no-extrapolation",
        action="store_true",
        help="leave probe/coverage combinations outside measured ranges undefined",
    )
    parser.add_argument("--max-temperature-span-k", type=float, default=1.0)
    parser.add_argument("--max-flow-relative-span", type=float, default=0.05)
    args = parser.parse_args(argv)
    bundle_path = SYNTHETIC_EXAMPLE if args.synthetic_example else args.neutral_bundle

    try:
        neutral = read_neutral_bundle(bundle_path)
        result = run_dispersive_from_neutral(
            neutral,
            homologous_probe_ids=args.homologous_probe_id,
            pressure_correction=not args.no_pressure_correction,
            ambient_pressure_pa=args.ambient_pressure_pa,
            extrapolate=not args.no_extrapolation,
            max_temperature_span_K=args.max_temperature_span_k,
            max_flow_relative_span=args.max_flow_relative_span,
        )
    except (NeutralContractError, ValueError) as exc:
        raise SystemExit(f"igc-dispersive: {exc}") from exc

    probe_properties = neutral.table("probe_properties.csv")
    record = {
        "settings": {
            "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "package_version": __version__,
            "python": platform.python_version(),
            "command_name": "igc-dispersive",
            "primary_retention_mode": "cofm",
            "sensitivity_retention_mode": "peak_max",
            "homologous_probe_selection": (
                "explicit" if args.homologous_probe_id else "all_carbon_numbered"
            ),
            "pressure_correction": not args.no_pressure_correction,
            "ambient_pressure_pa": args.ambient_pressure_pa,
            "extrapolate": not args.no_extrapolation,
            "max_temperature_span_K": args.max_temperature_span_k,
            "max_flow_relative_span": args.max_flow_relative_span,
            "gamma_d_warning_bounds_mJm2": GAMMA_D_WARNING_BOUNDS_MJ_M2,
            "gamma_d_critical_bounds_mJm2": GAMMA_D_CRITICAL_BOUNDS_MJ_M2,
            "reportability_rule": (
                "zero critical QC flags; at least three homologs with finite "
                "positive Dorris-Gray slope at every target coverage; no "
                "interpolated value outside its probe's measured coverage range; "
                "one detector gain across required injections"
            ),
        },
        "input": {
            "dataset_id": neutral.dataset_id,
            "manifest_fingerprint": _manifest_digest(neutral.manifest),
            "contract_name": neutral.manifest["contract_name"],
            "contract_version": neutral.manifest["contract_version"],
            "sample_id": result.sample_id,
            "specific_surface_area_m2_g": result.specific_surface_area_m2_g,
            "surface_area_source": result.surface_area_source,
            "flow_source_channels": result.flow_source_channels,
            "pressure_basis": result.pressure_basis,
            "pressure_roles": result.pressure_roles,
            "properties_sources": result.properties_sources,
            "calibration_sources": result.calibration_sources,
            "detector_gains": result.detector_gains,
            "probe_ids": sorted(result.injections["probe_id"].unique()),
            "carbon_numbers": sorted(
                int(value) for value in result.injections["carbon_number"].unique()
            ),
            "dead_time_injection_ids": neutral.table("injections.csv").loc[
                lambda table: table["role"] == "dead_time", "injection_id"
            ].astype(str).tolist(),
            "probe_property_records": probe_properties.loc[
                probe_properties["probe_id"].isin(result.injections["probe_id"]),
                ["probe_id", "probe_name", "cross_section_m2", "carbon_number",
                 "properties_source"],
            ].to_dict(orient="records"),
        },
        "result": {
            "reportable": result.reportable,
            "temperature_K": result.temperature_K,
            "dead_time_cofm_min": result.dead_time_cofm_min,
            "dead_time_peak_max_min": result.dead_time_peak_max_min,
            "coverage_points": len(result.gamma_d),
            "gamma_d_min_mJm2": float(result.gamma_d["gamma_d_mJm2"].min()),
            "gamma_d_max_mJm2": float(result.gamma_d["gamma_d_mJm2"].max()),
            "qc": result.qc,
        },
    }
    record = _json_safe(record)
    output = Path(args.output)
    if output.exists():
        raise SystemExit(f"igc-dispersive: output already exists: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix=".igc-dispersive-", dir=output.parent
    ) as temporary:
        staged = Path(temporary)
        result.injections.to_csv(
            staged / "dispersive_injections.csv", index=False
        )
        result.interpolated.to_csv(
            staged / "dispersive_interpolated_vn.csv", index=False
        )
        result.gamma_d.to_csv(staged / "dispersive_profile.csv", index=False)
        _write_figures(result, staged)
        (staged / "dispersive_run.json").write_text(
            json.dumps(record, indent=2, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        _write_readme(result, record, staged)
        staged.replace(output)
    print(f"Dispersive outputs written to {output}")
    print(f"Reportable profile: {'yes' if result.reportable else 'no'}")
    print(result.qc["summary"])


if __name__ == "__main__":
    main()
