"""CLI for converting cumulative Q3 data into D32 and SSA_dose."""

from __future__ import annotations

import argparse
from pathlib import Path

from igc_sea.analysis.geometric_ssa import analyze_distribution_path, summarize_replicates


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Calculate D(3,2) and optional geometric SSA_dose from cumulative Q3 data",
    )
    parser.add_argument("input", help="One Q3 CSV or a directory searched recursively")
    parser.add_argument("--output", "-o", default="geometric_ssa_results",
                        help="Output directory (default: geometric_ssa_results)")
    parser.add_argument("--density", type=float, default=None,
                        help="Particle density in g/cm3; omit to calculate D32 only")
    parser.add_argument("--density-basis", default="unspecified",
                        choices=["unspecified", "skeletal", "envelope"],
                        help="Provenance label for the supplied density")
    parser.add_argument("--first-bin-warning", type=float, default=10.0,
                        help="Flag when first populated bin exceeds this volume %% (default: 10)")
    args = parser.parse_args(argv)

    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    per_file = analyze_distribution_path(
        args.input,
        density_g_cm3=args.density,
        density_basis=args.density_basis,
        first_bin_warning_percent=args.first_bin_warning,
    )
    summary = summarize_replicates(per_file)
    per_file_path = output / "d32_per_measurement.csv"
    summary_path = output / "d32_by_formulation.csv"
    per_file.to_csv(per_file_path, index=False)
    summary.to_csv(summary_path, index=False)

    print(f"Processed {len(per_file)} distribution file(s) across "
          f"{per_file['formulation_id'].nunique(dropna=False)} formulation(s)")
    print(f"Per-measurement output: {per_file_path}")
    print(f"Formulation summary: {summary_path}")
    if args.density is None:
        print("Density not supplied: D32 calculated; SSA columns left blank.")
    flagged = int((per_file["qc_flags"] != "").sum())
    print(f"QC-flagged measurements: {flagged}/{len(per_file)}")
    failures = per_file.attrs.get("failures", [])
    if failures:
        print(f"Skipped incompatible or invalid CSV files: {len(failures)}")


if __name__ == "__main__":
    main()
