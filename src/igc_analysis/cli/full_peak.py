"""CLI for full-peak nonlinear inverse chromatography.

Fits a forward equilibrium-dispersive column model to the *complete*
baseline-corrected chromatograms of every probe injection, jointly across
blocks, and identifies an adsorption isotherm from peak shape rather than from
a single retention statistic.

Usage::

    igc-full-peak \\
      --neutral-bundle block1=/path/to/neutral_block1 \\
      --neutral-bundle block2=/path/to/neutral_block2 \\
      -o output/full_peak/

Outputs (all written to ``-o``): trace, transport, fit-summary, parameter,
prediction, residual and isotherm CSVs; a machine-readable JSON run record;
observed-vs-predicted, residual, isotherm and model-comparison figures as
vector PDF plus high-resolution PNG; and a Markdown README interpreting the
result.
"""

from __future__ import annotations

import argparse
import json
import platform
import re
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Keep PDF text as editable TrueType rather than outlines.
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42

from igc_analysis.analysis.full_peak import (
    CONDITION_LIMIT, CORRELATION_LIMIT, DEFAULT_N_CELLS,
    IDENTIFIABILITY_RSE_LIMIT, bracket_assignment_to_dataframe,
    build_trace_dataset_from_neutral, compare_models,
    compute_ssa_if_identifiable, predict_injection, recovered_isotherm,
    traces_to_dataframe, transport_sensitivity, transport_to_dataframe,
)
from igc_analysis.analysis.isotherm_models import MODELS, get_model
from igc_analysis.analysis.column_model import peak_moments

PNG_DPI = 300
NEUTRAL_LABEL_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$")


def _json_safe(value):
    """Convert NumPy scalars and non-finite floats to strict JSON values."""

    if isinstance(value, dict):
        return {key: _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def _save(fig, out_dir: Path, stem: str) -> None:
    """Save a figure as vector PDF and high-resolution PNG."""
    fig.savefig(out_dir / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(out_dir / f"{stem}.png", dpi=PNG_DPI, bbox_inches="tight")
    plt.close(fig)


def _plot_observed_vs_predicted(blocks, fits, best_name, out_dir, n_cells):
    model = get_model(best_name)
    fit = fits[best_name]
    n_inj = sum(len(b.injections) for b in blocks)
    ncols = 5
    nrows = int(np.ceil(n_inj / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.2 * ncols, 2.6 * nrows),
                             squeeze=False)
    k = 0
    for blk in blocks:
        for inj in blk.injections:
            ax = axes[k // ncols][k % ncols]
            res = predict_injection(inj, blk, model, fit.params, n_cells)
            ax.plot(inj.time_min, inj.c_out_mol_m3, lw=1.0, color="0.25",
                    label="observed")
            ax.plot(inj.time_min, res.c_out, lw=1.2, color="crimson",
                    ls="--", label="predicted")
            ax.set_title(f"{inj.block} inj{inj.injection_number}\n"
                         f"n={inj.n_injected_mol:.2e} mol", fontsize=8)
            ax.set_xlim(0, min(2.5, inj.time_min[-1]))
            ax.tick_params(labelsize=7)
            if k == 0:
                ax.legend(fontsize=7)
            k += 1
    for j in range(k, nrows * ncols):
        axes[j // ncols][j % ncols].axis("off")
    fig.suptitle(f"Observed vs predicted chromatograms — {best_name} model "
                 f"(normalised RMSE {fit.rmse_normalised:.4f})", fontsize=11)
    fig.supxlabel("time (min)", fontsize=9)
    fig.supylabel("outlet concentration (mol/m³)", fontsize=9)
    fig.tight_layout()
    _save(fig, out_dir, "full_peak_observed_vs_predicted")


def _plot_residuals(blocks, fits, best_name, out_dir, n_cells):
    model = get_model(best_name)
    fit = fits[best_name]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))

    ax = axes[0]
    for blk in blocks:
        for inj in blk.injections:
            res = predict_injection(inj, blk, model, fit.params, n_cells)
            ax.plot(inj.time_min, (inj.c_out_mol_m3 - res.c_out) / inj.peak_scale,
                    lw=0.8, alpha=0.8,
                    label=f"{inj.block} i{inj.injection_number}")
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xlim(0, 2.5)
    ax.set_xlabel("time (min)")
    ax.set_ylabel("normalised residual")
    ax.set_title(f"Residual traces — {best_name}")
    ax.legend(fontsize=6, ncol=2)

    ax = axes[1]
    names, vals = zip(*sorted(fit.per_injection_rmse.items()))
    ax.bar(range(len(vals)), vals, color="steelblue")
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, rotation=90, fontsize=7)
    ax.set_ylabel("normalised RMSE")
    ax.set_title("Per-injection fit quality")

    ax = axes[2]
    doses, shifts_obs, shifts_pred = [], [], []
    for blk in blocks:
        for inj in blk.injections:
            res = predict_injection(inj, blk, model, fit.params, n_cells)
            doses.append(inj.n_injected_mol)
            t0 = (inj.assigned_t0_min if inj.assigned_t0_min is not None
                  else blk.transport.t0_min)
            shifts_obs.append(peak_moments(inj.time_min, inj.c_out_mol_m3)[1]
                              - t0)
            shifts_pred.append(peak_moments(inj.time_min, res.c_out)[1]
                               - t0)
    ax.plot(doses, shifts_obs, "o", color="0.25", label="observed")
    ax.plot(doses, shifts_pred, "x", color="crimson", label="predicted")
    ax.set_xlabel("injected moles")
    ax.set_ylabel("net first moment  μ₁ − t₀ (min)")
    ax.set_title("Dose dependence of retention")
    ax.legend(fontsize=8)

    fig.tight_layout()
    _save(fig, out_dir, "full_peak_residual_diagnostics")


def _plot_bracket_assignments(assignments, out_dir):
    """Show methane-bracket drift and corrected retention versus dose/order."""
    blocks = list(assignments["block"].drop_duplicates())
    fig, axes = plt.subplots(len(blocks), 2, figsize=(10, 3.5 * len(blocks)),
                             squeeze=False)
    for row, block in enumerate(blocks):
        d = assignments[assignments["block"] == block].copy()
        order = np.arange(1, len(d) + 1)
        ax = axes[row, 0]
        ax.plot(order, d["assigned_t0_min"] * 60, "o-", color="tab:blue",
                label="assigned methane t₀")
        pre_t0 = pd.to_numeric(d["pre_bracket_t0_min"], errors="coerce").iloc[0]
        post_t0 = pd.to_numeric(d["post_bracket_t0_min"], errors="coerce").iloc[0]
        if np.isfinite(pre_t0):
            ax.axhline(pre_t0 * 60, color="0.4", ls="--",
                       label="pre-pair mean")
        if np.isfinite(post_t0):
            ax.axhline(post_t0 * 60, color="0.6", ls=":",
                       label="post-pair mean")
        ax.set_xlabel("probe injection order")
        ax.set_ylabel("t₀ (s)")
        ax.set_title(f"{block}: methane bracket assignment")
        ax.legend(fontsize=8)

        ax = axes[row, 1]
        ax.plot(d["n_injected_mol"], d["net_first_moment_min"] * 60,
                "o-", color="tab:orange")
        for i, rec in enumerate(d.itertuples(index=False), start=1):
            ax.annotate(str(i), (rec.n_injected_mol,
                                rec.net_first_moment_min * 60),
                        xytext=(4, 3), textcoords="offset points", fontsize=8)
        ax.set_xlabel("injected probe (mol)")
        ax.set_ylabel("probe first moment − assigned t₀ (s)")
        ax.set_title(f"{block}: drift-corrected retention")
    fig.tight_layout()
    _save(fig, out_dir, "full_peak_methane_bracket_diagnostics")


def _plot_isotherm(iso_df, fits, best_name, out_dir):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    axes[0].plot(iso_df["c_mol_m3"], iso_df["q_mol_g"], color="darkgreen", lw=1.6)
    axes[0].set_xlabel("gas-phase concentration c (mol/m³)")
    axes[0].set_ylabel("q (mol/g)")
    axes[0].set_title(f"Recovered isotherm q(c) — {best_name}")
    axes[1].plot(iso_df["c_mol_m3"], iso_df["dqdc_m3_g"], color="darkorange", lw=1.6)
    axes[1].set_xlabel("gas-phase concentration c (mol/m³)")
    axes[1].set_ylabel("dq/dc (m³/g)")
    axes[1].set_title("Isotherm slope dq/dc  (rising ⇒ cooperative)")
    for ax in axes:
        ax.grid(alpha=0.3)
    fig.tight_layout()
    _save(fig, out_dir, "full_peak_recovered_isotherm")


def _plot_model_comparison(table, out_dir):
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.0))
    m = table["model"]
    axes[0].bar(m, table["rmse_normalised"], color="steelblue")
    axes[0].set_ylabel("normalised RMSE")
    axes[0].set_title("In-sample fit")
    axes[1].bar(m, table["lodo_rmse"], color="indianred")
    axes[1].set_ylabel("held-out RMSE")
    axes[1].set_title("Leave-one-dose-out prediction")
    finite = table["aicc"].replace([np.inf, -np.inf], np.nan).dropna()
    base = finite.min() if len(finite) else 0.0
    axes[2].bar(m, table["aicc"] - base, color="seagreen")
    axes[2].set_ylabel("ΔAICc (vs best)")
    axes[2].set_title("AICc (autocorrelated residuals — guide only)")
    for ax in axes:
        ax.tick_params(axis="x", rotation=20)
        ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    _save(fig, out_dir, "full_peak_model_comparison")


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

def _markdown_table(df: pd.DataFrame) -> str:
    """Render a DataFrame as a Markdown table without needing `tabulate`."""
    def fmt(v):
        if isinstance(v, float):
            if not np.isfinite(v):
                return "n/a"
            return f"{v:.4g}"
        return str(v)

    header = "| " + " | ".join(str(c) for c in df.columns) + " |"
    sep = "|" + "|".join("---" for _ in df.columns) + "|"
    rows = ["| " + " | ".join(fmt(v) for v in row) + " |"
            for row in df.itertuples(index=False, name=None)]
    return "\n".join([header, sep] + rows)


def _scientific_interpretation_lines(fit, table, best_name) -> list[str]:
    """Build conservative, model-aware interpretation statements.

    Fit statistics compare candidate descriptions; they do not, by themselves,
    establish a molecular mechanism or diagnose a specific missing process.
    """
    lines: list[str] = []
    t = table.set_index("model")

    if "none" in t.index:
        r2_none = float(t.loc["none", "r_squared"])
        lines.append(
            f"- **Transport-only comparison.** The methane-calibrated, "
            f"no-adsorption control reaches R² = {r2_none:.3f}. Compare its "
            "held-out and residual diagnostics with the adsorption models; "
            "this statistic alone does not establish whether adsorption is "
            "required."
        )

    ranked = table.dropna(subset=["lodo_rmse"]).sort_values("lodo_rmse")
    if len(ranked) >= 2:
        b0, b1 = ranked.iloc[0], ranked.iloc[1]
        denominator = abs(float(b1["lodo_rmse"]))
        margin = (100.0 * (float(b1["lodo_rmse"]) - float(b0["lodo_rmse"]))
                  / denominator if denominator > 0 else float("nan"))
        margin_text = f"{margin:.1f}%" if np.isfinite(margin) else "an undefined margin"
        lines.append(
            f"- **Candidate-model ranking.** `{b0['model']}` gives the lowest "
            f"leave-one-dose-out error, by {margin_text} over "
            f"`{b1['model']}`. This is comparative evidence within the tested "
            "model set, not decisive mechanistic discrimination."
        )
        if np.isfinite(margin) and margin < 2.0:
            lines.append(
                f"  **Interpretive caution:** `{b0['model']}` and "
                f"`{b1['model']}` are not meaningfully separated by this "
                "dataset; parameters unique to either model should not be "
                "treated as demonstrated mechanisms."
            )

    if ("langmuir" in t.index
            and not bool(t.loc["langmuir", "all_params_identifiable"])):
        lines.append(
            "- **Langmuir identifiability.** At least one Langmuir parameter is "
            "not separately determined by these data. That diagnostic does not "
            "by itself establish linear-limit behavior, parameter correlation, "
            "or whether the measured range approached saturation."
        )

    if fit.cooperative:
        exponent = (float(fit.params[1]) if len(fit.params) > 1
                    else float("nan"))
        lines.append(
            "- **Concentration dependence.** Within the numerical best model, "
            f"the fitted exponent is {exponent:.3f} and the fitted `dq/dc` "
            "increases with concentration. This direction is model-conditional "
            "and must be assessed against dose order, methane-bracket "
            "sensitivity, residual structure, and held-out error."
        )

    best_r2 = float(t.loc[best_name, "r_squared"])
    if np.isfinite(best_r2):
        lines.append(
            f"- **Fit scope.** The numerical best model has R² = {best_r2:.3f}. "
            "Residual structure should be inspected directly; it is not "
            "mechanistically diagnostic of a particular kinetic, transport, "
            "or inlet-profile process."
        )
    else:
        lines.append(
            "- **Fit scope.** R² is unavailable for the numerical best model. "
            "Inspect the residual and held-out diagnostics before interpreting "
            "the recovered isotherm."
        )
    return lines


def _write_readme(out_dir, blocks, fits, table, best_name, ssa_verdict,
                  settings, sens=None) -> None:
    fit = fits[best_name]
    model = get_model(best_name)
    lines = []
    a = lines.append
    a("# Full-peak nonlinear inverse chromatography — results\n")
    a(f"Generated {settings['generated_utc']} with "
      f"`igc-analysis-pipeline` {settings['package_version']} "
      f"(commit-independent module `igc_analysis.analysis.full_peak`).\n")

    a("## What this is\n")
    mode = blocks[0].injections[0].transport_mode
    a("Every probe injection's **complete baseline-corrected trace** was fitted "
      "with a forward equilibrium-dispersive column model. Axial dispersion "
      "and inlet width were calibrated on the block's methane markers. "
      f"Dead time used `{mode}` assignment. See "
      "`docs/full_peak_architecture.md` for the equations.\n")
    if mode == "bracket_interpolated":
        a("The neutral bundle does not provide usable acquisition timestamps. "
          "Each probe t₀ was therefore interpolated by declared sequence position "
          "between the mean positions and mean first moments of the pre-block and "
          "post-block dead-time injections. This assumption is made visible in "
          "`methane_bracket_assignments.csv`.\n")

    a("## Dataset\n")
    for b in blocks:
        t0_sd_text = (f"{b.transport.t0_sd_min*60:.2f} s"
                      if np.isfinite(b.transport.t0_sd_min) else "unavailable")
        marker_word = "marker" if b.transport.n_markers == 1 else "markers"
        a(f"- **{b.block}** — {len(b.injections)} probe injections, "
          f"t₀ = {b.transport.t0_min:.4f} min "
          f"(SD {t0_sd_text} over "
          f"{b.transport.n_markers} methane {marker_word}), "
          f"N_eff = {b.transport.plate_number:.0f}, "
          f"V_void = {b.transport.void_volume_m3*1e6:.3f} mL.")
    pp0 = [float(np.max(i.pp0)) for b in blocks for i in b.injections]
    a(f"\nMeasured peak P/P₀ envelope: **{min(pp0):.4f} – {max(pp0):.4f}**.\n")

    a("## Model comparison\n")
    cols = ["model", "n_params", "rmse_normalised", "r_squared", "aicc",
            "lodo_rmse", "all_params_identifiable", "cooperative"]
    a(_markdown_table(table[cols]))
    a("")
    a("`lodo_rmse` is leave-one-dose-out held-out prediction error — the "
      "honest test, since in-sample RMSE always improves with more parameters. "
      "AICc is reported as a guide only: chromatographic residuals are strongly "
      "autocorrelated, so its absolute values overstate evidence.\n")

    a(f"## Numerically lowest-error model: `{best_name}`\n")
    a(f"{model.description}\n")
    if fit.n_params:
        a("| parameter | value | unit | conditional SE | rel. SE | identifiable |")
        a("|---|---|---|---|---|---|")
        for nm, v, u, se, rse in zip(fit.param_names, fit.params,
                                     fit.param_units, fit.std_errors,
                                     fit.rel_std_errors):
            a(f"| {nm} | {v:.5g} | {u} | {se:.3g} | "
              f"{'∞' if not np.isfinite(rse) else f'{rse:.1%}'} | "
              f"{'yes' if fit.identifiable[nm] else '**no**'} |")
        a("")
        a(f"Jacobian condition number {fit.condition_number:.3g}; maximum "
          f"parameter correlation {fit.max_abs_correlation:.4f} "
          f"(limits: {CONDITION_LIMIT:.0g} and {CORRELATION_LIMIT}). "
          f"A parameter is called identifiable only if its relative standard "
          f"error is below {IDENTIFIABILITY_RSE_LIMIT:.0%} *and* the fit is not "
          f"structurally confounded.\n")
        a("These Jacobian errors are conditional optimiser diagnostics, not "
          "full scientific uncertainty; chromatographic residuals are "
          "autocorrelated and transport/bracket sensitivity is reported "
          "separately.\n")
    if not fit.all_identifiable:
        a("> **Identifiability warning.** At least one parameter above is not "
          "separately determined by these data. Its numerical value is what the "
          "optimiser returned, not a measurement.\n")

    if sens is not None and not sens.empty:
        a("## Sensitivity to transport and methane-bracket assumptions\n")
        a("The isotherm parameters are only meaningful if they survive "
          "plausible changes in methane dead-time assignment, dispersion, and "
          "inlet width. Each row refits the numerical best model under one "
          "perturbation.\n")
        keep = [c for c in sens.columns
                if c == "scenario" or c.endswith("_pct_change")
                or c == "rmse_normalised"]
        a(_markdown_table(sens[keep]))
        a("")

    # --- Data-driven interpretation -------------------------------------
    a("## Scientific interpretation\n")
    for line in _scientific_interpretation_lines(fit, table, best_name):
        a(line)
    a("")

    a("## Surface area\n")
    if ssa_verdict.reportable:
        a(f"SSA = **{ssa_verdict.ssa_m2_g:.4f} m²/g** — {ssa_verdict.reason}\n")
    else:
        a(f"**No SSA is reported.** {ssa_verdict.reason}\n")
        a("The measured/predicted P/P₀ ratio is never converted into a surface "
          "area. Where no finite monolayer capacity is identifiable, the "
          "adsorption parameters, `q(c)` and `dq/dc` are the result and the "
          "geometric surface area remains undetermined.\n")

    a("## Limitations carried forward\n")
    for lim in settings["limitations"]:
        a(f"- {lim}")
    a("")
    a("## Files\n")
    a("`full_peak_traces.csv` (calibrated full traces), "
      "`methane_transport_summary.csv`, `full_peak_fit_summary.csv`, "
      "`methane_bracket_assignments.csv`, "
      "`full_peak_parameters.csv`, `full_peak_predictions.csv`, "
      "`full_peak_residuals.csv`, `recovered_isotherm.csv`, "
      "`full_peak_sensitivity.csv`, `full_peak_run.json`; figures as vector "
      "`.pdf` + 300 dpi `.png`.\n")

    (out_dir / "README.md").write_text("\n".join(lines))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def _parse_neutral_bundles(args) -> dict[str, str]:
    bundles: dict[str, str] = {}
    for spec in args.neutral_bundle or []:
        if "=" not in spec:
            raise SystemExit(f"--neutral-bundle expects LABEL=PATH, got {spec!r}")
        label, path = spec.split("=", 1)
        label = label.strip()
        if not NEUTRAL_LABEL_PATTERN.fullmatch(label):
            raise SystemExit(
                "--neutral-bundle label must be a neutral opaque identifier "
                "containing only letters, numbers, periods, underscores, or hyphens"
            )
        bundles[label] = path.strip()
    return bundles


def main(argv=None):
    p = argparse.ArgumentParser(
        description="Full-peak nonlinear inverse gas chromatography")
    p.add_argument(
        "--neutral-bundle",
        action="append",
        help="LABEL=PATH of an igc-neutral-data/0.2.0 bundle (repeatable)",
    )
    p.add_argument("--output", "-o", default="full_peak_results")
    p.add_argument("--probe", default="auto",
                   help="Probe-name override; default uses the bundle's sole analyte")
    p.add_argument(
        "--transport-mode", default="fixed_block_mean",
        choices=("fixed_block_mean", "bracket_interpolated", "bracket_pre",
                 "bracket_post"),
        help="Dead-time assignment. bracket_interpolated uses declared neutral "
             "sequence order when acquisition timestamps are absent.")
    p.add_argument("--models", default="none,henry,langmuir,freundlich",
                   help="Comma-separated isotherm models to compare")
    p.add_argument("--n-cells", type=int, default=DEFAULT_N_CELLS,
                   help=f"Column discretisation (default {DEFAULT_N_CELLS}); "
                        f"must be the same for calibration and fitting")
    p.add_argument("--n-starts", type=int, default=4,
                   help="Multistart count for the optimiser")
    p.add_argument("--no-lodo", action="store_true",
                   help="Skip leave-one-dose-out cross-validation (faster)")
    p.add_argument(
        "--lodo-models", default=None,
        help="Optional comma-separated subset for leave-one-dose-out, e.g. "
             "henry,freundlich. All models are still fitted in-sample.")
    p.add_argument("--cross-section", type=float, default=None,
                   help="Probe molecular cross-section (m²) for the gated SSA; "
                        "defaults to the declared probe property")
    args = p.parse_args(argv)

    neutral_bundles = _parse_neutral_bundles(args)
    if not neutral_bundles:
        p.error("no inputs given — use --neutral-bundle LABEL=PATH")

    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)
    model_names = tuple(m.strip() for m in args.models.split(",") if m.strip())
    for m in model_names:
        if m not in MODELS:
            p.error(f"unknown model {m!r}; available: {sorted(MODELS)}")
    lodo_models = (set(x.strip() for x in args.lodo_models.split(",") if x.strip())
                   if args.lodo_models else None)
    if lodo_models and not lodo_models <= set(model_names):
        p.error("--lodo-models must be a subset of --models")

    print(f"Full-peak inverse chromatography → {out_dir}")
    selected_inputs = neutral_bundles
    print(f"  blocks: {', '.join(neutral_bundles)}")
    print(f"  models: {', '.join(model_names)}  n_cells={args.n_cells}  "
          f"n_starts={args.n_starts}  LODO={'off' if args.no_lodo else 'on'}\n")

    blocks = build_trace_dataset_from_neutral(
        neutral_bundles,
        probe_override=args.probe,
        n_cells=args.n_cells,
        transport_mode=args.transport_mode,
        verbose=True,
    )

    # --- Tables ---
    traces = traces_to_dataframe(blocks)
    traces.to_csv(out_dir / "full_peak_traces.csv", index=False)
    transport = transport_to_dataframe(blocks)
    transport.to_csv(out_dir / "methane_transport_summary.csv", index=False)
    assignments = bracket_assignment_to_dataframe(blocks)
    assignments.to_csv(out_dir / "methane_bracket_assignments.csv", index=False)

    print("\n  Fitting models ...")
    fits, table = compare_models(blocks, model_names, n_cells=args.n_cells,
                                 n_starts=args.n_starts,
                                 do_lodo=not args.no_lodo,
                                 lodo_models=lodo_models, verbose=True)
    table.to_csv(out_dir / "full_peak_fit_summary.csv", index=False)

    # Select on held-out prediction when available, else AICc.
    key = "lodo_rmse" if table["lodo_rmse"].notna().any() else "aicc"
    best_name = table.sort_values(key).iloc[0]["model"]
    fit = fits[best_name]

    # Parameters table (every model, so nothing is hidden).
    prows = []
    for name, f in fits.items():
        for nm, v, u, se, rse in zip(f.param_names, f.params, f.param_units,
                                     f.std_errors, f.rel_std_errors):
            prows.append({"model": name, "parameter": nm, "value": v,
                          "unit": u, "std_error": se, "rel_std_error": rse,
                          "identifiable": f.identifiable[nm],
                          "condition_number": f.condition_number,
                          "max_abs_correlation": f.max_abs_correlation})
    pd.DataFrame(prows).to_csv(out_dir / "full_peak_parameters.csv", index=False)

    # Predictions and residuals for the selected model.
    model = get_model(best_name)
    pred_rows, resid_rows = [], []
    for blk in blocks:
        for inj in blk.injections:
            res = predict_injection(inj, blk, model, fit.params, args.n_cells)
            n = len(inj.time_min)
            pred_rows.append(pd.DataFrame({
                "block": inj.block, "injection": inj.injection_number,
                "time_min": inj.time_min,
                "c_observed_mol_m3": inj.c_out_mol_m3,
                "c_predicted_mol_m3": res.c_out,
                "model": best_name, "mass_balance": res.mass_balance}))
            resid_rows.append(pd.DataFrame({
                "block": inj.block, "injection": inj.injection_number,
                "time_min": inj.time_min,
                "residual_mol_m3": inj.c_out_mol_m3 - res.c_out,
                "residual_normalised":
                    (inj.c_out_mol_m3 - res.c_out) / inj.peak_scale,
                "model": best_name}))
    pd.concat(pred_rows, ignore_index=True).to_csv(
        out_dir / "full_peak_predictions.csv", index=False)
    pd.concat(resid_rows, ignore_index=True).to_csv(
        out_dir / "full_peak_residuals.csv", index=False)

    iso = recovered_isotherm(fit, blocks)
    iso.to_csv(out_dir / "recovered_isotherm.csv", index=False)

    # Sensitivity of the isotherm to the *fixed* transport assumptions.
    print("\n  Transport / window sensitivity ...")
    sens = transport_sensitivity(blocks, best_name, fit.params,
                                 n_cells=args.n_cells)
    sens.to_csv(out_dir / "full_peak_sensitivity.csv", index=False)

    # --- SSA guardrail ---
    if args.cross_section is not None:
        a_cross = args.cross_section
    else:
        from igc_analysis.io.neutral_data import read_neutral_bundle

        neutral = read_neutral_bundle(next(iter(neutral_bundles.values())))
        properties = neutral.table("probe_properties.csv")
        probe = blocks[0].injections[0].probe.casefold()
        matches = properties[
            properties["probe_name"].astype(str).str.casefold() == probe
        ]
        if len(matches) != 1 or str(matches.iloc[0]["cross_section_m2"]) == "":
            p.error("selected neutral probe has no unique declared cross_section_m2")
        a_cross = float(matches.iloc[0]["cross_section_m2"])
    ssa_verdict = compute_ssa_if_identifiable(fit, blocks, a_cross)

    # --- Figures ---
    print("\n  Writing figures ...")
    _plot_observed_vs_predicted(blocks, fits, best_name, out_dir, args.n_cells)
    _plot_residuals(blocks, fits, best_name, out_dir, args.n_cells)
    _plot_isotherm(iso, fits, best_name, out_dir)
    _plot_model_comparison(table, out_dir)
    _plot_bracket_assignments(assignments, out_dir)

    # --- Machine-readable run record ---
    from igc_analysis import __version__ as pkg_version
    input_mode = "igc-neutral-data/0.2.0"
    limitations = [
        "Each input mapping is treated as one independently characterised block.",
        "Methane and liquid-probe inlet profiles may differ.",
        "Packed-bed length is not required; the model uses dead-time-derived "
        "void time/volume and an effective plate number instead.",
        "Scientific interpretation still requires a declared study design, "
        "replication structure, and reportability review.",
    ]
    if any(not injection.source_methane_pre or not injection.source_methane_post
           for block in blocks for injection in block.injections):
        limitations.append(
            "At least one block lacks dead-time markers on both sides of its probe block."
        )
    if args.transport_mode == "bracket_interpolated":
        limitations.append(
            "Bracket interpolation uses neutral sequence position when trustworthy "
            "per-injection timestamps are absent."
        )

    neutral_input_provenance = {}
    from igc_analysis.io.neutral_data import read_neutral_bundle

    for label, path in neutral_bundles.items():
        neutral = read_neutral_bundle(path)
        neutral_input_provenance[label] = {
            "dataset_id": neutral.dataset_id,
            "contract_version": neutral.contract_version,
            "source_fingerprint": neutral.manifest.get("source_fingerprint"),
        }

    settings = {
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "package_version": pkg_version,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "command_name": "igc-full-peak",
        "input_mode": input_mode,
        "input_labels": list(selected_inputs),
        "neutral_inputs": neutral_input_provenance,
        "models": list(model_names),
        "n_cells": args.n_cells,
        "n_starts": args.n_starts,
        "lodo": not args.no_lodo,
        "lodo_models": sorted(lodo_models) if lodo_models else "all",
        "probe_override": args.probe,
        "transport_mode": args.transport_mode,
        "cross_section_m2": a_cross,
        "identifiability_rse_limit": IDENTIFIABILITY_RSE_LIMIT,
        "correlation_limit": CORRELATION_LIMIT,
        "condition_limit": CONDITION_LIMIT,
        "limitations": limitations,
    }
    record = {
        "settings": settings,
        "selection_criterion": key,
        "selected_model": best_name,
        "transport": transport.to_dict(orient="records"),
        "model_comparison": table.to_dict(orient="records"),
        "parameter_bounds": {n: {"names": get_model(n).param_names,
                                 "lower": list(get_model(n).bounds[0]),
                                 "upper": list(get_model(n).bounds[1])}
                             for n in model_names},
        "convergence": {n: {"converged": bool(f.converged),
                            "n_starts": f.n_starts,
                            "n_successful_starts": f.n_successful_starts,
                            "message": f.message}
                        for n, f in fits.items()},
        "identifiability": {n: f.identifiable for n, f in fits.items()},
        "qc": {
            "mass_balance_mean": float(fit.mass_balance_mean),
            "mass_balance_min": float(fit.mass_balance_min),
            "all_params_identifiable": bool(fit.all_identifiable),
            "cooperative": bool(fit.cooperative),
            "n_injections": int(sum(len(b.injections) for b in blocks)),
            "blocks_used": [b.block for b in blocks],
        },
        "ssa": {"reportable": ssa_verdict.reportable,
                "ssa_m2_g": ssa_verdict.ssa_m2_g,
                "reason": ssa_verdict.reason},
    }
    (out_dir / "full_peak_run.json").write_text(
        json.dumps(_json_safe(record), indent=2, allow_nan=False)
    )

    _write_readme(out_dir, blocks, fits, table, best_name, ssa_verdict,
                  settings, sens=sens)

    # --- Console summary ---
    print(f"\n{'=' * 72}")
    print(table[["model", "n_params", "rmse_normalised", "r_squared",
                 "lodo_rmse", "all_params_identifiable",
                 "cooperative"]].to_string(index=False))
    print(f"\nSelected model (by {key}): {best_name}")
    for nm, v, u, rse in zip(fit.param_names, fit.params, fit.param_units,
                             fit.rel_std_errors):
        flag = "" if fit.identifiable[nm] else "   [NOT IDENTIFIABLE]"
        rse_s = "inf" if not np.isfinite(rse) else f"{rse:.1%}"
        print(f"  {nm} = {v:.5g} {u}  (rel SE {rse_s}){flag}")
    print(f"\nSSA: {'REPORTED' if ssa_verdict.reportable else 'NOT REPORTED'} "
          f"— {ssa_verdict.reason}")
    print(f"\nOutputs written to {out_dir}")


if __name__ == "__main__":
    main()
