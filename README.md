# igc-analysis-pipeline

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21760976.svg)](https://doi.org/10.5281/zenodo.21760976)

A source-neutral Python toolkit for inverse gas chromatography analysis. Public
workflows accept documented tabular inputs and do not contain protected IGC
acquisition-file access, acquisition-format interpretation, credentials, or
third-party comparison adapters. A generic particle-size CSV reader remains
for the geometric surface-area dosing utility.

## Current release boundary

The implemented end-to-end chromatography workflows are corrected BET surface
area and full-peak nonlinear analysis through `igc-neutral-data/0.2.0`. The
repository also provides a geometric surface-area dosing utility.

Calculation modules for dispersive surface energy, acid-base analysis,
retention, peak detection, and quality control remain available for development
and unit testing. Their former source-coupled command-line workflows have been
removed. They will be exposed again only after neutral-contract consumers and
source-attributed property requirements are complete.

## Installation

Python 3.10 or newer is required.

```bash
git clone https://github.com/Brunaugh-Lab/igc-analysis-pipeline.git
cd igc-analysis-pipeline
python -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[dev]"
```

On Windows PowerShell, activate the environment with:

```powershell
.venv\Scripts\Activate.ps1
```

No acquisition-system runtime or source-specific library is required by this
repository.

Lab members should begin with the [student quick start](docs/student_quick_start.md),
which connects private extraction, neutral-bundle validation, analysis, QC
review, record retention, and citation without assuming a particular computer
or local folder layout.

## BET surface-area analysis

Run the corrected BET workflow on one validated neutral bundle:

```bash
igc-bet \
  --neutral-bundle /path/to/neutral_bundle \
  --output output/bet
```

The default calculation preserves the validated conventions from the final
pre-split workflow: eluted peak-apex concentration, calibrated injected amount,
per-injection measured column flow, direct absolute-pressure or pressure-drop
James--Martin correction, matched dead-time convention, sensitivity analyses,
and a strict reportability gate. The legacy loop-concentration convention is
available only as an explicit comparison with `--concentration loop`.

Outputs include injection, isotherm, and linearization CSVs; vector and raster
diagnostic figures; a strict-JSON run record; and a short interpretation README.
A numerical fit is not automatically reportable: use the reportability verdict
and review every QC message.

Run the packaged nondegenerate closed-form synthetic example:

```bash
igc-bet \
  --synthetic-example \
  --output output/synthetic-bet
```

This fixture is generated from a known Type II BET isotherm and exercises a
finite fit and positive reportability gate. It demonstrates synthetic numerical
recovery, not equivalence on governed experimental data.

See `docs/bet_architecture.md` for equations, readiness requirements, and the
verification boundary.

## Full-peak analysis

Each input must be a validated `igc-neutral-data/0.2.0` bundle. Supply one
opaque label per independently characterized acquisition block:

```bash
igc-full-peak \
  --neutral-bundle block-001=/path/to/neutral_block_001 \
  --neutral-bundle block-002=/path/to/neutral_block_002 \
  --transport-mode bracket_interpolated \
  --output output/full_peak
```

The command validates hashes, row counts, headers, units, controlled values,
ordering, keys, calibration consistency, and scientific provenance before any
analysis runs. It then:

1. preserves and baseline-corrects each complete detector trace;
2. applies the declared area-to-amount calibration;
3. converts the trace to outlet concentration while conserving injected mass;
4. characterizes dead-time transport independently within each block;
5. fits candidate adsorption models jointly across injections;
6. reports parameter identifiability, model comparison, residuals, sensitivity,
   and a gated surface-area verdict.

Outputs include analysis-ready CSV tables, vector and raster figures, a
machine-readable run record, and a Markdown interpretation summary. Local input
paths are not copied into the run record.

See `docs/full_peak_architecture.md` for the model and
`docs/neutral_data_contract.md` for the input boundary.

### Runnable synthetic example

From a source checkout, the bundled synthetic fixture exercises the complete
command without using experimental data:

```bash
igc-full-peak \
  --neutral-bundle synthetic=src/igc_analysis/contracts/igc-neutral-data/0.2.0/examples/synthetic_peak_shape \
  --models none,henry \
  --n-cells 40 \
  --n-starts 1 \
  --no-lodo \
  --output output/synthetic-smoke
```

The reduced model set, grid, and disabled cross-validation make this a quick
installation check. Remove those speed-oriented options for a scientific run
and review every generated diagnostic before interpreting the result.

## Neutral contract

The bundled contract lives at:

```text
src/igc_analysis/contracts/igc-neutral-data/0.2.0/
```

It includes:

- `schema.json` — logical schema and controlled values;
- `FIELD_DICTIONARY.md` — field meaning, units, and scientific boundaries;
- `validate_bundle.py` — dependency-free validator;
- `MIGRATION.md` — compatibility rules;
- `examples/synthetic_peak_shape/` — fully synthetic example bundle.

Validate a bundle directly with:

```bash
python src/igc_analysis/contracts/igc-neutral-data/0.2.0/validate_bundle.py \
  /path/to/neutral_bundle
```

Contract `0.1.0` was an experimental baseline and is not accepted by the
current reader. Inputs are never silently coerced between contract versions.

## Other commands

Geometric surface-area dosing consumes cumulative volume-distribution data:

```bash
igc-ssa-dose /path/to/distribution.csv --density 1.2
```

Use each command's `--help` output for its current options.

## Scientific guardrails

- Structural contract validity does not establish scientific reportability.
- Each bundle currently represents one experiment and one acquisition block.
- Multiple blocks are never pooled for dead-time characterization.
- Study relationships and replicate structure must be declared outside the
  core acquisition bundle until a versioned study-design contract exists.
- Peak asymmetry is not uniquely diagnostic of energetic heterogeneity;
  transport, packing, diffusion, baseline, dose nonlinearity, and detector
  behavior remain alternative explanations.
- A fitted adsorption model does not automatically justify a monolayer capacity
  or specific surface area. Those outputs are gated by model structure and
  parameter identifiability.
- Probe properties and calibration parameters must retain their declared
  provenance. Public analysis does not silently substitute source-specific
  defaults.

## Development

```bash
python -m pytest
ruff check .
python -m build
```

Tests use synthetic inputs unless a governed integration test is explicitly
enabled outside this repository. Real experimental bundles must not be added to
Git.

Before contributing, read `CONTRIBUTING.md`. The repository includes a CI
boundary check that rejects tracked files under the reserved local-data paths;
`.gitignore` is not the only protection against accidental recommitment.

GitHub Actions verifies the locked environment, runs the test suite on Python
3.10 through 3.14, enforces release-blocking lint checks, builds and installs
both distribution formats, runs the BET and full-peak wheels on synthetic data,
and audits the locked runtime dependencies for known vulnerabilities.

## Repository structure

```text
src/igc_analysis/
├── analysis/       Scientific calculations and QC
├── cli/            Supported command-line entry points
├── contracts/      Versioned neutral input contract and synthetic fixture
├── io/             Neutral-bundle and particle-size readers
└── plotting/       Reusable figure generation
tests/               Synthetic and calculation-level tests
docs/                Architecture and contract documentation
```

## Release status

The repository owner approved release under the MIT License and public
visibility after the content, history, package, and reproducibility checks in
the release checklist pass.

Scientific-reference provenance remains an ongoing documentation task rather
than permission to overinterpret results. The completed release evidence and
post-release checks are recorded in `docs/public_release_checklist.md`.

## Citation

If you use this software for an analysis, cite the archived version you used:

- version `v2026.8.4`: <https://doi.org/10.5281/zenodo.21796032>
- version `v2026.7.31`: <https://doi.org/10.5281/zenodo.21760977>
- all-versions concept DOI: <https://doi.org/10.5281/zenodo.21760976>

Use the version DOI in manuscripts that depend on a specific software release.
Use the concept DOI on websites or when referring to the evolving project
generally. Machine-readable citation metadata are provided in `CITATION.cff`.
Release `v2026.8.4` adds the corrected source-neutral BET workflow.

## License

Copyright (c) 2026 Ashlee Brunaugh. This project is distributed under the MIT
License; see `LICENSE` for the complete terms.
