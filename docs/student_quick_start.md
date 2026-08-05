# Student quick start

This guide is the standard Brunaugh Lab workflow for running the public IGC
analysis pipeline. It deliberately does not document acquisition-file access.
Lab members who need to create neutral bundles must also have access to the
lab's private extraction repository.

## 1. Use a fixed software release

Do not change software versions partway through a study. Use the exact release
selected for the study. Release `v2026.8.6` supports corrected BET, dispersive
surface energy, acid/base characterization, and explicit batch execution from
validated neutral bundles. Its version DOI will be recorded in this guide after
Zenodo archives the release.

Install Python 3.10 or newer, create an isolated environment, and install the
pinned release:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install "igc-analysis-pipeline @ git+https://github.com/Brunaugh-Lab/igc-analysis-pipeline.git@v2026.8.6"
```

Confirm the installed version:

```bash
python -c "import igc_analysis; print(igc_analysis.__version__)"
```

## 2. Check the installation with synthetic data

This check uses no experimental data:

```bash
igc-bet --synthetic-example --output output/synthetic-bet
```

The command should finish successfully and produce a reportable synthetic BET
fit. This confirms that the software runs; it does not validate an experimental
dataset.

Also check dispersive analysis:

```bash
igc-dispersive --synthetic-example --output output/synthetic-dispersive
```

The resulting profile should be reportable and close to 40.0 mJ/m².

Also check acid/base analysis:

```bash
igc-acid-base --synthetic-example --output output/synthetic-acid-base
```

The resulting profile should be reportable, with Ka close to 0.03 and Kb close
to 0.05. These are synthetic installation checks, not experimental validation.

Finally, check the batch command:

```bash
igc-report --synthetic-example --output output/synthetic-batch
```

This runs packaged BET, dispersive, acid/base, and reduced full-peak examples
as four explicit jobs. It checks installation and orchestration, not a study
design or experimental reportability.

## 3. Create and validate the neutral bundle

Use the private lab extraction workflow to convert the acquisition data into an
`igc-neutral-data/0.2.0` bundle. Keep the raw acquisition data, neutral bundle,
and derived output outside both Git repositories.

Before analysis, validate the bundle with the validator supplied by the private
workflow or this repository. Do not manually edit a bundle to make validation
pass. Resolve the source data or metadata issue instead.

## 4. Run the selected analysis

For corrected BET surface area:

```bash
igc-bet \
  --neutral-bundle /path/to/neutral_bundle \
  --output /path/to/results/bet
```

For full-peak nonlinear analysis, follow the command and model guidance in the
[main README](../README.md#full-peak-analysis). Use one opaque label for each
independently characterized acquisition block.

For dispersive surface energy:

```bash
igc-dispersive \
  --neutral-bundle /path/to/neutral_bundle \
  --output /path/to/results/dispersive
```

If the command reports more than three carbon-numbered analytes or duplicate
carbon numbers, inspect `probe_properties.csv` and repeat
`--homologous-probe-id OPAQUE_ID` for each member of the intended homologous
series. Do not substitute names or guess the series; ask the study lead if the
IDs are not documented.

Review the actual-versus-target coverage table and every extrapolation flag.
Any extrapolation makes the profile non-reportable. Center-of-mass retention is
primary; peak maximum is a sensitivity result.
Treat a detector-gain warning as a calibration review requirement rather than
an instruction to rescale the output manually. The pipeline marks that profile
non-reportable until the review is resolved.

For acid/base characterization, copy the documented opaque homologous and
polar probe IDs from the study analysis plan and declare each explicitly:

```bash
igc-acid-base \
  --neutral-bundle /path/to/neutral_bundle \
  --homologous-probe-id OPAQUE_HOMOLOG_ID_1 \
  --homologous-probe-id OPAQUE_HOMOLOG_ID_2 \
  --homologous-probe-id OPAQUE_HOMOLOG_ID_3 \
  --polar-probe-id OPAQUE_POLAR_ID_1 \
  --polar-probe-id OPAQUE_POLAR_ID_2 \
  --polar-probe-id OPAQUE_POLAR_ID_3 \
  --output /path/to/results/acid-base
```

Do not choose probes by name or fill missing property values from memory. The
neutral bundle must contain the study-approved, source-attributed property set.
Ka and Kb are convention-dependent descriptors; peak maximum is a sensitivity
result, and van Oss components are not part of this workflow.

For several independent analyses, create an explicit
`igc-analysis-batch/0.1.0` manifest and run:

```bash
igc-report \
  --manifest /path/to/batch.json \
  --output /path/to/results/batch-001
```

Every job accepts exactly one neutral bundle. The command does not infer
replicates, pool acquisitions, create cross-bundle fits, or produce a combined
scientific verdict. Review the manifest and each child result separately. See
[`batch_reporting_architecture.md`](batch_reporting_architecture.md) for the
manifest schema and allowed settings.

## 5. Review before reporting

Opening the final CSV is not sufficient. Review:

1. the run record and every QC message;
2. the diagnostic figures and residuals;
3. the reportability verdict;
4. the input-bundle identity and software version; and
5. any sensitivity analysis relevant to the result.

A numerical BET fit, surface-area value, dispersive profile, or acid/base
profile must not be reported when the pipeline marks it non-reportable. Ask the
study lead when a warning, extrapolation, probe policy, or model choice is
unclear.

## 6. Preserve the analysis record

Keep the validated neutral bundle, complete output directory, run record, and
study-specific interpretation together in the study's governed data location.
Record the exact software release and version DOI in the analysis note or
manuscript. Do not commit experimental inputs or outputs to this repository.

For a manuscript, cite the version DOI listed in the [citation section of the
README](../README.md#citation). For a website or general reference to the
evolving project, use the concept DOI.
