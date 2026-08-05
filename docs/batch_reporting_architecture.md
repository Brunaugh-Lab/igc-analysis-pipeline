# Source-neutral batch-reporting architecture

Status: verified for release `v2026.8.6`. Entry point: `igc-report`.
Independent audit, adversarial rollback/state-isolation tests, installed wheel
and source-distribution runs, and the complete release gate are passed.

## Purpose and boundary

`igc-report` runs an explicit list of supported chromatography workflows and
collects their existing outputs into one transactional directory. It accepts
only validated `igc-neutral-data/0.2.0` bundles. Inputs and jobs are declared
explicitly rather than discovered from directories or inferred from filenames.

The batch manifest is an execution plan, not a study-design model. Bundle and
job IDs are opaque bookkeeping identifiers. The orchestrator does not infer
replicates, treatment groups, before/after relationships, statistical units,
or permission to pool acquisitions. Those relationships require a separately
reviewed study design keyed by immutable neutral `dataset_id` values.

## Manifest 0.1.0

The machine-readable schema is packaged at
`contracts/igc-analysis-batch/0.1.0/schema.json`. Relative bundle paths are
resolved from the manifest directory, allowing an analysis folder to
move without rewriting every entry.

```json
{
  "schema_version": "igc-analysis-batch/0.1.0",
  "batch_id": "batch-001",
  "bundles": [
    {"bundle_id": "bundle-001", "path": "neutral/bundle-001"}
  ],
  "jobs": [
    {
      "job_id": "bet-001",
      "analysis": "bet",
      "bundle_ids": ["bundle-001"],
      "settings": {
        "retention": "peak_max",
        "concentration": "eluted"
      }
    }
  ]
}
```

Supported analysis values are `bet`, `dispersive`, `acid_base`, and
`full_peak`. Every job requires exactly one bundle. The standalone
`igc-full-peak` command retains its reviewed multi-block interface, but a batch
manifest cannot use it to declare or create a cross-bundle fit.

Every setting is allow-listed and forwarded to the corresponding public
command. Important keys include:

- BET: `probe`, `retention`, `concentration`, `origin`, `p0_min`, `p0_max`,
  `ambient_pressure_pa`, `pressure_correction`, and `sensitivity`;
- dispersive: `homologous_probe_ids`, `ambient_pressure_pa`,
  `pressure_correction`, `extrapolate`, `max_temperature_span_k`, and
  `max_flow_relative_span`;
- acid/base: the dispersive settings plus at least three explicit
  `polar_probe_ids`; and
- full peak: `probe`, `transport_mode`, `models`, `n_cells`, `n_starts`,
  `lodo`, `lodo_models`, and `cross_section_m2`.

Unknown keys, duplicate IDs, duplicate dataset aliases, missing bundle
references, unsupported analyses, and invalid setting types fail closed.

## Execution and outputs

All bundles are structurally validated before analysis begins. Jobs then call
the existing supported command paths, so batch execution cannot bypass their
scientific readiness checks, source-attributed property requirements, QC, or
reportability rules.

The requested output path must not exist. Every job runs under a temporary
staging directory; any failed job removes the entire staged batch. A successful
run atomically creates:

```text
batch-output/
├── README.md
├── batch_run.json
├── batch_summary.csv
├── bet-001/
├── dispersive-001/
├── acid-base-001/
└── full-peak-001/
```

Each job directory is the output of its underlying command using the declared
settings. Batch execution contains process-global plotting state so job order
does not change another command's figures. The batch run record retains
schema/package versions, dataset IDs, neutral-manifest
digests, job-to-bundle mappings, job-relative result locations, and the
workflow-specific reportability scope. It omits local input paths.

There is deliberately no batch-level scientific verdict. A BET SSA,
dispersive profile, acid/base profile, and full-peak recovered SSA have
different reportability meanings. The summary exposes each verdict and scope
without combining them.

Every batch job accepts exactly one neutral bundle. In particular,
`igc-report` does not create a joint full-peak fit across bundles. A reviewed
study that needs the standalone `igc-full-peak` multi-block interface must
declare and justify that relationship outside this batch manifest.

## Synthetic verification

The installed command provides a path-free smoke test using only packaged
synthetic bundles:

```bash
igc-report --synthetic-example --output output/synthetic-batch
```

This runs BET, dispersive, acid/base, and a reduced full-peak model comparison.
It verifies installation and orchestration, not experimental reportability or
a study design.
