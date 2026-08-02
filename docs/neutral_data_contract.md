# Neutral data ingestion boundary

The public analysis package accepts `igc-neutral-data/0.2.0` bundles. Protected-file access, source-format parsing, source table joins, credentials, and extraction-specific tests belong outside this repository.

## Validation before analysis

`igc_sea.io.read_neutral_bundle()` runs the bundled dependency-free validator before exposing any table. Validation covers manifest hashes and row counts, exact headers, controlled units and values, ordering, identifiers, foreign keys, trace continuity, calibration consistency, and cross-table scientific provenance.

Structural validity is necessary but not sufficient for an analysis. Each workflow must separately check that its required probes, calibration, conditions, dead-time relationship, sample metadata, and study design are present.

## Full-peak path

`igc-full-peak --neutral-bundle LABEL=PATH` is the first end-to-end consumer. It:

- preserves the complete raw detector trace and acquisition sequence;
- obtains sample mass, flow, temperature, pressure drop, probe identity, saturation pressure, detector gain, and calibration from declared neutral fields;
- integrates according to the calibration's declared area unit;
- characterizes dead-time transport separately for each supplied block;
- records neutral dataset identifiers, contract versions, and source fingerprints without copying local source paths into the run record.

One bundle currently represents one experiment and one acquisition block. Multiple bundles may be analyzed jointly only when the caller explicitly supplies separate block labels; a future study-design package must define experimental relationships and replication.

## Compatibility

Contract `0.1.0` is an immutable experimental baseline but is not accepted by the current reader. Contract `0.2.0` introduced breaking calibration and condition-provenance changes. Inputs are never silently coerced between versions.

Source-specific readers are outside this repository. Additional public analysis
entry points will be enabled only after they consume the neutral contract.
