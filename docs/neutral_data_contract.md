# Neutral data ingestion boundary

The public analysis package accepts `igc-neutral-data/0.2.0` bundles. Protected-file access, source-format parsing, source table joins, credentials, and extraction-specific tests belong outside this repository.

## Validation before analysis

`igc_analysis.io.read_neutral_bundle()` runs the bundled dependency-free validator before exposing any table. Validation covers manifest hashes and row counts, exact headers, controlled units and values, ordering, identifiers, foreign keys, trace continuity, calibration consistency, and cross-table scientific provenance.

Structural validity is necessary but not sufficient for an analysis. Each workflow must separately check that its required probes, calibration, conditions, dead-time relationship, sample metadata, and study design are present.

## Implemented consumer paths

The implemented end-to-end consumers are:

- `igc-bet --neutral-bundle PATH`, for corrected BET analysis of one block;
- `igc-dispersive --neutral-bundle PATH`, for coverage-resolved Dorris--Gray
  analysis of one homologous probe series;
- `igc-full-peak --neutral-bundle LABEL=PATH`, for nonlinear full-trace
  analysis of one or more independently characterized blocks.

All validate the neutral bundle before analysis.

### BET readiness

`igc-bet` requires exactly one experiment, one acquisition block, one uniquely
selected analyte probe, at least one dead-time trace, calibrated probe traces,
sample mass, molecular cross-section, per-injection saturation vapor pressure,
temperature, and a uniquely attributed flow channel for measured flow. Pressure
correction additionally requires either absolute inlet pressure or pressure drop;
outlet pressure may be declared or supplied explicitly as the ambient-pressure
CLI setting. When both absolute pressures and pressure drop are declared, they
must agree. Legacy loop concentration additionally requires
`injection_loop_volume_m3`.

The command retains source-neutral injection IDs and property, calibration,
flow-channel, vapor-pressure, and pressure-basis provenance in its outputs. A
structurally valid bundle may still produce a non-reportable BET result.

### Dispersive readiness

`igc-dispersive` requires exactly one experiment and one block; at least three
homologous analytes with unique carbon numbers and at least three calibrated
coverage points each; complete probe and dead-time traces; sample mass;
supplied SSA and its source; probe cross-sections and property sources;
calibration sources; positive nominal coverages; and stable measured
temperature and flow. The dispersive pressure provenance gate is stricter than
the general contract and the current BET consumer.

If more than three analytes declare carbon numbers, the caller must explicitly
select the homologous series with repeatable opaque `--homologous-probe-id`
values. Pressure correction requires measured inlet pressure or measured
pressure drop; outlet pressure may be measured or use the configured ambient
pressure, and the resolved roles are recorded. The dispersive consumer requires
declared detector gain for all required injections. Multiple gains are
preserved and warned, not silently normalized, and make the profile
non-reportable until the declared calibrations are reviewed across them.

The consumer calculates actual coverage from calibrated amount and the supplied
SSA. It never substitutes a source-specific default. It records all
interpolation and extrapolation decisions, and it marks incomplete three-probe
fits non-reportable.

### Full-peak readiness

The full-peak consumer:

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
