# IGC neutral data contract 0.2.0

Status: experimental specification. The logical model may change before 1.0.0.

This directory defines the source-neutral handoff between an authorized private extraction adapter and scientific analysis software. It contains no source-format names, access mechanism, credentials, study data, or unpublished scientific results.

## Contract profiles

Version 0.2.0 defines one required profile:

- `trace-core`: a single normalized experiment with complete ordered detector traces, explicit injection roles and components, independently declared conditions, column context, probe properties, and calibration metadata.

Analysis-specific readiness is stricter than structural validity. A structurally valid bundle may still be insufficient for BET, surface-energy, peak-shape, full-peak, or other future analyses.

## Required bundle files

```text
bundle/
├── manifest.json
├── experiment.csv
├── columns.csv
├── conditions.csv
├── injections.csv
├── injection_components.csv
├── traces.csv
├── probe_properties.csv
└── calibration.csv
```

Every bundle represents exactly one experiment in version 0.2.0. Relationships among bundles—such as before/after, treatment/control, matrix/sample, process level, batch, or held-out prediction—belong in a separate future study package keyed by immutable `dataset_id` values. Multiple acquisition blocks are not pooled into one bundle by the current adapter.

## Normative rules

- CSV files use UTF-8, a comma delimiter, a header row, `.` as the decimal separator, and an empty field for a missing optional value.
- A bundle contains only `manifest.json` and the CSV files declared by its profile. Undeclared files or directories are rejected.
- Literal sentinel values such as `NA`, `N/A`, `null`, `None`, `-999`, or `.` are invalid.
- Times use seconds, temperatures use kelvin, pressures use pascals, amounts use moles, lengths use metres, and mass uses kilograms or grams exactly as stated in the field name.
- `dataset_id`, `experiment_id`, `column_id`, `injection_id`, and related identifiers are neutral opaque identifiers. They must not encode source table names, source filenames, credentials, commercial identifiers, collaborator names, or protected study paths.
- Identifier syntax and path-shaped strings are validated mechanically. Free-text provenance and descriptions still require human disclosure review before a bundle is shared outside its governed location.
- `sequence_index` records acquisition order. Injection rows are physically stored in contiguous increasing `sequence_index` order; order must not be reconstructed from filenames.
- An injection role and its chemical components are separate concepts. Version 0.2.0 accepts multiple components even though public co-injection analysis is not yet required.
- `conditions.csv` may contain both measured and target values for the same quantity. Their `value_role` and `measurement_basis` must remain explicit.
- `traces.csv` preserves the complete exported point order and sampling resolution. Rows are physically stored in contiguous increasing `point_index` order within each injection/channel trace. `signal_raw` is mandatory.
- Resampling, smoothing, baseline correction, integration, peak maximum, center of mass, width, asymmetry, tailing, clipping QC, transport fitting, retention volume, coverage, isotherm calculations, and surface-energy calculations are analysis transformations. They are not substitutes for the raw trace.
- Optional `signal_corrected` values are diagnostic only. If present, `preprocessing_method` and `preprocessing_version` are required and `signal_raw` remains mandatory.
- Target coverage and actual coverage are not interchangeable. Actual coverage is recalculated from calibrated amount, probe cross-section, sample mass, and explicitly sourced specific surface area.
- Temperature-dependent saturation vapor pressure is evaluated per injection component and includes its source and model identifier.
- Property and calibration values must declare provenance. Public analysis software must not silently replace missing values with unpublished defaults.

## Files in this specification

- `schema.json` is the machine-readable logical schema.
- `FIELD_DICTIONARY.md` is the human-readable field and scientific-meaning reference.
- `validate_bundle.py` performs structural and cross-table validation using only the Python standard library.
- `examples/synthetic_peak_shape/` is a fully synthetic valid fixture with bracketing dead-time injections and asymmetric probe traces.

Validate the fixture from this directory:

```bash
python3 validate_bundle.py examples/synthetic_peak_shape
```

The validator confirms contract conformance, not scientific reportability or model validity.
