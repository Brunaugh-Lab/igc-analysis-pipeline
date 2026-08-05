# IGC neutral data contract 0.2.0 field dictionary

This dictionary describes scientific meaning. `schema.json` is authoritative for column names, primitive types, required fields, controlled values, keys, and numeric bounds.

## `manifest.json`

| Field | Meaning |
|---|---|
| `contract_name` | Must equal `igc-neutral-data`. |
| `contract_version` | Must equal `0.2.0`. |
| `profile` | Required table/validation profile; version 0.2.0 defines `trace-core`. |
| `dataset_id` | Immutable opaque identifier shared with `experiment.csv`. |
| `created_at` | ISO 8601 creation timestamp with timezone. |
| `adapter_version` | Version of the data-preparation software that emitted the normalized bundle. |
| `source_fingerprint` | Optional lowercase SHA-256 fingerprint of the source observations; it contains no filename or path. |
| `files` | Object keyed by required filename, with `row_count` and lowercase SHA-256 digest for each CSV. |

## `experiment.csv`

One row identifies the normalized acquisition experiment and its sample-level quantitative inputs.

| Field | Meaning |
|---|---|
| `dataset_id` | Must match the manifest. |
| `experiment_id` | Opaque identifier referenced by all injections. |
| `sample_id` | Neutral sample-registry join key. Public fixtures use synthetic identifiers. |
| `sample_mass_g` | Sample mass used in calculations, in grams. |
| `injection_loop_volume_m3` | Injection-loop volume in cubic metres when a concentration convention uses it. |
| `column_id` | Packed-column identifier in `columns.csv`. |
| `acquisition_started_at` | ISO 8601 timestamp with timezone when available. |
| `specific_surface_area_m2_g` | Externally supplied SSA used for coverage calculations; it is not a BET output calculated from this bundle. |
| `surface_area_source` | Citation, measurement identifier, or explicit modeled-estimate designation. Required when SSA is present. |

## `columns.csv`

One row describes the packed column and the experimental role it played.

| Field | Meaning |
|---|---|
| `column_id` | Opaque packed-column identifier. |
| `packing_replicate_id` | Identifier for an independently packed column; not automatically a biological replicate. |
| `column_role` | `sample`, `matrix_control`, or `reference`. |
| `internal_diameter_m` | Internal column diameter in metres, used when carrier velocity is calculated. |
| `packed_bed_length_m` | Packed-bed length in metres when known. |
| `conditioning_description` | Neutral protocol identifier or non-sensitive description. |
| `sample_batch_id` | Join key to an external sample or process registry. |
| `density_kg_m3` | Supplied or measured density in kilograms per cubic metre. |
| `density_basis` | `bulk_packed`, `envelope`, `skeletal`, or `modeled`; required whenever density is present. |

## `conditions.csv`

Each row is one condition value for one injection. Multiple rows may preserve before/after observations, measured values, targets, conversions, or interpolation.

| Field | Meaning |
|---|---|
| `condition_id` | Unique condition-observation identifier. |
| `injection_id` | Injection to which the condition applies. |
| `quantity` | `column_temperature`, `flow_standard`, `flow_column`, `pressure_inlet`, `pressure_outlet`, `pressure_drop`, or `relative_humidity`. |
| `value` | Numeric condition value. |
| `unit` | `K`, `m3_s`, `Pa`, or unitless `fraction`, as appropriate for the quantity. |
| `value_role` | `measured` or `target`. |
| `measured_at` | ISO 8601 observation timestamp when available. It is normally absent for a declared target. |
| `measurement_basis` | `direct`, `before`, `after`, `mean_before_after`, `interpolated`, `converted`, `converted_before`, `converted_after`, or `declared_target`. The converted-before/after values identify reconstructed observations that retain the source observation timing basis. |
| `source_channel` | Optional neutral acquisition-channel identifier, required by the adapter for measured flow values when channel selection occurred. |

Allowed quantity/unit pairs are fixed: temperature/K; either flow/m3_s; any pressure/Pa; relative humidity/fraction. Relative humidity spans 0–1 inclusive.

`flow_standard` is a volumetric flow referenced to exactly 273.15 K;
`flow_column` is the volumetric flow at the declared column temperature.
`pressure_inlet` and `pressure_outlet` are absolute pressures. `pressure_drop`
is the nonnegative inlet-minus-outlet differential. Analyses may cross-check
these redundant pressure declarations and reject inconsistent values.

## `injections.csv`

Each row is one acquisition event, independent of how many chemical components it contains.

| Field | Meaning |
|---|---|
| `experiment_id` | Parent experiment. |
| `injection_id` | Opaque injection identifier. |
| `block_id` | Explicit scientific/acquisition block for bracketing and joint fits. |
| `sequence_index` | Zero-based acquisition order, unique within the experiment. |
| `acquired_at` | ISO 8601 injection timestamp with timezone when available. |
| `role` | `probe`, `dead_time`, `blank`, or `reference`. |
| `target_coverage_fraction` | Declared acquisition target as a nonnegative fraction. It is not calibrated actual coverage. |
| `detector_gain` | Positive detector gain when relevant to response calibration or clipping assessment. |
| `detector_channel` | Stable neutral channel identifier that matches its trace rows. |
| `clipping_observed` | Optional adapter-observed acquisition saturation flag. Public analysis performs its own clipping QC. |

## `injection_components.csv`

One or more rows declare the chemical components of each injection.

| Field | Meaning |
|---|---|
| `injection_id` | Parent injection event. |
| `component_index` | Zero-based order within the injection. |
| `probe_id` | Probe/property identity in `probe_properties.csv`. |
| `component_role` | `analyte`, `dead_time_marker`, `carrier_component`, or `reference_component`. |
| `target_amount_mol` | Declared target amount in moles, not a calibrated result. |
| `calibration_id` | Calibration used by public analysis to convert integrated detector area to moles. |
| `saturation_vapor_pressure_Pa` | Vapor pressure evaluated for this component at the injection temperature. |
| `vapor_pressure_source` | Citation or reference identifier for the value/model. |
| `vapor_pressure_model_id` | Explicit equation, convention, or model identifier used in evaluation. |

For a probe component, vapor-pressure fields are required only when a requested analysis uses $P/P_0$. All three vapor-pressure fields occur together.

## `traces.csv`

Each row is one original detector observation. The compound primary key is injection, detector channel, and point index.

| Field | Meaning |
|---|---|
| `injection_id` | Parent injection. |
| `point_index` | Zero-based original point order; physically stored in contiguous increasing order within each injection/channel trace. |
| `time_s` | Nonnegative elapsed time in seconds; strictly increasing within a trace. |
| `detector_channel` | Must match the injection's declared channel in version 0.2.0. |
| `signal_raw` | Original finite detector response at exported resolution. |
| `signal_unit` | Explicit physical unit or `arbitrary_unit`; never imply calibration that did not occur. |
| `signal_corrected` | Optional diagnostic corrected response. It never replaces `signal_raw`. |
| `preprocessing_method` | Required when corrected response is present. |
| `preprocessing_version` | Required when corrected response is present. |

## `probe_properties.csv`

Each row is one source-attributed probe-property set.

| Field | Meaning |
|---|---|
| `probe_id` | Stable property identifier. |
| `probe_name` | Neutral chemical name. |
| `molar_mass_g_mol` | Molar mass in grams per mole. |
| `cross_section_m2` | Molecular cross-sectional area in square metres. |
| `gamma_l_d_mJ_m2` | Dispersive liquid surface tension in millijoules per square metre. |
| `donor_number_kJ_mol` | Donor-number term in kilojoules per mole under the selected convention. |
| `acceptor_number_kJ_mol` | Acceptor-number term in kilojoules per mole under the selected convention. |
| `carbon_number` | Carbon number used by a declared homologous-series method. |
| `properties_source` | Citation, reference identifier, or explicit user-supplied designation. |

## `calibration.csv`

Each row declares one detector-area-to-amount model. Parameters are interpreted only through `calibration_model`:

$$
\begin{aligned}
\text{linear: } & n = p_0 + p_1 A \\
\text{quadratic: } & n = p_0 + p_1 A + p_2 A^2 \\
\text{power law: } & n = p_0 A^{p_1}
\end{aligned}
$$

where $A$ uses `area_unit` and $n$ uses `amount_unit`.

| Field | Meaning |
|---|---|
| `calibration_id` | Stable calibration identifier. |
| `probe_id` | Probe to which the calibration applies. |
| `calibration_model` | `linear`, `quadratic`, or `power_law`. |
| `parameter_0` | Model parameter $p_0$. |
| `parameter_1` | Model parameter $p_1$. |
| `parameter_2` | Model parameter $p_2$; required only for a quadratic model. |
| `area_unit` | Exact detector-area convention consumed by the parameters. For `uV_min`, integrate a `uV` trace using minutes (or divide an area integrated in `uV_s` by 60) before applying the calibration. |
| `amount_unit` | `mol` in version 0.2.0. |
| `calibration_source` | Calibration identifier or explicit provenance. |

## Scientific-readiness boundaries

Structural validity does not establish scientific readiness. Analysis software must additionally check, at minimum:

- peak shape: sufficient pre/post-peak baseline, full trace, detector context, dead-time relationship, and clipping status;
- BET: adequate qualifying injections, calibrated amount, pressure convention, saturation pressure, mass, and model/QC settings;
- dispersive surface energy: homologous probe set, cross-sections, calibrated amount, target/actual coverage inputs, SSA provenance, and stable conditions;
- acid-base surface energy: declared polar-probe inclusion policy and source-attributed donor/acceptor properties;
- full-peak fitting: dead-time bracketing, transport inputs, calibrated response, block structure, and sensitivity/cross-validation plan;
- RH/temperature/flow maps: explicit condition provenance, conditioning history, independently packed columns/batches, and a separate declared study design.

Compressibility-corrected retention analysis additionally requires inlet and outlet absolute pressures, or pressure drop plus one absolute pressure from which the other can be reconstructed. Pressure drop alone is not scientifically sufficient.

Version 0.2.0 still represents one experiment and one adapter-selected acquisition block per bundle. Study-level relationships and joint multi-block fitting require a future study package; analysis must not silently pool separate bundles or blocks.
