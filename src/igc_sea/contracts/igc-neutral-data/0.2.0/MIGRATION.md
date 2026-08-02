# Migration from 0.1.0 to 0.2.0

Version 0.2.0 is intentionally breaking. Version 0.1.0 remains unchanged.

- Calibration columns `coefficient_0`, `coefficient_1`, and `coefficient_2` become `parameter_0`, `parameter_1`, and `parameter_2`.
- `calibration_model` adds `power_law`, interpreted as $n=p_0A^{p_1}$.
- `conditions.csv` adds the optional neutral `source_channel` field and requires it for measured flow values.
- `measurement_basis` adds `converted_before` and `converted_after` for reconstructed observations tied to before/after source readings.
- `manifest.json` may include an optional lowercase SHA-256 `source_fingerprint`.

Consumers must select a contract version explicitly. Do not pass a 0.2.0 bundle to the 0.1.0 validator or silently rename headers.
