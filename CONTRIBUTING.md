# Contributing

Contributions to the IGC Analysis Pipeline should be scientifically explicit,
reproducible, and portable across supported environments.

## Repository contents

Do not commit:

- study datasets or generated analysis results;
- machine-specific paths;
- personal identifiers or unpublished study labels; or
- files unrelated to the documented software and synthetic examples.

Use synthetic fixtures with clearly fictional identifiers. The directories
`governed_data/`, `local_bundles/`, and `data/examples/` are reserved for local
work and must never contain tracked files. Do not bypass these protections with
forced staging.

## Development checks

Create an isolated environment and install the development dependencies as
described in `README.md`, then run:

```bash
python -m pytest
ruff check .
python scripts/check_public_release_boundary.py
python -m build
```

All tests and examples in this repository must run from packaged synthetic
fixtures without network access.

## Scientific changes

For changes that affect reported quantities or interpretation:

- state the scientific assumption and applicable units;
- add synthetic tests for normal and failure paths;
- preserve provenance fields and explicit reportability gates;
- distinguish technical observations from independent replicates; and
- avoid mechanistic conclusions that are not directly tested.

Do not silently relax validation, identifiability, mass-balance, or
surface-area-reporting guardrails.

## Before opening a pull request

Confirm that the working tree contains only intended source, documentation, and
synthetic fixtures. Review the complete diff and run all development checks.
Describe the scientific and software consequences of the change in the pull
request.
