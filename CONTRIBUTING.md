# Contributing

This repository contains source-neutral inverse gas chromatography analysis.
Contributions must preserve the boundary between public analysis and private,
authorized data extraction.

## Data and source boundaries

Do not commit:

- experimental or governed datasets;
- acquisition files or source-specific exports;
- credentials, access methods, or extraction implementation details;
- local paths, collaborator identifiers, or unpublished study labels; or
- outputs derived from non-public data.

Use synthetic fixtures with clearly fictional identifiers. The directories
`governed_data/`, `local_bundles/`, and `data/examples/` are reserved for local
work and must never contain tracked files. Do not bypass these protections with
forced staging.

Source-specific extraction and mapping tests belong in the authorized private
extraction repository. Public contributions should begin at the versioned
neutral data contract.

## Development checks

Create an isolated environment and install the development dependencies as
described in `README.md`, then run:

```bash
python -m pytest
ruff check .
python scripts/check_public_release_boundary.py
python -m build
```

All tests and examples in this repository must run without private files,
network access, or a source-specific acquisition runtime.

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
If a contribution may expose private data or source-specific implementation,
stop and contact a repository maintainer privately rather than opening a public
issue or pull request.
