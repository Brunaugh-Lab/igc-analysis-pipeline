# Release quality and reproducibility

Each tagged release of the IGC Analysis Pipeline must satisfy the following
criteria before publication.

## Scientific behavior

- Supported commands validate their declared input contract before analysis.
- Numerical outputs retain the settings, software version, input identifiers,
  and scientific-property provenance needed to reproduce the calculation.
- QC and reportability are evaluated independently for each workflow.
- Synthetic recovery cases exercise expected numerical behavior and
  adversarial cases exercise failure and abstention paths.
- Interpretation remains conditional on model structure, identifiability, and
  the reported QC diagnostics.

## Software verification

- The complete test suite passes on every supported Python version.
- Release-blocking lint and dependency checks pass.
- Wheel and source distributions build successfully.
- Both distribution formats install in isolated environments.
- Packaged synthetic examples complete through their installed command-line
  entry points.
- Generated run records and artifacts contain portable identifiers rather than
  machine-specific paths.

## Publication

- The release version is updated consistently in package and citation metadata.
- An immutable Git tag identifies the exact tested commit.
- GitHub publishes release notes describing scientific capabilities and changes.
- Zenodo archives the tagged version and assigns a version-specific DOI.
- The README and `CITATION.cff` identify the current version DOI.

Existing tags are never moved. Changes that can alter scientific results receive
a new release, new tag, and new version DOI.
