# Public-release checklist

This checklist is the visibility-change gate for `igc-analysis-pipeline`.
Passing CI is necessary but does not by itself authorize publication.

## Current blocking decisions

- [ ] GitHub Support has confirmed server-side removal of the superseded Git
  objects and cached views that contained governed experimental data.
- [ ] Direct requests for each affected old commit and file view no longer
  resolve.
- [ ] The code owner has confirmed authority to apply a public software license,
  and the approved `LICENSE` file and package metadata are present.
- [x] The public repository and distribution name are
  `igc-analysis-pipeline`; the Python import name is `igc_analysis`.
- [ ] Maintainer metadata and a citation record have been approved and added.

Do not change repository visibility while any item above is unresolved.

## Content boundary

- [ ] Only source-neutral analysis begins in this repository.
- [ ] No acquisition-file access, source-specific mapping, credentials, or
  private extraction tests are tracked.
- [ ] No governed or experimental datasets or derived study outputs are
  tracked.
- [ ] Every bundled example is synthetic and uses fictional, opaque identifiers.
- [ ] No local filesystem paths, collaborator identifiers, or unpublished study
  labels appear in tracked files or built distributions.
- [ ] `python scripts/check_public_release_boundary.py` passes.

## Scientific and software verification

- [ ] The locked test suite and lint checks pass.
- [ ] Both source and wheel distributions build successfully.
- [ ] The installed wheel completes the synthetic full-peak workflow.
- [ ] The neutral bundle validator rejects malformed hashes, units, keys,
  ordering, and provenance.
- [ ] Generated interpretation remains conditional on the numerical and
  identifiability checks it reports.
- [ ] Scientific references and decision thresholds have appropriate provenance
  or are explicitly identified as conventions.

## Final visibility-change procedure

1. Fetch the remote and confirm local `main` matches `origin/main`.
2. Confirm GitHub advertises only the intended branch and approved tags.
3. Re-run the content-boundary scan against tracked files and both built
   distributions.
4. Run the complete locked CI-equivalent test, lint, build, installed-wheel
   synthetic workflow, and dependency audit.
5. Confirm the GitHub Support purge using direct old-object requests; do not rely
   only on the visible branch history.
6. Record the reviewed commit and evidence in the repository architecture note.
7. Change visibility only as a separate, explicit owner-approved action.

## Post-release verification

- [ ] Clone without authentication into a new temporary directory.
- [ ] Run the documented installation and synthetic example from that clone.
- [ ] Inspect the public repository file list, release artifacts, and metadata.
- [ ] Confirm the private extraction repository and its implementation remain
  inaccessible to unauthenticated users.
