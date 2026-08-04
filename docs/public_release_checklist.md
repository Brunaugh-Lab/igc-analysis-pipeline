# Public-release checklist

This checklist is the visibility-change gate for `igc-analysis-pipeline`.
Passing CI is necessary but does not by itself authorize publication.

## Current blocking decisions

- [x] GitHub Support confirmed server-side cache clearance and garbage
  collection for superseded, unreachable objects on 2026-08-02.
- [x] Authenticated direct requests for every affected former commit and file
  view return `404` or `File not found`. Exact identifiers are retained only in
  the lab's private governance record.
- [x] The code owner explicitly approved release under MIT without a separate
  University sign-off. `LICENSE`, SPDX package metadata, and the README notice
  identify Ashlee Brunaugh as the copyright holder.
- [x] The public repository and distribution name are
  `igc-analysis-pipeline`; the Python import name is `igc_analysis`.
- [x] Sole-author and maintainer metadata are present in `pyproject.toml` and
  `CITATION.cff`, including Ashlee Brunaugh's ORCID and U-M contact.

Do not change repository visibility while any item above is unresolved.

## Content boundary

- [x] Only source-neutral analysis begins in this repository.
- [x] No acquisition-file access, source-specific mapping, credentials, or
  private extraction tests are tracked.
- [x] No governed or experimental datasets or derived study outputs are
  tracked.
- [x] Every bundled example is synthetic and uses fictional, opaque identifiers.
- [x] No local filesystem paths, collaborator identifiers, or unpublished study
  labels appear in tracked files or built distributions.
- [x] `python scripts/check_public_release_boundary.py` passes. CI rejects
  reserved data paths, non-synthetic tabular/database files, absolute user
  paths, and known source-specific markers.

## Scientific and software verification

- [x] The locked test suite and lint checks pass.
- [x] Both source and wheel distributions build successfully.
- [x] The installed wheel completes the synthetic full-peak workflow.
- [x] The neutral bundle validator rejects malformed hashes, units, keys,
  ordering, and provenance.
- [x] Generated interpretation remains conditional on the numerical and
  identifiability checks it reports.
- [x] The supported public commands were reviewed for scientific-reference and
  decision-threshold provenance. Development-only calculation modules remain
  explicitly outside the advertised public workflows until their scientific
  documentation is complete.

## Verification evidence — 2026-08-02

- Reachable-history scans found no absolute user paths, protected source-format
  markers, removed-data identifiers, private-key headers, or common token
  patterns. The only tracked tabular data objects are the declared synthetic
  neutral-contract fixture.
- `264` tests passed; Ruff and lock integrity passed.
- `CITATION.cff` passed the Citation File Format `1.2.0` schema validator.
- Source and wheel distributions contain the MIT license; the source
  distribution also contains `CITATION.cff`.
- Extracted distribution scans found no private paths or source-specific
  markers and no tabular data outside the synthetic fixture.
- The installed wheel completed the synthetic full-peak workflow; installed
  wheel and source-distribution command smoke tests passed.
- The exact locked runtime dependency audit reported no known vulnerabilities.

The evidence above applies to archived release `v2026.7.31`. The separately
dated gate below covers the BET command added in release candidate `v2026.8.4`.

## BET candidate gate — 2026-08-04

- [x] Resolve every P1/P2 finding from the initial independent read-only audit.
- [x] Demonstrate a nondegenerate synthetic Type II fit through the installed
  wheel and a synthetic Type III abstention in the workflow regression suite.
- [x] Obtain a second independent read-only review of the remediated diff.
  Claude Code reported no P0/P1 findings and confirmed every original P0--P2
  item fixed by execution. Its two P2 hardening recommendations were resolved
  with single-source package versioning and a full-peak clipping regression.
- [x] Complete governed legacy-versus-neutral equivalence against the final
  corrected historical workflow without committing governed data or outputs.
  The 2026-08-04 positive-control comparison agreed to floating-point precision
  across the numerical fit, pressure correction, QC, classification, and
  reportability; sensitive evidence remains outside this repository.
- [x] Re-run the complete suite, lint, lock, working-tree boundary, artifact,
  installed-wheel, and dependency checks. On 2026-08-04, 285 tests passed;
  Ruff, lock integrity, and diff checks passed; all 79 tracked and candidate
  paths passed the content boundary; both 2026.8.4 release distributions built and
  passed artifact scans; the isolated wheel recovered the packaged Type II
  example with clean QC; and the locked runtime audit found no known
  vulnerabilities. Repeat the boundary check against the staged tree before
  commit and rely on GitHub CI for release evidence.
- [x] Update `CITATION.cff`, tag the immutable release, and record the new
  version DOI only after the release exists. Release `v2026.8.4` points to
  verified commit `c9c7465`; GitHub Actions run `30939032517` passed; Zenodo
  archived the release at version DOI `10.5281/zenodo.21796032`; and the concept
  DOI remains `10.5281/zenodo.21760976`.

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

## Dispersive candidate gate — 2026-08-04

- [x] Reconstruct the production dispersive call path from final corrected
  historical commit `52d1cf9`, including calibrated actual coverage, explicit
  SSA, matched CoM/peak-max dead times, measured conditions, James--Martin
  correction, production `extrapolate=True`, Dorris--Gray fitting, and
  $W_{cohesion}^d=2\gamma_s^d$.
- [x] Implement `igc-dispersive` as a neutral `0.2.0` consumer with explicit
  property/calibration/SSA provenance, transactional outputs, strict JSON,
  PDF/PNG diagnostics, and no local input paths.
- [x] Add a reproducible closed-form fixture with nonzero pressure correction
  and asymmetric peaks. The final isolated wheel recovered 39.89--40.04 mJ/m²
  from a true 40.0 mJ/m² CoM profile; peak maximum differed by approximately
  2 mJ/m² as an actual sensitivity calculation.
- [x] Complete two independent read-only Claude Code audits. The first audit's
  apparent extrapolation blocker was resolved by reading the removed production
  orchestrator; the production call explicitly enabled extrapolation. The
  public candidate retains that numerical comparison but makes any extrapolated
  profile non-reportable. All subsequent P1/P2 findings were remediated.
- [x] Exercise malformed provenance, target-only conditions, dead-time flow
  drift, clipping, co-injected analytes, changing detector gain, sparse and
  out-of-range coverage, missing cross-section, nonuniform sampling, empty-fit
  plotting, generator reproducibility, and partial-output prevention.
- [x] Re-run local verification after remediation: 308 tests passed; Ruff,
  lock, diff, neutral-validator, and working-tree disclosure-boundary checks
  passed; source and wheel distributions built; and the isolated wheel
  completed the packaged example with clean QC and positive reportability.
- [ ] Complete an authorized governed legacy-versus-neutral equivalence run.
  Keep the neutral bundle and comparison outputs outside Git.
- [ ] Commit and push the candidate, obtain passing GitHub CI, tag a new
  immutable release, archive it in Zenodo, and update the student guide from
  "next tagged release" to the exact version and DOI.

## Post-release verification

- [x] Clone without authentication into a new temporary directory.
- [x] Run the documented installation and synthetic example from that clone.
- [x] Inspect the public repository file list, release artifacts, and metadata.
- [x] Confirm the private extraction repository and its implementation remain
  inaccessible to unauthenticated users.

Verified on 2026-08-02 at commit `e530475`: GitHub reports the repository as
public, a clone with credential helpers disabled succeeded, the cloned boundary
check passed, and the documented synthetic full-peak workflow completed. An
unauthenticated GitHub API request returned `200` for this repository and `404`
for the private extraction repository. GitHub Actions run `30752878464` passed
the complete Python 3.10--3.14 release matrix before the visibility change.
