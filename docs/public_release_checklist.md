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
- [x] Exercise malformed provenance, target-only conditions, measured-pressure
  enforcement, explicit homolog selection, dead-time flow drift, clipping,
  co-injected analytes, changing detector gain, sparse and
  out-of-range coverage, missing cross-section, nonuniform sampling, empty-fit
  plotting, generator reproducibility, and partial-output prevention.
- [x] Re-run local verification after remediation: 313 tests passed; Ruff,
  lock, diff, neutral-validator, and working-tree disclosure-boundary checks
  passed; source and wheel distributions built; and the isolated wheel
  completed the packaged example with clean QC and positive reportability.
- [x] Complete an authorized governed legacy-versus-neutral equivalence run,
  with the neutral bundle and comparison outputs kept outside Git. All 24
  homologous injections matched for actual coverage; peak-maximum retention
  volume agreed to floating-point precision; center-of-mass retention volume
  differed by at most 0.00072 mL/g (0.0037% relative) after time-weighting
  correction; and the nine-point primary surface-energy profile differed by at
  most 0.001 mJ/m2 (0.0028% relative). This was one authorized acquisition,
  not a claim of universal equivalence.
- [x] Commit and push the candidate and obtain passing GitHub CI. Checkpoint
  `6f9dd6a` passed run `30962983337`, including Python 3.10--3.14, packaging,
  installed synthetic workflows, lint, boundary, and dependency jobs.
- [x] Tag immutable release `v2026.8.5`, archive it in Zenodo, and record the
  exact version DOI. The tag points to verified commit `c45daec`; release-gate
  CI run `30964023144` passed; Zenodo version DOI is
  `10.5281/zenodo.21799047`; and the concept DOI remains
  `10.5281/zenodo.21760976`.

## Acid/base candidate gate — 2026-08-04

- [x] Reconstruct the final historical Schultz--Gutmann path: all selected
  probes mapped from calibrated actual coverage to a common target grid,
  Schultz reference lines at each coverage, polar-probe $\Delta G_{sp}$, and
  regression-derived Ka/Kb. Van Oss was historically opt-in and remains
  outside this source-neutral milestone.
- [x] Implement `igc-acid-base` with explicit opaque homologous and polar probe
  selections, source-attributed property requirements, CoM retention as
  primary, peak maximum as the historical sensitivity calculation,
  transactional outputs, strict JSON, figures, and interpretation boundaries.
- [x] Add a deterministic detector-trace fixture that recovers Ka=0.03 and
  Kb=0.05 with at least three probes and nontrivial CoM/peak-maximum
  sensitivity. The installed wheel completed the fixture with positive
  reportability.
- [x] Complete an authorized governed comparison outside Git. With identical
  interpolated values and the identical historical property set, the new and
  historical calculation layers matched Ka, Kb, and R² exactly. No governed
  input, detailed result, or comparison artifact was retained in this
  repository.
- [x] Add and test a public input-unit/plausibility gate for declared liquid
  dispersive tensions. The public consumer rejects inconsistent values rather
  than silently converting them.
- [x] Resolve every P0--P2 finding from an independent public-repository audit.
  Claude Code identified an empty-result crash, misleading invalid-retention
  count, missing Schultz-line QC, undeclared direct NumPy dependency, and
  missing failure-path tests. All were remediated, and the final narrow audit
  reported no P0--P2 findings.
- [x] Re-run the complete suite, lint, lock, boundary, artifact, installed-wheel,
  and dependency gates after audit remediation. On 2026-08-04, 327 tests
  passed; Ruff, lock integrity, diff checks, and the public content boundary
  passed; fresh wheel and source distributions built; installed dispersive and
  acid/base synthetic workflows were reportable; and the locked runtime audit
  found no known vulnerabilities.
- [x] Commit and push the final candidate and obtain passing GitHub CI.
  Checkpoint `6f9dd6a` passed run `30962983337`, including the complete Python
  3.10--3.14 matrix and every release gate.
- [x] Include the command in immutable release `v2026.8.5`, archive it in
  Zenodo, and record version DOI `10.5281/zenodo.21799047`.

## Neutral batch-report candidate gate — 2026-08-05

- [x] Recover only the safe orchestration requirements from the final
  historical workflow. Source-folder scanning, filename-derived replicate
  grouping, source coupling, auto-pooling, and study-design inference were not
  restored.
- [x] Implement `igc-report` with packaged `igc-analysis-batch/0.1.0`
  manifests, one validated neutral bundle per job, allow-listed settings, and
  the existing BET, dispersive, acid/base, and full-peak validators.
- [x] Publish transactionally only after every job succeeds. Existing,
  nested, or concurrently changed output locations fail closed without
  deleting unexpected content.
- [x] Keep reportability workflow-specific. The batch record contains no local
  paths, does not combine scientific verdicts, and labels full-peak QC
  `NOT_COMBINED` rather than inventing a pass.
- [x] Complete independent Claude Code audit, remediation, re-audit, and final
  focused check. The initial plotting-state, pooling-scope, output-safety,
  QC-summary, numeric-domain, and test-coverage findings were resolved; the
  final check reported no P0--P2 findings.
- [x] Pass the complete local gate at `ee3cba6`: 352 tests, Ruff, lock and diff
  checks, public-boundary scan, source/wheel builds, installed four-workflow
  runs from both distributions, and persisted-output path scans.
- [x] Obtain passing GitHub Actions run `31008560968`, including Python
  3.10--3.14, dependency audit, lock/lint/public boundary, distribution builds,
  and installed synthetic workflows.
- [ ] Tag immutable release `v2026.8.6`, archive it in Zenodo, record the exact
  version DOI after the archive exists, and verify a credential-free install
  of the tag before directing students to batch reporting.

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

Release `v2026.8.5` was verified on 2026-08-05 at commit `c45daec`. GitHub
Actions run `30964023144` passed the complete Python 3.10--3.14 matrix and all
packaging, installed-workflow, lint, boundary, and dependency jobs. A fresh
credential-free clone of the immutable tag installed as version `2026.8.5` and
completed the packaged BET, dispersive, and acid/base synthetic workflows with
positive reportability. Zenodo published version DOI
`10.5281/zenodo.21799047` under concept DOI `10.5281/zenodo.21760976`.
