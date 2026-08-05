# Source-neutral dispersive surface-energy architecture

Status: released in `v2026.8.5` and included in `v2026.8.6`. Entry point:
`igc-dispersive`.
Synthetic closed-form recovery, independent review, and regression verification
are complete.

## Boundary

The workflow accepts one validated `igc-neutral-data/0.2.0` bundle containing
the detector traces, experimental conditions, calibration, and probe properties
required by the calculation.

## Required declarations

- exactly one experiment and one acquisition block;
- at least three selected analyte probes with unique carbon numbers; when a
  bundle also contains other carbon-numbered analytes, supply the homologous
  series as repeatable opaque `--homologous-probe-id` values;
- at least three calibrated coverage points for each probe;
- complete detector traces and at least one dead-time trace;
- sample mass, supplied specific surface area, and SSA provenance;
- probe molecular cross-section, carbon number, and property provenance;
- a compatible area-to-moles calibration and provenance for every probe;
- a positive declared detector gain for every required injection;
- a positive nominal coverage for every probe injection;
- stable, measured column temperature and flow; and
- absolute pressures, or pressure drop plus declared ambient absolute pressure,
  when pressure correction is enabled.

The inlet pressure or pressure drop must be measured. Outlet pressure may be
measured or use the configured ambient pressure. Pressure roles and the
resolved basis are retained in the run record. Detector-gain variation is
retained and flagged for review because each area-to-amount calibration must be
valid for the gain used by its injection; the profile remains non-reportable
until that review is resolved.

The supplied SSA is an experimental input used to convert calibrated amount to
actual fractional coverage. The workflow does not silently substitute an
instrument default or estimate SSA from the dispersive experiment.

## Calculation path

1. Validate the neutral bundle and per-file hashes.
2. Characterize peak-maximum and center-of-mass dead times separately.
3. Baseline-correct every probe trace and apply its declared calibration.
4. Calculate monolayer capacity for each probe,

   $$
   n_m=\frac{SSA\,m}{N_A a_{probe}},
   \qquad
   \theta_{actual}=\frac{n_{probe}}{n_m}.
   $$

5. Calculate the specific net retention volume using measured column flow and
   the James-Martin pressure-gradient factor,

   $$
   V_N=\frac{jF_{col}(t_R-t_0)}{m}.
   $$

6. Map each probe's $V_N(\theta_{actual})$ curve to the nominal coverage grid
   with piecewise-linear interpolation. Extrapolation is enabled by default to
   preserve the corrected historical convention, but every extrapolated value
   is flagged and makes the profile non-reportable. `--no-extrapolation` leaves
   out-of-range values undefined. Incomplete three-probe fits are also
   non-reportable.
7. At each coverage, fit the homologous series,

   $$
   RT\ln(V_N)=\Delta G_{CH_2}\,n_C+b,
   $$

   and calculate the Dorris-Gray dispersive surface energy,

   $$
   \gamma_s^d=\frac{(\Delta G_{CH_2})^2}
   {4N_A^2a_{CH_2}^2\gamma_{CH_2}(T)}.
   $$

8. Report center-of-mass retention as primary and peak maximum as a sensitivity
   result. Dispersive work of cohesion is reported as
   $W_{cohesion}^d=2\gamma_s^d$.
9. Run count, physical-bound, profile-shape, retention-order, clipping,
   extrapolation, and retention-definition checks. A numerical profile is not
   automatically reportable.

## Outputs and interpretation boundary

The command writes injection-level and interpolated retention tables, the
coverage-resolved profile, PDF and PNG figures, a strict-JSON run record, and a
short README. The record retains dataset and sample IDs, manifest fingerprint,
SSA/property/calibration provenance, flow channels, pressure basis, software
version, settings, QC, and reportability. Local input paths are omitted.

The profile-shape label is descriptive. A declining, flat, increasing, or
U-shaped profile is not a unique mechanistic assignment; transport, packing,
peak shape, calibration, dose nonlinearity, and extrapolation remain plausible
contributors.

## Verification boundary

The packaged fixture begins with complete synthetic detector traces and a known
40.0 mJ/m² Dorris-Gray profile. It exercises calibration, actual-coverage
recovery, matched dead-time subtraction, measured conditions, pressure
correction, interpolation/extrapolation, both retention definitions, QC, CLI
outputs, and strict provenance records.

The fixture uses nonzero pressure drop and asymmetric probe peaks, so the
James--Martin arithmetic is nondegenerate and the center-of-mass and
peak-maximum results differ. Separate adversarial tests verify fail-closed
extrapolation, target-only conditions, dead-time condition drift, clipping,
missing provenance, incomplete homologous coverage, and empty-fit plotting.
Explicit homolog selection, measured-pressure enforcement, and detector-gain
variation are also covered.
