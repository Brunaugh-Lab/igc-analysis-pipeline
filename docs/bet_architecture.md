# Source-neutral BET architecture

Status: released in `v2026.8.4` and included in `v2026.8.6`; independent review,
synthetic recovery, and regression verification are complete. Entry point:
`igc-bet`.

## Boundary

The workflow accepts one validated `igc-neutral-data/0.2.0` bundle. It does not
read acquisition files or infer source schemas. Extraction and source-specific
normalization remain outside this repository.

## Required declarations

- exactly one experiment and one acquisition block;
- one uniquely selected analyte probe and at least one dead-time marker;
- complete detector traces and a compatible area-to-moles calibration;
- sample mass and probe molecular cross-section;
- evaluated saturation vapor pressure for every probe injection;
- column temperature and either column flow or standard flow from one declared
  measured source channel, with target values used only as explicit fallback;
- absolute inlet/outlet pressures or pressure drop plus an explicit ambient
  absolute pressure when pressure correction is enabled;
- injection-loop volume only for the legacy loop-concentration comparison.

`flow_standard` is referenced to 273.15 K. Inlet and outlet pressures are
absolute; pressure drop is inlet minus outlet. Redundant declarations must
agree within the documented numerical tolerance or the workflow stops.

## Corrected calculation path

1. Validate the complete bundle and its per-file hashes.
2. Characterize peak-maximum and center-of-mass dead times from the declared
   dead-time traces.
3. Baseline-correct each probe trace and apply its declared calibration.
4. Use the eluted peak-apex concentration by default:

   $$
   c_{apex}=\frac{S_{apex}(n_{inj}/A)}{F_{col}}
   $$

5. Compute relative pressure from the evaluated per-injection vapor pressure:

   $$
   P/P_0=\frac{c_{apex}RT}{P_{sat}}
   $$

6. Apply the James-Martin pressure-gradient correction:

   $$
   V_N=jF_{col}(t_R-t_0),\qquad
   j=\frac{3}{2}\frac{(p_i/p_o)^2-1}{(p_i/p_o)^3-1}
   $$

7. Integrate $V_N/m=dq/dc$, perform the adaptive BET linearization, and run
   Rouquerol, retention, isotherm-shape, sensitivity, and saturation checks.
8. Classify isotherm shape and issue a strict reportability verdict. A numeric
   fit alone is never permission to report an SSA.

Peak-maximum retention subtracts peak-maximum dead time; center-of-mass
retention subtracts center-of-mass dead time. The legacy loop-concentration and
alternative origin conventions are explicit comparison modes, not silent
defaults.

## Outputs and provenance

The command writes injection, isotherm, and BET-linearization tables; PDF and
PNG diagnostics; a strict-JSON run record; and a conditional README. Outputs
record the neutral dataset ID, contract version, manifest digest, source
fingerprint when supplied, injection IDs, dead-time IDs, sample ID, declared
property/calibration/vapor-pressure sources, flow channel, pressure basis,
software version, conventions, QC, sensitivity, and reportability.

Local input paths are deliberately omitted.

## Verification

Unit tests, closed-form recovery, adversarial QC cases, installed-distribution
smoke tests, and regression comparisons cover the numerical fit, pressure
correction, classification, and reportability path. Application to a new study
still requires review of that study's inputs, diagnostics, and QC verdict.

The packaged closed-form fixture deliberately isolates the BET calculation. It
uses zero pressure drop, symmetric probe peaks, and identical dead-time markers,
so James--Martin sensitivity, retention-convention sensitivity, and
dead-time-stability QC are exercised by separate adversarial tests rather than
by the packaged recovery example.
