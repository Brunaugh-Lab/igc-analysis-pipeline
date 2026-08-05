# Source-neutral acid/base architecture

## Scope

`igc-acid-base` provides a coverage-resolved Schultz reference line and
Gutmann Ka/Kb regression from a validated `igc-neutral-data/0.2.0` bundle. It
preserves the corrected shared retention path used by `igc-dispersive`:
calibrated actual coverage, declared SSA and cross-sections, measured flow and
temperature, James--Martin pressure correction, and matched dead-time
definitions.

The command does not select probes by chemical name. The user explicitly
declares at least three opaque IDs for the homologous reference series and at
least three opaque IDs for the polar-probe inclusion set.

## Calculation

At each target coverage, the homologous series defines the Schultz line

$$
RT\ln(V_N) = m\left(a\sqrt{\gamma_L^d}\right)+b.
$$

For each declared polar probe, the vertical displacement from that line is

$$
\Delta G_{sp}=RT\ln(V_N)-\left[m\left(a\sqrt{\gamma_L^d}\right)+b\right].
$$

The declared donor number, $DN$, and modified acceptor number, $AN^*$, enter
the Gutmann regression

$$
\frac{\Delta G_{sp}}{AN^*}=K_a\frac{DN}{AN^*}+K_b.
$$

Center-of-mass retention is the primary calculation. Peak maximum is emitted
as a sensitivity calculation because it was the historical acid/base
convention and can differ when peaks are asymmetric.

## Readiness and reportability

In addition to structural contract validity, the workflow requires:

- exactly one experiment and acquisition block;
- at least three calibrated coverage points for every selected probe;
- at least three homologs with unique declared carbon numbers;
- at least three explicitly included polar probes;
- source-attributed cross-section and dispersive liquid surface tension for
  every selected probe;
- source-attributed donor and modified acceptor numbers for every polar probe;
- measured per-injection flow and column temperature; and
- the pressure measurements required by the selected correction mode.

The consumer applies a 1--100 mJ/m² unit/plausibility gate to declared liquid
dispersive tensions. It does not silently convert values that appear to have
been supplied in J/m².

A profile is reportable only when there are no critical QC flags, every
coverage has a regression-derived finite Ka, Kb, and R² from at least three
polar probes, every selected-probe value is interpolated within its measured
coverage range, and the required injections use one detector gain. Negative
$\Delta G_{sp}$, negative Ka/Kb, or R² below 0.5 remain visible as review
warnings; R² below 0.3 is critical. The upstream Schultz line is independently
gated at R² 0.98 (warning) and 0.95 (critical), and its implied dispersive
surface energy uses the same expected and critical bounds as `igc-dispersive`.

Nonpositive net retention is retained as an undefined probe/mode value rather
than converted to a positive quantity. Insufficient positive observations are
a critical QC condition, so other coverages remain inspectable but the profile
is non-reportable.

## Interpretation boundary

Ka and Kb depend on the selected probe set, property convention, coverage,
retention definition, and regression quality. They are descriptive
donor/acceptor parameters, not a unique molecular mechanism. Differences may
also reflect transport, packing, calibration, peak shape, or coverage mapping.

Van Oss components are deliberately unavailable in the source-neutral command.
The current contract does not represent the required liquid acid/base
components with explicit source attribution, and silently reusing internal
name-based defaults would violate the public boundary.
