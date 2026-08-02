# Full-peak nonlinear inverse chromatography — architecture

Status: implemented. Entry point `igc-full-peak` / `igc_analysis.analysis.full_peak`.

## Motivation

Conventional retention-based workflows compress each chromatogram into scalar
descriptors — peak-maximum or centre-of-mass retention time and peak area — and
then builds an isotherm from net retention volumes. That discards the peak
*shape*, which is where finite-concentration information lives: at nonlinear
loading, the isotherm curvature deforms the peak (fronting vs tailing) and
shifts its moments in a dose-dependent way.

This module keeps the complete baseline-corrected FID trace `S(t)` and fits a
forward column-transport model to all injections simultaneously, so an
adsorption isotherm is *identified* from peak shape rather than assumed from a
single retention statistic.

## Design principles

1. **Validate, then analyze.** The public path reads the source-neutral
   `igc-neutral-data/0.2.0` contract through `io.neutral_data`. Baseline and
   physical transformations remain public analysis operations.
2. **Transport is calibrated on methane, then frozen.** Adsorption parameters
   must not be free to absorb transport broadening.
3. **No invented geometry.** The packed-bed length is unrecorded, so the model
   is parameterised by methane-derived void time and an apparent plate number.
4. **A number is not a result.** Parameters are reported only with an
   identifiability verdict; SSA is gated structurally (see below).

## Module layout

| Module | Responsibility |
|---|---|
| `analysis/isotherm_models.py` | `q(c)`, `dq/dc`, bounds, units, `c→0` safeguards, `has_finite_capacity` |
| `analysis/column_model.py` | Equilibrium-dispersive forward solver; methane transport characterisation |
| `io/neutral_data.py` | Contract validation and source-neutral table loading |
| `analysis/full_peak.py` | Neutral trace construction, joint inverse fit, model comparison, SSA guardrail |
| `cli/full_peak.py` | `igc-full-peak` CLI and all file outputs |

## 1. Trace dataset

Each neutral trace becomes calibrated long-form rows. The detector responds to
molar flow of probe, so the trace is converted to an outlet concentration by
**conserving injected mass**:

```
n_dot(t) = n_injected * S_corr(t) / ∫ S_corr dt        [mol/min]
c_out(t) = n_dot(t) / F_column                          [mol/m³]
P/P0(t)  = c_out(t) * R * T / P_sat(T)                  [-]
```

`n_injected` comes from the declared calibration model and exact `area_unit`;
`F_column` is the declared measured column flow, or an explicit conversion
from declared standard flow and temperature; `P_sat` is the source-attributed
per-injection value in the component table.

This normalisation makes the recovered *trace shape* independent of a constant
FID gain when `n_injected` is held fixed. It does **not** by itself prove
end-to-end gain invariance because `n_injected` is obtained from the calibrated
peak area. `FIDGain` is therefore retained as QC metadata; differing gains need
a gain-specific calibration audit rather than an assumed cancellation.

**Negative baseline residuals.** Raw negative excursions are *retained* in the
trace table (columns `signal_corrected_uV`, for diagnostics) but the physical
quantities (`molar_flow`, `c_out`, `pp0`) are computed from a
nonnegative-clipped copy. Clipping is applied before normalisation so that mass
conservation refers to the same nonnegative signal used physically. Documented
in `build_trace_dataset_from_neutral`.

## 2. Methane transport characterisation

Methane is unretained, so its peak is the system response with `dq/dc = 0`.

### Bracketed dead-time assignment

The CLI supports `fixed_block_mean`, `bracket_interpolated`, `bracket_pre`, and
`bracket_post` transport modes. `bracket_interpolated` uses the contract's
explicit `sequence_index`, identifies dead-time markers before and after the
probe block, and assigns each probe an injection-specific `t0`. If trustworthy
acquisition timestamps are unavailable, interpolation uses sequence position
between the mean positions of the pre- and post-marker groups. The assignment,
interpolation fraction, bracket drift, and neutral marker IDs are exported in
`methane_bracket_assignments.csv`.

The effective plate number and inlet width remain block-level methane-derived
quantities. Only `t0` is interpolated; the assumption is explicit and is tested
against pre-only, post-only, and uniform-shift sensitivity scenarios.
Per **block** (never pooled across blocks — the blocks differ measurably) we
estimate:

- `t0` — void time (first moment / fitted);
- `t_inj` — effective rectangular inlet pulse width (injection + extra-column);
- `N` — apparent plate number lumping axial dispersion and all non-adsorptive
  broadening.

`t0` and `N` are fitted to all methane peaks in the block with the same forward
solver used for the probe, with no adsorption. Reproducibility (SD across
markers) and the between-block shift are reported. These parameters are then
**fixed** when fitting the selected probe; sensitivity to perturbing them is
reported.

Void volume follows without any bed length:

```
V_void = F_column * t0
```

## 3. Forward column model

Equilibrium-dispersive model. Column mass balance over a slice, with `c` in
mol/m³ (gas) and `q` in mol/g (adsorbed on sample mass `m`):

```
∂c/∂t + (m/V_void) ∂q/∂t + u ∂c/∂z = D_ax ∂²c/∂z²
```

The phase ratio is exactly `m/V_void`, both measured — this is why no bed
length is required. With `∂q/∂t = (dq/dc) ∂c/∂t` and

```
k'(c) = (m/V_void) * dq/dc(c)          [dimensionless]
```

nondimensionalising by `ζ = z/L`, `τ = t/t0` (`t0 = L/u`) and using
`N = uL/(2 D_ax)` gives the solved form, in which `L` and `u` cancel:

```
[1 + k'(c)] ∂c/∂τ + ∂c/∂ζ = (1/(2N)) ∂²c/∂ζ²
```

- Inlet BC: `c(0, τ)` = rectangular pulse of width `t_inj/t0` carrying the
  injected moles (Danckwerts-type inlet; the pulse is the measured dose).
- Outlet BC: zero-gradient.
- Discretisation: `N_z` uniform cells, first-order upwind advection (physically
  correct sign, no spurious oscillation) plus central diffusion; explicit
  time-marching with a CFL-limited step.
- Concentrations are clipped to `≥ 0` each step.
- Mass balance (out/in) is returned as a diagnostic on every solve.

`L` and `u` never appear separately — only `t0` (measured) and `N` (fitted on
methane). Physical parameters (`t0`, `V_void`, `m`, `F`, `T`) are kept distinct
from effective nuisance parameters (`N`, `t_inj`) in the outputs.

## 4. Adsorption models

| Model | `q(c)` | `dq/dc` | Finite capacity? |
|---|---|---|---|
| `none` | 0 | 0 | n/a (transport-only control) |
| `henry` | `K_H c` | `K_H` | No |
| `langmuir` | `q_s K c/(1+Kc)` | `q_s K/(1+Kc)²` | **Yes** (`q_s`) |
| `freundlich` | `K_F c^n` | `n K_F c^(n-1)` | No |

**Freundlich exponent convention (important).** `q = K_F c^n` with

- `n < 1` → classical favourable/concave isotherm, `dq/dc` *decreases* with `c`
  (retention weakens at higher loading; peak elutes earlier);
- `n = 1` → Henry;
- `n > 1` → **convex / concentration-strengthened (cooperative)**: `dq/dc`
  *increases* with `c`, so retention strengthens with loading and the peak
  elutes *later* at higher dose.

`n > 1` is therefore the anti-Langmuir/cooperative branch this dataset probes.

Safeguards: `c` is floored at a small positive `c_floor` before evaluating
`c^(n-1)` so `n < 1` cannot produce an infinite slope at the origin; all models
return finite, nonnegative `dq/dc`.

**No model is assumed to yield a monolayer capacity.** Only `langmuir` has a
structural `q_s`.

## 5. Joint inverse fit

All completed probe injections are fitted **simultaneously** with:

- one shared adsorption parameter set (the column is one packed bed);
- per-injection measured dose, flow and temperature;
- per-block transport (`t0`, `t_inj`, `N`) held fixed from methane;
- **amplitude-normalised residuals** — each injection's residual is divided by
  that injection's own peak scale, so large-dose peaks cannot dominate the
  objective purely through amplitude;
- explicit parameter bounds and multiple random starts (`scipy.least_squares`,
  trust-region reflective);
- convergence status recorded per start; the best objective is retained.

Conditional optimiser uncertainty: parameter covariance from the Jacobian at the optimum
(`σ² (JᵀJ)⁻¹`), giving standard errors and a correlation matrix. A parameter is
declared **unidentifiable** when its relative standard error exceeds a
threshold or the Jacobian is rank-deficient/ill-conditioned; the report then
says so instead of quoting the optimiser's number.

Because chromatographic residuals are autocorrelated, these standard errors do
not represent full scientific uncertainty. Bracket/transport sensitivity and
leave-one-injection-out prediction must be interpreted separately.

## 6. Model comparison

Beyond in-sample R²: RMSE, AICc, BIC, mass-balance error, parameter
identifiability, and **leave-one-dose-out** prediction (refit without one
injection, predict it). A model is only preferred if it wins on held-out
prediction as well as in-sample fit.

## 7. SSA guardrail (mandatory)

`compute_ssa_if_identifiable()` returns a surface area **only if all** hold:

1. the selected model structurally has a finite capacity (`has_finite_capacity`);
2. that capacity parameter is identifiable (finite, positive, relative SE below
   threshold);
3. the measured concentration range actually approaches saturation
   (a meaningful fraction of `q_s` is reached).

Otherwise it returns `None` with an explicit reason. The measured/predicted
`P/P0` ratio is **never** converted into an SSA correction. For a cooperative
response the outputs report the isotherm parameters, `q(c)`, `dq/dc` and a
low-concentration affinity measure, and state that geometric SSA is undetermined.

## References

- Guiochon, Felinger, Shirazi & Katti, *Fundamentals of Preparative and
  Nonlinear Chromatography*, 2nd ed., Academic Press (2006) — equilibrium-
  dispersive model, ch. 6; inverse method, ch. 3.
- Rouquerol, Rouquerol & Sing, *Adsorption by Powders and Porous Solids* (1999).
- Thielmann, *J. Chromatogr. A* **1037** (2004) 115–131 — IGC finite concentration.
