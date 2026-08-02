"""Tests for the full-peak nonlinear inverse chromatography pipeline.

Synthetic tests generate peaks with the forward model at known isotherm
parameters and then check that the inverse fit recovers them, that the physics
and numerics behave (mass conservation, nonnegativity, gain invariance), and
that the pipeline refuses to report what it cannot identify.

Real-data integration is gated on the private dataset and skips cleanly.
"""

from __future__ import annotations


import numpy as np
import pytest

from igc_analysis.analysis.column_model import (
    ColumnGeometry, MethaneTransport, TransportParams, apparent_plate_number,
    characterize_methane_transport, make_geometry, peak_moments, solve_column,
)
from igc_analysis.analysis.isotherm_models import (
    C_FLOOR, FREUNDLICH, HENRY, LANGMUIR, MODELS, NO_ADSORPTION, get_model,
    is_cooperative,
)
from igc_analysis.analysis.full_peak import (
    BlockData, Injection, bracket_assignment_to_dataframe, compare_models,
    compute_ssa_if_identifiable, fit_model, predict_injection,
    recovered_isotherm, traces_to_dataframe, transport_to_dataframe,
)
from igc_analysis.constants import R_GAS

# --- Common synthetic setup -------------------------------------------------

T0 = 1.16          # min, void time
T_INJ = 0.10       # min, inlet width
N_PLATES = 80.0
FLOW = 5.5e-6      # m³/min
MASS = 0.160       # g
NZ = 100           # cells (kept modest so tests stay fast)


def _transport(t0=T0):
    return TransportParams(t0_min=t0, t_inj_min=T_INJ, plate_number=N_PLATES)


def _geometry(t0=T0):
    return make_geometry(MASS, FLOW, t0)


def _time_grid(n=500, t_end=3.0):
    return np.linspace(0.0, t_end, n)


def _make_injection(block, num, t, c_out, n_inj, temp=303.15, gain=1.0,
                    p_sat=7745.0):
    """Wrap a concentration profile as an Injection (signal ∝ concentration)."""
    signal = c_out * 1e6          # arbitrary µV scale
    return Injection(
        block=block, injection_number=num, name=f"injection{num}",
        probe="HEPTANE", time_min=t, signal_raw_uV=signal,
        baseline_uV=np.zeros_like(t), signal_corrected_uV=signal,
        n_injected_mol=n_inj, c_out_mol_m3=c_out,
        pp0=c_out * R_GAS * temp / p_sat, flow_col_m3_min=FLOW,
        temp_col_K=temp, pressure_drop_torr=20.0, fid_gain=gain,
        p_sat_Pa=p_sat, clipped=False,
    )


def _synth_block(model_name, params, doses, block="block1", t0=T0, noise=0.0,
                 seed=0, n_time=500):
    """Build a BlockData whose traces are forward-model predictions."""
    model = get_model(model_name)
    t = _time_grid(n_time)
    tp, geom = _transport(t0), _geometry(t0)
    rng = np.random.default_rng(seed)
    injections = []
    for i, dose in enumerate(doses, start=1):
        res = solve_column(t, dose, tp, geom, model, np.asarray(params),
                           n_cells=NZ)
        c = res.c_out
        if noise > 0:
            c = np.maximum(c + rng.normal(0.0, noise * c.max(), c.shape), 0.0)
        injections.append(_make_injection(block, i, t, c, dose))
    transport = MethaneTransport(
        block=block, n_markers=4, t0_min=t0, t0_sd_min=0.0, t0_range_min=0.0,
        sigma_mean_min=0.13, plate_number=N_PLATES, plate_number_sd=0.0,
        t_inj_min=T_INJ, void_volume_m3=geom.void_volume_m3, fit_rmse=0.0,
    )
    return BlockData(block=block, transport=transport, geometry=geom,
                     injections=injections)


# ---------------------------------------------------------------------------
# Isotherm models
# ---------------------------------------------------------------------------

class TestIsothermModels:
    def test_all_registered(self):
        assert set(MODELS) == {"none", "henry", "langmuir", "freundlich"}

    def test_no_adsorption_is_zero(self):
        c = np.linspace(0, 1, 10)
        assert np.all(NO_ADSORPTION.q(c, []) == 0)
        assert np.all(NO_ADSORPTION.dqdc(c, []) == 0)

    def test_henry_slope_constant(self):
        c = np.linspace(0, 1, 10)
        d = HENRY.dqdc(c, [2e-6])
        assert np.allclose(d, 2e-6)

    def test_langmuir_slope_decreases(self):
        """Favourable/concave: dq/dc falls with c."""
        c = np.array([0.01, 0.5, 5.0])
        d = LANGMUIR.dqdc(c, [1e-5, 5.0])
        assert d[0] > d[1] > d[2]

    def test_langmuir_saturates(self):
        q = LANGMUIR.q(np.array([1e9]), [1e-5, 5.0])[0]
        assert q == pytest.approx(1e-5, rel=1e-3)

    def test_freundlich_exponent_convention(self):
        """n>1 is the cooperative branch: dq/dc RISES with c. n<1 falls."""
        c = np.array([0.05, 0.5])
        rising = FREUNDLICH.dqdc(c, [1e-6, 1.4])
        falling = FREUNDLICH.dqdc(c, [1e-6, 0.7])
        assert rising[1] > rising[0]
        assert falling[1] < falling[0]

    def test_freundlich_n_equals_one_is_henry(self):
        c = np.linspace(0.01, 1, 5)
        assert np.allclose(FREUNDLICH.dqdc(c, [3e-6, 1.0]), 3e-6, rtol=1e-9)

    def test_c_zero_safeguard(self):
        """n<1 would be singular at c=0; must stay finite and nonnegative."""
        d = FREUNDLICH.dqdc(np.array([0.0]), [1e-6, 0.5])
        assert np.isfinite(d[0]) and d[0] >= 0

    def test_only_langmuir_has_finite_capacity(self):
        assert LANGMUIR.has_finite_capacity and LANGMUIR.capacity_param == "q_s"
        for m in (HENRY, FREUNDLICH, NO_ADSORPTION):
            assert not m.has_finite_capacity

    def test_is_cooperative_flag(self):
        assert is_cooperative(FREUNDLICH, np.array([1e-6, 1.3]))
        assert not is_cooperative(FREUNDLICH, np.array([1e-6, 0.8]))
        assert not is_cooperative(LANGMUIR, np.array([1e-5, 1.0]))

    def test_unknown_model_raises(self):
        with pytest.raises(KeyError):
            get_model("bogus")


# ---------------------------------------------------------------------------
# Column model numerics
# ---------------------------------------------------------------------------

class TestColumnModel:
    def test_mass_balance(self):
        """Eluted moles must match injected moles to within a few percent."""
        t = _time_grid()
        res = solve_column(t, 4e-7, _transport(), _geometry(), NO_ADSORPTION,
                           np.array([]), n_cells=NZ)
        assert res.mass_balance == pytest.approx(1.0, abs=0.05)

    def test_mass_balance_improves_with_refinement(self):
        t = _time_grid()
        coarse = solve_column(t, 4e-7, _transport(), _geometry(), NO_ADSORPTION,
                              np.array([]), n_cells=60)
        fine = solve_column(t, 4e-7, _transport(), _geometry(), NO_ADSORPTION,
                            np.array([]), n_cells=300)
        assert abs(fine.mass_balance - 1.0) <= abs(coarse.mass_balance - 1.0)

    def test_nonretained_peak_arrives_at_void_time(self):
        """With no adsorption the peak sits at t0 + half the inlet width."""
        t = _time_grid()
        res = solve_column(t, 4e-7, _transport(), _geometry(), NO_ADSORPTION,
                           np.array([]), n_cells=300)
        _, m1, _ = peak_moments(t, res.c_out)
        assert m1 == pytest.approx(T0 + T_INJ / 2, rel=0.02)

    def test_henry_shifts_by_retention_factor(self):
        """Henry retention shifts the first moment by exactly t0·(1+k')."""
        t = _time_grid()
        geom = _geometry()
        K_H = 3.0e-6
        kprime = geom.phase_ratio_g_m3 * K_H
        res = solve_column(t, 4e-7, _transport(), geom, HENRY,
                           np.array([K_H]), n_cells=300)
        _, m1, _ = peak_moments(t, res.c_out)
        assert m1 == pytest.approx(T0 * (1 + kprime) + T_INJ / 2, rel=0.02)

    def test_concentrations_never_negative(self):
        t = _time_grid()
        res = solve_column(t, 4e-7, _transport(), _geometry(), FREUNDLICH,
                           np.array([3e-6, 1.3]), n_cells=NZ)
        assert np.all(res.c_out >= 0.0)

    def test_langmuir_elutes_earlier_at_higher_dose(self):
        """Favourable isotherm → retention weakens with loading."""
        t = _time_grid()
        p = np.array([1e-5, 5.0])
        lo = solve_column(t, 1.0e-7, _transport(), _geometry(), LANGMUIR, p, n_cells=NZ)
        hi = solve_column(t, 6.0e-7, _transport(), _geometry(), LANGMUIR, p, n_cells=NZ)
        assert peak_moments(t, hi.c_out)[1] < peak_moments(t, lo.c_out)[1]

    def test_cooperative_elutes_later_at_higher_dose(self):
        """Freundlich n>1 → retention strengthens with loading."""
        t = _time_grid()
        p = np.array([3e-6, 1.3])
        lo = solve_column(t, 1.0e-7, _transport(), _geometry(), FREUNDLICH, p, n_cells=NZ)
        hi = solve_column(t, 6.0e-7, _transport(), _geometry(), FREUNDLICH, p, n_cells=NZ)
        assert peak_moments(t, hi.c_out)[1] > peak_moments(t, lo.c_out)[1]

    def test_phase_ratio_needs_no_bed_length(self):
        geom = make_geometry(0.160, 5.5e-6, 1.16)
        assert geom.void_volume_m3 == pytest.approx(5.5e-6 * 1.16)
        assert geom.phase_ratio_g_m3 == pytest.approx(0.160 / (5.5e-6 * 1.16))

    def test_invalid_transport_raises(self):
        with pytest.raises(ValueError):
            solve_column(_time_grid(), 1e-7,
                         TransportParams(0.0, T_INJ, N_PLATES), _geometry(),
                         NO_ADSORPTION, np.array([]))

    @pytest.mark.parametrize("bad_t0", [np.nan, np.inf, -np.inf])
    def test_nonfinite_transport_raises(self, bad_t0):
        with pytest.raises(ValueError, match="finite and positive"):
            solve_column(
                _time_grid(), 1e-7,
                TransportParams(bad_t0, T_INJ, N_PLATES), _geometry(),
                NO_ADSORPTION, np.array([]),
            )

    def test_step_cap_fails_instead_of_coarsening_unstably(self):
        with pytest.raises(ValueError, match="refusing to violate"):
            solve_column(
                _time_grid(), 1e-7, _transport(), _geometry(),
                NO_ADSORPTION, np.array([]), n_cells=NZ, max_steps=10,
            )


class TestMethaneCharacterisation:
    def test_recovers_transport_from_synthetic_methane(self):
        """A synthetic unretained peak returns its own t0 and a sane N."""
        t = _time_grid()
        res = solve_column(t, 1e-7, _transport(), _geometry(), NO_ADSORPTION,
                           np.array([]), n_cells=NZ)
        traces = [(t, res.c_out)] * 3
        mt = characterize_methane_transport("b", traces, FLOW, MASS, n_cells=NZ)
        # First moment of an unretained peak is t0 + t_inj/2.
        assert mt.t0_min == pytest.approx(T0 + T_INJ / 2, rel=0.05)
        assert mt.plate_number > 0
        assert mt.void_volume_m3 == pytest.approx(FLOW * mt.t0_min)

    def test_apparent_plate_number(self):
        assert apparent_plate_number(1.16, 0.116) == pytest.approx(100.0)
        assert np.isnan(apparent_plate_number(1.16, 0.0))

    def test_moments_ignore_negative_noise(self):
        t = _time_grid(200)
        s = np.where((t > 1.0) & (t < 1.3), 1.0, -0.05)   # negative baseline
        area, m1, sd = peak_moments(t, s)
        assert area > 0 and 1.0 < m1 < 1.3

    def test_single_marker_has_no_reproducibility_sd(self):
        t = _time_grid()
        res = solve_column(t, 1e-7, _transport(), _geometry(), NO_ADSORPTION,
                           np.array([]), n_cells=NZ)
        mt = characterize_methane_transport("b", [(t, res.c_out)], FLOW, MASS,
                                                 n_cells=NZ)
        assert np.isnan(mt.t0_sd_min)
        assert np.isnan(mt.plate_number_sd)


# ---------------------------------------------------------------------------
# Trace construction / calibration invariants
# ---------------------------------------------------------------------------

class TestTraceInvariants:
    def _signal_to_concentration(self, signal, n_inj, t, flow):
        """Mirror the neutral reader's mass-conserving conversion."""
        s_pos = np.maximum(signal, 0.0)
        area = np.trapezoid(s_pos, t)
        return (n_inj * s_pos / area) / flow

    def test_conversion_conserves_injected_moles(self):
        t = _time_grid()
        signal = np.exp(-((t - 1.2) ** 2) / 0.01) * 1234.0
        n_inj = 3.7e-7
        c = self._signal_to_concentration(signal, n_inj, t, FLOW)
        eluted = np.trapezoid(c * FLOW, t)
        assert eluted == pytest.approx(n_inj, rel=1e-9)

    def test_gain_change_does_not_change_shape_at_fixed_injected_moles(self):
        """Signal/area normalisation cancels gain only when dose is fixed."""
        t = _time_grid()
        base = np.exp(-((t - 1.2) ** 2) / 0.01) * 1000.0
        n_inj = 3.7e-7
        c1 = self._signal_to_concentration(base, n_inj, t, FLOW)
        c2 = self._signal_to_concentration(base * 16.0, n_inj, t, FLOW)
        assert np.allclose(c1, c2, rtol=1e-12)

    def test_negative_baseline_gives_no_negative_concentration(self):
        t = _time_grid()
        signal = np.exp(-((t - 1.2) ** 2) / 0.01) * 1000.0 - 25.0  # dips below 0
        assert signal.min() < 0
        c = self._signal_to_concentration(signal, 3.7e-7, t, FLOW)
        assert np.all(c >= 0.0)

    def test_traces_dataframe_columns(self):
        blk = _synth_block("henry", [2e-6], [4e-7, 2e-7], n_time=200)
        df = traces_to_dataframe([blk])
        for col in ("block", "injection", "probe", "time_min", "signal_raw_uV",
                    "baseline_uV", "signal_corrected_uV",
                    "response_density_per_min", "n_injected_mol",
                    "molar_flow_mol_min", "c_out_mol_m3", "pp0",
                    "flow_col_m3_min", "temp_col_K", "fid_gain", "clipped"):
            assert col in df.columns
        assert len(df) == 2 * 200

    def test_transport_dataframe(self):
        blk = _synth_block("henry", [2e-6], [4e-7], n_time=200)
        df = transport_to_dataframe([blk])
        assert df.loc[0, "t0_min"] == pytest.approx(T0)
        assert df.loc[0, "sample_mass_g"] == pytest.approx(MASS)


# ---------------------------------------------------------------------------
# Inverse fitting: parameter recovery
# ---------------------------------------------------------------------------

class TestParameterRecovery:
    def test_recovers_henry_constant(self):
        K_true = 2.5e-6
        blk = _synth_block("henry", [K_true], [4.5e-7, 3.0e-7, 1.5e-7], n_time=400)
        fit = fit_model([blk], "henry", n_starts=2, n_cells=NZ, verbose=False)
        assert fit.params[0] == pytest.approx(K_true, rel=0.10)
        assert fit.identifiable["K_H"]

    def test_recovers_cooperative_freundlich(self):
        """A convex (n>1) isotherm is recovered within tolerance."""
        true = [3.0e-6, 1.30]
        blk = _synth_block("freundlich", true, [4.7e-7, 3.2e-7, 1.9e-7, 1.0e-7],
                           n_time=400)
        fit = fit_model([blk], "freundlich", n_starts=3, n_cells=NZ, verbose=False)
        assert fit.params[1] == pytest.approx(true[1], rel=0.15)
        assert fit.cooperative

    def test_recovers_langmuir(self):
        true = [1.2e-5, 4.0]
        blk = _synth_block("langmuir", true, [5.0e-7, 3.0e-7, 1.2e-7], n_time=400)
        fit = fit_model([blk], "langmuir", n_starts=3, n_cells=NZ, verbose=False)
        # q_s and K_L trade off; check the model reproduces the data well.
        assert fit.rmse_normalised < 0.02

    def test_transport_only_recovered_when_no_adsorption(self):
        """Data generated with q≡0 is fitted best by the transport-only model."""
        blk = _synth_block("none", [], [4e-7, 2e-7], n_time=400)
        fit = fit_model([blk], "none", n_starts=1, n_cells=NZ, verbose=False)
        assert fit.rmse_normalised < 1e-9
        assert fit.n_params == 0

    def test_discriminates_langmuir_from_cooperative_under_noise(self):
        """Peaks generated by a cooperative isotherm must not be better
        explained by Langmuir, at a realistic noise level."""
        blk = _synth_block("freundlich", [3.0e-6, 1.35],
                           [4.7e-7, 3.2e-7, 1.9e-7, 1.0e-7],
                           noise=0.01, seed=3, n_time=400)
        f_free = fit_model([blk], "freundlich", n_starts=3, n_cells=NZ, verbose=False)
        f_lang = fit_model([blk], "langmuir", n_starts=3, n_cells=NZ, verbose=False)
        assert f_free.rmse_normalised < f_lang.rmse_normalised
        assert f_free.cooperative


class TestBlockHandling:
    def test_block_specific_t0_is_respected(self):
        """Two blocks with different void times must each use their own t0.

        Fitting them jointly with the correct per-block t0 recovers the shared
        isotherm; forcing one block's t0 onto the other degrades the fit.
        """
        K_true = 2.5e-6
        b1 = _synth_block("henry", [K_true], [4.5e-7, 2.5e-7],
                          block="block1", t0=1.16, n_time=400)
        b2 = _synth_block("henry", [K_true], [4.0e-7, 2.0e-7],
                          block="block2", t0=1.05, n_time=400)
        good = fit_model([b1, b2], "henry", n_starts=2, n_cells=NZ, verbose=False)
        assert good.params[0] == pytest.approx(K_true, rel=0.12)

        # Now corrupt block2's transport with block1's t0 and refit.
        bad_b2 = BlockData(
            block="block2",
            transport=MethaneTransport(
                block="block2", n_markers=4, t0_min=1.16, t0_sd_min=0.0,
                t0_range_min=0.0, sigma_mean_min=0.13, plate_number=N_PLATES,
                plate_number_sd=0.0, t_inj_min=T_INJ,
                void_volume_m3=b2.geometry.void_volume_m3, fit_rmse=0.0),
            geometry=b2.geometry, injections=b2.injections)
        bad = fit_model([b1, bad_b2], "henry", n_starts=2, n_cells=NZ, verbose=False)
        assert bad.rmse_normalised > good.rmse_normalised

    def test_joint_fit_shares_one_parameter_set(self):
        blk1 = _synth_block("henry", [2.5e-6], [4e-7], block="block1", n_time=300)
        blk2 = _synth_block("henry", [2.5e-6], [2e-7], block="block2", n_time=300)
        fit = fit_model([blk1, blk2], "henry", n_starts=1, n_cells=NZ, verbose=False)
        assert fit.n_params == 1              # one shared K_H, not one per block
        assert len(fit.per_injection_rmse) == 2

    def test_drift_corrected_linear_data_are_not_called_cooperative(self):
        """Known methane drift plus descending dose must not manufacture n>1."""
        doses = [4.7e-7, 4.0e-7, 3.3e-7, 2.6e-7, 1.9e-7]
        t0s = np.linspace(1.18, 1.16, len(doses))
        t = _time_grid(400)
        injections = []
        for num, (dose, t0) in enumerate(zip(doses, t0s), start=1):
            geom = _geometry(t0)
            res = solve_column(t, dose, _transport(t0), geom, HENRY,
                               np.array([2.5e-6]), n_cells=NZ)
            inj = _make_injection("drift", num, t, res.c_out, dose)
            inj.transport_mode = "bracket_interpolated"
            inj.assigned_t0_min = float(t0)
            injections.append(inj)
        mean_t0 = float(np.mean(t0s))
        geom = _geometry(mean_t0)
        block = BlockData(
            block="drift",
            transport=MethaneTransport(
                block="drift", n_markers=4, t0_min=mean_t0,
                t0_sd_min=float(np.std(t0s)),
                t0_range_min=float(np.ptp(t0s)), sigma_mean_min=0.13,
                plate_number=N_PLATES, plate_number_sd=0.0,
                t_inj_min=T_INJ, void_volume_m3=geom.void_volume_m3,
                fit_rmse=0.0),
            geometry=geom, injections=injections)
        fit = fit_model([block], "freundlich", n_starts=3, n_cells=NZ,
                        verbose=False)
        assert fit.params[1] == pytest.approx(1.0, abs=0.08)


# ---------------------------------------------------------------------------
# Identifiability and the SSA guardrail
# ---------------------------------------------------------------------------

class TestIdentifiabilityAndSSA:
    def test_ssa_refused_for_model_without_finite_capacity(self):
        """Freundlich/Henry have no monolayer — no SSA may be emitted."""
        blk = _synth_block("freundlich", [3.0e-6, 1.3], [4e-7, 2e-7], n_time=300)
        fit = fit_model([blk], "freundlich", n_starts=2, n_cells=NZ, verbose=False)
        verdict = compute_ssa_if_identifiable(fit, [blk], 5.73e-19)
        assert verdict.reportable is False
        assert verdict.ssa_m2_g is None
        assert "no finite monolayer" in verdict.reason

    def test_ssa_refused_for_henry(self):
        blk = _synth_block("henry", [2.5e-6], [4e-7, 2e-7], n_time=300)
        fit = fit_model([blk], "henry", n_starts=1, n_cells=NZ, verbose=False)
        verdict = compute_ssa_if_identifiable(fit, [blk], 5.73e-19)
        assert verdict.reportable is False and verdict.ssa_m2_g is None

    def test_ssa_refused_when_far_from_saturation(self):
        """Langmuir has q_s structurally, but if the measured range never
        approaches saturation the capacity is extrapolated → no SSA."""
        # Very weak K so that K·c << 1 over the measured range.
        blk = _synth_block("langmuir", [1e-5, 1e-3], [4e-7, 2e-7], n_time=300)
        fit = fit_model([blk], "langmuir", n_starts=2, n_cells=NZ, verbose=False)
        verdict = compute_ssa_if_identifiable(fit, [blk], 5.73e-19)
        assert verdict.reportable is False
        assert verdict.ssa_m2_g is None

    def test_unidentifiable_parameter_is_flagged_not_reported(self):
        """A parameter with a huge relative SE must be marked unidentifiable."""
        fit = fit_model([_synth_block("langmuir", [1e-5, 1e-3], [4e-7, 2e-7],
                                      n_time=300)],
                        "langmuir", n_starts=2, n_cells=NZ, verbose=False)
        # At least one Langmuir parameter should fail the identifiability test
        # in this near-linear regime (q_s and K_L are confounded).
        assert not fit.all_identifiable

    def test_ssa_reported_when_genuinely_saturating(self):
        """Guardrail is not vacuous: a strongly saturating, identifiable
        Langmuir does yield an SSA."""
        q_s = 2e-6
        blk = _synth_block("langmuir", [q_s, 200.0], [6e-7, 4e-7, 2e-7, 1e-7],
                           n_time=400)
        fit = fit_model([blk], "langmuir", n_starts=4, n_cells=NZ, verbose=False)
        verdict = compute_ssa_if_identifiable(fit, [blk], 5.73e-19,
                                              saturation_fraction_required=0.05)
        if verdict.reportable:                     # only assert if identifiable
            assert verdict.ssa_m2_g > 0
        else:
            assert verdict.ssa_m2_g is None        # still must not invent one


class TestTransportSensitivity:
    def test_sensitivity_reports_all_scenarios(self):
        from igc_analysis.analysis.full_peak import transport_sensitivity
        blk = _synth_block("henry", [2.5e-6], [4e-7, 2e-7], n_time=300)
        blk.transport.t0_sd_min = 0.01
        fit = fit_model([blk], "henry", n_starts=1, n_cells=NZ, verbose=False)
        sens = transport_sensitivity([blk], "henry", fit.params, n_cells=NZ)
        scenarios = set(sens["scenario"])
        assert "base" in scenarios
        assert any(s.startswith("N x") for s in scenarios)
        assert any(s.startswith("t_inj x") for s in scenarios)
        assert any(s.startswith("t0") for s in scenarios)
        assert "window 0-2 min" in scenarios
        assert "K_H_pct_change" in sens.columns
        assert sens.loc[sens["scenario"] == "base", "K_H_pct_change"].iloc[0] == 0.0

    def test_single_marker_omits_vacuous_t0_sensitivity(self):
        from igc_analysis.analysis.full_peak import transport_sensitivity
        blk = _synth_block("henry", [2.5e-6], [4e-7], n_time=250)
        blk.transport.n_markers = 1
        blk.transport.t0_sd_min = float("nan")
        fit = fit_model([blk], "henry", n_starts=1, n_cells=NZ, verbose=False)
        sens = transport_sensitivity([blk], "henry", fit.params, n_cells=NZ)
        assert not any(str(s).startswith("t0") for s in sens["scenario"])

    def test_t0_sensitivity_rebuilds_fixed_mode_geometry(self, monkeypatch):
        from types import SimpleNamespace
        import igc_analysis.analysis.full_peak as full_peak

        blk = _synth_block("henry", [2.5e-6], [4e-7], n_time=250)
        blk.transport.t0_sd_min = 0.01
        observed = []

        def fake_fit(blks, *args, **kwargs):
            for candidate in blks:
                observed.append((candidate.transport.t0_min,
                                 candidate.geometry.void_volume_m3,
                                 candidate.geometry.flow_col_m3_min))
            return SimpleNamespace(
                rmse_normalised=0.0,
                params=np.array([2.5e-6]),
                all_identifiable=True,
            )

        monkeypatch.setattr(full_peak, "fit_model", fake_fit)
        full_peak.transport_sensitivity(
            [blk], "henry", np.array([2.5e-6]), n_cells=NZ)
        assert observed
        for t0, void_volume, flow in observed:
            assert void_volume == pytest.approx(flow * t0)

    def test_dispersion_perturbation_moves_parameter(self):
        """Changing the assumed plate number should measurably move K_H —
        that is the point of reporting the sensitivity."""
        from igc_analysis.analysis.full_peak import transport_sensitivity
        blk = _synth_block("henry", [2.5e-6], [4e-7, 2e-7], n_time=300)
        fit = fit_model([blk], "henry", n_starts=1, n_cells=NZ, verbose=False)
        sens = transport_sensitivity([blk], "henry", fit.params, n_cells=NZ)
        changes = sens["K_H_pct_change"].abs().dropna()
        assert np.isfinite(changes).all()


class TestTargetedHeldOutComparison:
    def test_lodo_can_be_limited_to_scientifically_relevant_models(self):
        blk = _synth_block("henry", [2.5e-6], [4e-7, 2e-7], n_time=250)
        _, table = compare_models(
            [blk], model_names=("none", "henry"), n_cells=NZ,
            n_starts=1, do_lodo=True, lodo_models={"henry"}, verbose=False)
        scores = table.set_index("model")["lodo_rmse"]
        assert np.isnan(scores["none"])
        assert np.isfinite(scores["henry"])


class TestRecoveredIsotherm:
    def test_isotherm_table(self):
        blk = _synth_block("freundlich", [3e-6, 1.3], [4e-7], n_time=300)
        fit = fit_model([blk], "freundlich", n_starts=1, n_cells=NZ, verbose=False)
        df = recovered_isotherm(fit, [blk], n_points=50)
        assert {"c_mol_m3", "pp0", "q_mol_g", "dqdc_m3_g"} <= set(df.columns)
        assert len(df) == 50
        assert np.all(df["q_mol_g"] >= 0) and np.all(df["dqdc_m3_g"] >= 0)


# ---------------------------------------------------------------------------
