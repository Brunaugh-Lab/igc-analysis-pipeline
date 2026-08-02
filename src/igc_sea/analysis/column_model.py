"""Equilibrium-dispersive column model for finite-concentration IGC.

Solves the nonlinear chromatography mass balance forward in time so that a
complete peak shape can be predicted from an assumed adsorption isotherm, and
therefore *inverted* to identify that isotherm from measured peaks.

Model
-----
Over a column slice, with ``c`` the gas-phase concentration (mol/m³) and ``q``
the adsorbed amount per gram of sample (mol/g):

.. math::

    \\frac{\\partial c}{\\partial t}
    + \\frac{m}{V_{void}} \\frac{\\partial q}{\\partial t}
    + u \\frac{\\partial c}{\\partial z}
    = D_{ax} \\frac{\\partial^2 c}{\\partial z^2}

The phase ratio is exactly ``m / V_void`` — both measured — which is why this
implementation needs **no packed-bed length**.  Writing
``∂q/∂t = (dq/dc) ∂c/∂t`` and defining the dimensionless local retention factor

.. math::

    k'(c) = \\frac{m}{V_{void}} \\frac{dq}{dc}(c)

then nondimensionalising with ``ζ = z/L`` and ``τ = t/t_0`` (``t_0 = L/u``),
and using the apparent plate number ``N = uL/(2 D_ax)``:

.. math::

    [1 + k'(c)] \\frac{\\partial c}{\\partial \\tau}
    + \\frac{\\partial c}{\\partial \\zeta}
    = \\frac{1}{2N} \\frac{\\partial^2 c}{\\partial \\zeta^2}

``L`` and ``u`` cancel; only the **measured** void time ``t_0`` and the
**effective** plate number ``N`` survive.  ``N`` is a nuisance parameter that
lumps axial dispersion, injection profile and extra-column broadening; it is
calibrated on methane and then held fixed.

Numerics
--------
Method of lines on ``N_z`` uniform cells: first-order upwind for advection
(monotone, no spurious oscillation at sharp fronts) and central differences for
dispersion, marched explicitly with a CFL-limited step.  Concentrations are
clipped at zero every step, so a noisy or negative inlet can never produce a
negative physical concentration.  Every solve returns a mass-balance
diagnostic (eluted moles / injected moles).

References
----------
Guiochon, Felinger, Shirazi & Katti, *Fundamentals of Preparative and Nonlinear
Chromatography*, 2nd ed. (2006), ch. 6 (equilibrium-dispersive model).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from igc_sea.analysis.isotherm_models import IsothermModel


@dataclass
class TransportParams:
    """Non-adsorptive transport parameters for one block.

    Attributes
    ----------
    t0_min : float
        Void (dead) time, minutes — **physical**, from methane.
    t_inj_min : float
        Effective rectangular inlet pulse width, minutes — *effective*
        nuisance parameter (injection + extra-column contribution).
    plate_number : float
        Apparent plate number N — *effective* nuisance parameter lumping axial
        dispersion and all non-adsorptive broadening.
    """

    t0_min: float
    t_inj_min: float
    plate_number: float


@dataclass
class ColumnGeometry:
    """Measured column quantities, none of which require a bed length.

    Attributes
    ----------
    sample_mass_g : float
        Adsorbent mass (g).
    flow_col_m3_min : float
        Carrier volumetric flow at column temperature (m³/min).
    void_volume_m3 : float
        ``F_column × t0`` (m³).
    """

    sample_mass_g: float
    flow_col_m3_min: float
    void_volume_m3: float

    @property
    def phase_ratio_g_m3(self) -> float:
        """m / V_void in g/m³ — the coefficient multiplying dq/dc."""
        if self.void_volume_m3 <= 0:
            return 0.0
        return self.sample_mass_g / self.void_volume_m3


@dataclass
class SolveResult:
    """Outcome of a forward solve."""

    time_min: np.ndarray        # outlet time grid (min)
    c_out: np.ndarray           # outlet concentration (mol/m³)
    mass_balance: float         # eluted / injected moles (1.0 = perfect)
    n_steps: int


def make_geometry(sample_mass_g: float, flow_col_m3_min: float,
                  t0_min: float) -> ColumnGeometry:
    """Build :class:`ColumnGeometry` from measured flow and methane void time."""
    return ColumnGeometry(
        sample_mass_g=sample_mass_g,
        flow_col_m3_min=flow_col_m3_min,
        void_volume_m3=flow_col_m3_min * t0_min,
    )


def solve_column(
    time_min: np.ndarray,
    n_injected_mol: float,
    transport: TransportParams,
    geometry: ColumnGeometry,
    model: IsothermModel,
    params: np.ndarray,
    n_cells: int = 120,
    cfl: float = 0.4,
    max_steps: int = 400_000,
) -> SolveResult:
    """Predict the outlet concentration profile for one injection.

    Parameters
    ----------
    time_min : np.ndarray
        Output time grid (minutes), ascending; the measured trace grid.
    n_injected_mol : float
        Moles injected from the declared calibration.
    transport : TransportParams
        Block transport parameters (``t0``, ``t_inj``, ``N``).
    geometry : ColumnGeometry
        Measured mass, flow and void volume.
    model : IsothermModel
        Adsorption isotherm supplying ``dq/dc``.
    params : np.ndarray
        Isotherm parameters.
    n_cells : int
        Spatial discretisation of ζ ∈ [0, 1].
    cfl : float
        Courant safety factor for the explicit step (< 1).
    max_steps : int
        Hard cap on time steps (guards against a pathological parameter set).

    Returns
    -------
    SolveResult
        Outlet concentration interpolated onto ``time_min``, plus the
        mass-balance diagnostic.
    """
    time = np.asarray(time_min, dtype=float)
    parameters = np.asarray(params, dtype=float)
    t0 = float(transport.t0_min)
    N = float(transport.plate_number)
    F = float(geometry.flow_col_m3_min)
    phase_ratio = float(geometry.phase_ratio_g_m3)
    dose = float(n_injected_mol)

    if time.ndim != 1 or time.size < 2:
        raise ValueError("time_min must be a one-dimensional grid with at least two points")
    if not np.all(np.isfinite(time)) or time[0] < 0 or np.any(np.diff(time) <= 0):
        raise ValueError("time_min must contain finite, nonnegative, strictly increasing values")
    if not np.isfinite(t0) or not np.isfinite(N) or t0 <= 0 or N <= 0:
        raise ValueError("t0 and plate_number must be finite and positive")
    if not np.isfinite(F) or F <= 0:
        raise ValueError("flow_col_m3_min must be finite and positive")
    if not np.isfinite(phase_ratio) or phase_ratio < 0:
        raise ValueError("phase ratio must be finite and nonnegative")
    if not np.isfinite(dose) or dose < 0:
        raise ValueError("n_injected_mol must be finite and nonnegative")
    if not np.all(np.isfinite(parameters)):
        raise ValueError("isotherm parameters must be finite")
    if not isinstance(n_cells, (int, np.integer)) or n_cells < 2:
        raise ValueError("n_cells must be an integer of at least 2")
    if not np.isfinite(cfl) or not 0 < cfl < 1:
        raise ValueError("cfl must be finite and between 0 and 1")
    if not isinstance(max_steps, (int, np.integer)) or max_steps < 1:
        raise ValueError("max_steps must be a positive integer")

    # --- Inlet pulse (dimensionless time) ---------------------------------
    # A rectangular pulse of width t_inj carrying n_injected moles at flow F:
    #   c_in = n_injected / (F * t_inj)
    declared_t_inj = float(transport.t_inj_min)
    if not np.isfinite(declared_t_inj) or declared_t_inj < 0:
        raise ValueError("t_inj_min must be finite and nonnegative")
    t_inj = max(declared_t_inj, 1e-6)
    tau_inj = t_inj / t0
    c_in = dose / (F * t_inj)
    if not np.isfinite(c_in):
        raise ValueError("inlet concentration is non-finite")

    # --- Grid --------------------------------------------------------------
    nz = int(n_cells)
    dzeta = 1.0 / nz
    D = 1.0 / (2.0 * N)          # dimensionless dispersion coefficient

    # Explicit stability: advection (speed 1) and diffusion.
    dtau_adv = cfl * dzeta
    dtau_dif = cfl * 0.5 * dzeta * dzeta / D if D > 0 else np.inf
    dtau = float(min(dtau_adv, dtau_dif))

    tau_end = float(time[-1]) / t0
    n_steps = int(np.ceil(tau_end / dtau))
    if n_steps > max_steps:
        raise ValueError(
            f"stable explicit solve requires {n_steps} steps, exceeding "
            f"max_steps={max_steps}; refusing to violate the stability limit"
        )

    c = np.zeros(nz + 1)         # cell-centred + inlet ghost at index 0
    out_tau = np.empty(n_steps + 1)
    out_c = np.empty(n_steps + 1)
    out_tau[0] = 0.0
    out_c[0] = 0.0

    eluted = 0.0                 # dimensionless integral of outlet c dτ

    for step in range(1, n_steps + 1):
        tau = step * dtau

        # Inlet boundary: rectangular pulse.
        c[0] = c_in if tau <= tau_inj else 0.0

        # Local retention factor from the isotherm slope.
        kprime = phase_ratio * np.asarray(model.dqdc(c[1:], parameters), dtype=float)
        denom = 1.0 + kprime
        if not np.all(np.isfinite(denom)) or np.any(denom <= 0):
            raise ValueError("isotherm produced a non-finite or non-positive retardation factor")

        # Upwind advection: (c_i - c_{i-1}) / dzeta
        adv = (c[1:] - c[:-1]) / dzeta

        # Central diffusion with zero-gradient outlet.
        c_ext = np.empty(nz + 2)
        c_ext[:nz + 1] = c
        c_ext[nz + 1] = c[nz]                      # ∂c/∂ζ = 0 at outlet
        dif = (c_ext[2:] - 2.0 * c_ext[1:nz + 1] + c_ext[:nz]) / (dzeta * dzeta)

        c[1:] = c[1:] + dtau * (-adv + D * dif) / denom
        if not np.all(np.isfinite(c)):
            raise FloatingPointError("column solve produced non-finite concentrations")
        np.maximum(c, 0.0, out=c)                  # no negative concentrations

        out_tau[step] = tau
        out_c[step] = c[nz]
        eluted += c[nz] * dtau

    # Mass balance: ∫ c_out F dt = F t0 ∫ c_out dτ  should equal n_injected.
    eluted_mol = eluted * F * t0
    mass_balance = eluted_mol / dose if dose > 0 else 0.0

    c_interp = np.interp(time / t0, out_tau, out_c,
                         left=0.0, right=0.0)
    if not np.all(np.isfinite(c_interp)) or not np.isfinite(mass_balance):
        raise FloatingPointError("column solve returned a non-finite result")
    return SolveResult(time_min=time, c_out=c_interp,
                       mass_balance=float(mass_balance), n_steps=n_steps)


# ---------------------------------------------------------------------------
# Methane-based transport characterisation
# ---------------------------------------------------------------------------

def peak_moments(time_min: np.ndarray, signal: np.ndarray) -> tuple[float, float, float]:
    """Return ``(area, first_moment, sigma)`` of a nonnegative trace.

    The first moment is the centre of mass; ``sigma`` is the square root of the
    second central moment.  Both are computed on the nonnegative part of the
    signal, so baseline noise cannot bias them negative.
    """
    s = np.maximum(np.asarray(signal, dtype=float), 0.0)
    t = np.asarray(time_min, dtype=float)
    area = float(np.trapezoid(s, t))
    if area <= 0:
        return 0.0, float("nan"), float("nan")
    m1 = float(np.trapezoid(s * t, t) / area)
    m2 = float(np.trapezoid(s * (t - m1) ** 2, t) / area)
    return area, m1, float(np.sqrt(max(m2, 0.0)))


def apparent_plate_number(first_moment: float, sigma: float) -> float:
    """Apparent plate number ``N = (μ₁/σ)²`` from peak moments."""
    if sigma <= 0 or not np.isfinite(sigma):
        return float("nan")
    return float((first_moment / sigma) ** 2)


@dataclass
class MethaneTransport:
    """Transport characterisation for one block, derived from methane peaks."""

    block: str
    n_markers: int
    t0_min: float                 # mean first moment (void time)
    t0_sd_min: float              # reproducibility across markers
    t0_range_min: float
    sigma_mean_min: float
    plate_number: float           # mean apparent N
    plate_number_sd: float
    t_inj_min: float              # effective inlet width (fitted)
    void_volume_m3: float
    fit_rmse: float               # methane forward-model fit quality

    def to_params(self) -> TransportParams:
        return TransportParams(t0_min=self.t0_min, t_inj_min=self.t_inj_min,
                               plate_number=self.plate_number)


def characterize_methane_transport(
    block: str,
    traces: list[tuple[np.ndarray, np.ndarray]],
    flow_col_m3_min: float,
    sample_mass_g: float,
    n_cells: int = 120,
) -> MethaneTransport:
    """Estimate block transport parameters from that block's methane peaks.

    Methane is unretained, so its peak is the system response with
    ``dq/dc = 0``.  Moments give ``t0`` and an apparent plate number directly;
    the effective inlet width ``t_inj`` is then chosen so the forward
    (no-adsorption) solve best reproduces the observed methane peaks.

    Parameters
    ----------
    block : str
        Block label (kept per block — never pooled across blocks).
    traces : list of (time_min, corrected_signal)
        The block's baseline-corrected methane chromatograms.
    flow_col_m3_min, sample_mass_g : float
        Measured conditions for this block.
    n_cells : int
        Spatial discretisation for the calibration solves.

    Returns
    -------
    MethaneTransport
    """
    from igc_sea.analysis.isotherm_models import NO_ADSORPTION

    m1s, sigmas, Ns = [], [], []
    for t, s in traces:
        _, m1, sd = peak_moments(t, s)
        if np.isfinite(m1) and np.isfinite(sd):
            m1s.append(m1)
            sigmas.append(sd)
            Ns.append(apparent_plate_number(m1, sd))

    if not m1s:
        raise ValueError(f"{block}: no usable methane peaks")

    t0 = float(np.mean(m1s))
    sigma_mean = float(np.mean(sigmas))
    N_mean = float(np.mean(Ns))

    geom = make_geometry(sample_mass_g, flow_col_m3_min, t0)

    # Choose the effective inlet width that best reproduces the methane peaks.
    # The total observed variance is shared between the rectangular inlet
    # (variance t_inj²/12) and column dispersion (variance t0²/N); scanning
    # t_inj and refitting N keeps that split explicit rather than arbitrary.
    best = None
    for frac in np.linspace(0.0, 0.9, 19):
        var_total = sigma_mean ** 2
        var_inlet = frac * var_total
        t_inj = float(np.sqrt(12.0 * var_inlet)) if var_inlet > 0 else 1e-6
        var_col = max(var_total - var_inlet, 1e-12)
        N_eff = (t0 ** 2) / var_col
        tp = TransportParams(t0_min=t0, t_inj_min=t_inj, plate_number=N_eff)

        rmse_acc = []
        for t, s in traces:
            area = float(np.trapezoid(np.maximum(s, 0.0), t))
            if area <= 0:
                continue
            # Unit "dose": compare shapes after area normalisation.
            res = solve_column(t, 1.0, tp, geom, NO_ADSORPTION,
                               np.array([]), n_cells=n_cells)
            pred = res.c_out
            pa = float(np.trapezoid(pred, t))
            if pa <= 0:
                continue
            obs_n = np.maximum(s, 0.0) / area
            pred_n = pred / pa
            scale = max(float(np.max(obs_n)), 1e-30)
            rmse_acc.append(float(np.sqrt(np.mean(((obs_n - pred_n) / scale) ** 2))))
        if not rmse_acc:
            continue
        rmse = float(np.mean(rmse_acc))
        if best is None or rmse < best[0]:
            best = (rmse, t_inj, N_eff)

    if best is None:
        best = (float("nan"), 1e-6, N_mean)

    rmse, t_inj_best, N_best = best

    return MethaneTransport(
        block=block,
        n_markers=len(m1s),
        t0_min=t0,
        t0_sd_min=(float(np.std(m1s, ddof=1)) if len(m1s) >= 2
                   else float("nan")),
        t0_range_min=float(np.max(m1s) - np.min(m1s)),
        sigma_mean_min=sigma_mean,
        plate_number=float(N_best),
        plate_number_sd=(float(np.std(Ns, ddof=1)) if len(Ns) >= 2
                         else float("nan")),
        t_inj_min=float(t_inj_best),
        void_volume_m3=geom.void_volume_m3,
        fit_rmse=float(rmse),
    )
