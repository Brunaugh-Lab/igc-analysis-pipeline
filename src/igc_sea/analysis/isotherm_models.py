"""Adsorption isotherm models for finite-concentration inverse chromatography.

Each model supplies the two quantities the column-transport model needs —
the adsorbed amount ``q(c)`` and its slope ``dq/dc`` — plus parameter bounds,
units and numerical safeguards near ``c = 0``.

Units convention (used consistently across the full-peak pipeline):

- ``c`` : gas-phase concentration, mol/m³
- ``q`` : adsorbed amount per gram of sample, mol/g
- ``dq/dc`` : m³/g

**Finite capacity and SSA.** Only a model with a structural saturation
capacity can support a monolayer/SSA calculation.  ``has_finite_capacity``
records that, and ``capacity_param`` names the parameter.  A cooperative or
power-law isotherm has no monolayer, and the pipeline must not manufacture
one — see :func:`igc_sea.analysis.full_peak.compute_ssa_if_identifiable`.

References
----------
Guiochon, Felinger, Shirazi & Katti, *Fundamentals of Preparative and
Nonlinear Chromatography*, 2nd ed. (2006), ch. 3.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

# Concentration floor (mol/m³) used before evaluating power laws, so that a
# Freundlich exponent n < 1 cannot produce an infinite slope at c = 0.
C_FLOOR = 1e-12


@dataclass
class IsothermModel:
    """A named adsorption isotherm.

    Attributes
    ----------
    name : str
        Short identifier, e.g. ``"langmuir"``.
    param_names : list[str]
        Ordered parameter names.
    param_units : list[str]
        Units for each parameter, for reporting.
    bounds : tuple[list[float], list[float]]
        ``(lower, upper)`` optimiser bounds, same order as ``param_names``.
    initial : list[float]
        Default starting values.
    has_finite_capacity : bool
        True only if the model has a structural saturation capacity, i.e. a
        monolayer is defined.  Gates any SSA calculation.
    capacity_param : str | None
        Name of the saturation-capacity parameter, if any.
    description : str
        One-line physical description.
    """

    name: str
    param_names: list[str]
    param_units: list[str]
    bounds: tuple[list[float], list[float]]
    initial: list[float]
    has_finite_capacity: bool
    capacity_param: str | None
    description: str
    _q: object = field(repr=False, default=None)
    _dqdc: object = field(repr=False, default=None)

    @property
    def n_params(self) -> int:
        return len(self.param_names)

    def q(self, c: np.ndarray, params: np.ndarray) -> np.ndarray:
        """Adsorbed amount q(c) in mol/g (nonnegative)."""
        c = np.maximum(np.asarray(c, dtype=float), 0.0)
        out = self._q(c, np.asarray(params, dtype=float))
        return np.maximum(np.nan_to_num(out, nan=0.0, posinf=0.0), 0.0)

    def dqdc(self, c: np.ndarray, params: np.ndarray) -> np.ndarray:
        """Isotherm slope dq/dc in m³/g (finite, nonnegative)."""
        c = np.maximum(np.asarray(c, dtype=float), 0.0)
        out = self._dqdc(c, np.asarray(params, dtype=float))
        return np.maximum(np.nan_to_num(out, nan=0.0, posinf=0.0), 0.0)

    def henry_constant(self, params: np.ndarray,
                       c_ref: float = C_FLOOR) -> float:
        """Low-concentration affinity dq/dc as c → 0 (m³/g).

        Evaluated at ``c_ref`` rather than exactly 0 so power-law models with
        a singular or vanishing origin slope still return a finite reference
        value.  For Freundlich with n > 1 this tends to 0, which is itself the
        diagnostic that there is no Henry region.
        """
        return float(self.dqdc(np.array([c_ref]), params)[0])


# ---------------------------------------------------------------------------
# Model definitions
# ---------------------------------------------------------------------------

def _none_q(c, p):
    return np.zeros_like(c)


def _none_dqdc(c, p):
    return np.zeros_like(c)


NO_ADSORPTION = IsothermModel(
    name="none",
    param_names=[],
    param_units=[],
    bounds=([], []),
    initial=[],
    has_finite_capacity=False,
    capacity_param=None,
    description="No adsorption — transport-only control model (q ≡ 0).",
    _q=_none_q,
    _dqdc=_none_dqdc,
)


def _henry_q(c, p):
    return p[0] * c


def _henry_dqdc(c, p):
    return np.full_like(c, p[0])


HENRY = IsothermModel(
    name="henry",
    param_names=["K_H"],
    param_units=["m^3/g"],
    bounds=([0.0], [1.0]),
    initial=[1e-6],
    has_finite_capacity=False,
    capacity_param=None,
    description="Linear (Henry) adsorption: q = K_H·c; dq/dc constant, so "
                "retention is dose-independent.",
    _q=_henry_q,
    _dqdc=_henry_dqdc,
)


def _langmuir_q(c, p):
    qs, K = p[0], p[1]
    return qs * K * c / (1.0 + K * c)


def _langmuir_dqdc(c, p):
    qs, K = p[0], p[1]
    denom = 1.0 + K * c
    return qs * K / (denom * denom)


LANGMUIR = IsothermModel(
    name="langmuir",
    param_names=["q_s", "K_L"],
    param_units=["mol/g", "m^3/mol"],
    bounds=([1e-9, 1e-6], [1e-1, 1e6]),
    initial=[1e-5, 1e-2],
    has_finite_capacity=True,
    capacity_param="q_s",
    description="Langmuir (favourable, concave): finite monolayer q_s; dq/dc "
                "DECREASES with c, so peaks elute earlier at higher dose.",
    _q=_langmuir_q,
    _dqdc=_langmuir_dqdc,
)


def _freundlich_q(c, p):
    KF, n = p[0], p[1]
    cc = np.maximum(c, C_FLOOR)
    return KF * cc ** n


def _freundlich_dqdc(c, p):
    KF, n = p[0], p[1]
    cc = np.maximum(c, C_FLOOR)
    return n * KF * cc ** (n - 1.0)


FREUNDLICH = IsothermModel(
    name="freundlich",
    param_names=["K_F", "n"],
    param_units=["mol/g / (mol/m^3)^n", "-"],
    # n spans both branches so the fit can choose; n>1 is the cooperative one.
    bounds=([1e-12, 0.2], [1e2, 3.0]),
    initial=[1e-6, 1.0],
    has_finite_capacity=False,
    capacity_param=None,
    description="Freundlich power law q = K_F·c^n. n<1 favourable/concave; "
                "n=1 Henry; n>1 CONVEX / concentration-strengthened "
                "(cooperative) — dq/dc rises with c, peaks elute later at "
                "higher dose. No monolayer capacity exists for any n.",
    _q=_freundlich_q,
    _dqdc=_freundlich_dqdc,
)


MODELS: dict[str, IsothermModel] = {
    m.name: m for m in (NO_ADSORPTION, HENRY, LANGMUIR, FREUNDLICH)
}


def get_model(name: str) -> IsothermModel:
    """Look up an isotherm model by name.

    Raises
    ------
    KeyError
        If the model is not registered.
    """
    key = name.strip().lower()
    if key not in MODELS:
        raise KeyError(
            f"Unknown isotherm model {name!r}. Available: {sorted(MODELS)}")
    return MODELS[key]


def is_cooperative(model: IsothermModel, params: np.ndarray) -> bool:
    """True if the fitted isotherm is concentration-strengthened (convex).

    Defined as ``dq/dc`` increasing with ``c``: the Freundlich branch ``n > 1``.
    Langmuir and Henry are never cooperative by construction.
    """
    if model.name == "freundlich":
        return bool(params[1] > 1.0)
    return False
