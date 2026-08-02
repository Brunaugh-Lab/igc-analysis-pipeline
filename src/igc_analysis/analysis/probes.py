"""Source-neutral probe-property data structures.

Probe identities and physical properties enter public analyses through a
validated neutral bundle. This module intentionally contains no acquisition
format readers, filename rules, or source-schema knowledge.
"""

from __future__ import annotations

import math
from dataclasses import dataclass


SUPPORTED_ISOTHERM_PROBES = ("HEXANE", "HEPTANE", "OCTANE")


class ProbeError(ValueError):
    """Raised when declared probe properties are missing or nonphysical."""


@dataclass(frozen=True)
class ProbeProperties:
    """Physical properties declared for an isotherm probe.

    The five vapour-pressure coefficients use the explicit equation

    ``ln(P_sat[Pa]) = C1 + C2/T + C3*ln(T) + C4*T**C5``

    with temperature in kelvin. The property source must be retained in the
    neutral bundle that supplied these values.
    """

    name: str
    cross_section_m2: float
    antoine_c1: float
    antoine_c2: float
    antoine_c3: float
    antoine_c4: float
    antoine_c5: float
    carbon_number: float | None = None
    molecular_mass: float | None = None
    boiling_pt_C: float | None = None

    @property
    def antoine(self) -> dict[str, float]:
        return {
            "C1": self.antoine_c1,
            "C2": self.antoine_c2,
            "C3": self.antoine_c3,
            "C4": self.antoine_c4,
            "C5": self.antoine_c5,
        }

    def p_sat(self, temperature_K: float) -> float:
        """Evaluate saturation vapour pressure in pascals."""

        ln_psat = (
            self.antoine_c1
            + self.antoine_c2 / temperature_K
            + self.antoine_c3 * math.log(temperature_K)
            + self.antoine_c4 * temperature_K**self.antoine_c5
        )
        return math.exp(ln_psat)

    def validate(self) -> None:
        """Raise :class:`ProbeError` when a required value is unusable."""

        if not (
            self.cross_section_m2
            and self.cross_section_m2 > 0
            and math.isfinite(self.cross_section_m2)
        ):
            raise ProbeError(
                f"{self.name}: cross_section_m2={self.cross_section_m2!r} "
                "must be finite and positive"
            )
        for key, value in self.antoine.items():
            if value is None or not math.isfinite(value):
                raise ProbeError(
                    f"{self.name}: vapour-pressure coefficient {key}={value!r} "
                    "must be finite"
                )
        try:
            pressure = self.p_sat(303.15)
        except (ValueError, OverflowError) as exc:
            raise ProbeError(
                f"{self.name}: saturation-pressure evaluation failed: {exc}"
            ) from exc
        if not (math.isfinite(pressure) and pressure > 0):
            raise ProbeError(
                f"{self.name}: saturation pressure at 303.15 K is nonphysical"
            )


@dataclass
class ProbeSelection:
    """A declared probe choice and its source-neutral provenance."""

    probe: str
    source: str
