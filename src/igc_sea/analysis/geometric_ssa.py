"""Geometric surface-area dosing surrogate from volume distributions."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from igc_sea.io.particle_size import VolumeDistribution, read_cumulative_q3


@dataclass(frozen=True)
class GeometricSSAResult:
    source_file: str
    formulation_id: str
    product: str
    measurement_time: str
    dispersing_system: str
    d32_um: float
    d32_bin_lower_um: float
    d32_bin_upper_um: float
    cumulative_start_percent: float
    cumulative_end_percent: float
    observed_volume_percent: float
    first_bin_volume_percent: float
    n_intervals: int
    density_g_cm3: float | None
    density_basis: str
    ssa_geo_m2_g: float | None
    ssa_geo_bin_low_m2_g: float | None
    ssa_geo_bin_high_m2_g: float | None
    ssa_dose_m2_g: float | None
    qc_flags: str

    def as_dict(self) -> dict:
        return asdict(self)


def _d32_from_volume(diameter_um: np.ndarray, volume: np.ndarray) -> float:
    total = float(np.sum(volume))
    denominator = float(np.sum(volume / diameter_um))
    return total / denominator if total > 0 and denominator > 0 else float("nan")


def calculate_d32(
    distribution: VolumeDistribution,
    *,
    density_g_cm3: float | None = None,
    density_basis: str = "unspecified",
    first_bin_warning_percent: float = 10.0,
    endpoint_tolerance_percent: float = 0.5,
) -> GeometricSSAResult:
    """Calculate D(3,2) from cumulative Q3 bin increments.

    The central estimate assigns each bin's volume to its geometric midpoint.
    Lower/upper values assign it to the bin edges and quantify bin-position
    sensitivity. ``ssa_dose_m2_g`` is the lower geometric SSA bound (upper-edge
    D32) and is conservative only with respect to PSD bin placement.
    """

    x = distribution.boundaries_um
    q = distribution.cumulative_q3_percent
    if np.any(~np.isfinite(x)) or np.any(~np.isfinite(q)):
        raise ValueError("Q3 boundaries and cumulative values must be finite")
    if np.any(x <= 0) or np.any(np.diff(x) <= 0):
        raise ValueError("Q3 boundaries must be positive and strictly increasing")
    if np.any(q < -endpoint_tolerance_percent) or np.any(q > 100 + endpoint_tolerance_percent):
        raise ValueError("Cumulative Q3 values fall outside 0-100%")

    increments_percent = np.diff(q)
    if np.any(increments_percent < -1e-8):
        raise ValueError("Cumulative Q3 must be monotonic nondecreasing")
    increments_percent = np.clip(increments_percent, 0.0, None)
    observed = float(np.sum(increments_percent))
    if observed <= 0:
        raise ValueError("Cumulative Q3 contains no positive volume increments")

    midpoint = np.sqrt(x[:-1] * x[1:])
    d32_mid = _d32_from_volume(midpoint, increments_percent)
    d32_lower = _d32_from_volume(x[:-1], increments_percent)
    d32_upper = _d32_from_volume(x[1:], increments_percent)

    first_positive = np.flatnonzero(increments_percent > 0)
    first_bin = float(increments_percent[first_positive[0]]) if len(first_positive) else 0.0
    flags: list[str] = []
    if abs(float(q[0])) > endpoint_tolerance_percent:
        flags.append("Q3_START_NOT_ZERO")
    if abs(float(q[-1]) - 100.0) > endpoint_tolerance_percent:
        flags.append("Q3_END_NOT_100")
    if first_bin > first_bin_warning_percent:
        flags.append("LOW_END_SENSITIVE")

    if density_g_cm3 is not None and density_g_cm3 <= 0:
        raise ValueError("density_g_cm3 must be positive")
    ssa_mid = 6.0 / (density_g_cm3 * d32_mid) if density_g_cm3 else None
    # Larger assigned diameters produce the lower geometric SSA estimate.
    ssa_low = 6.0 / (density_g_cm3 * d32_upper) if density_g_cm3 else None
    ssa_high = 6.0 / (density_g_cm3 * d32_lower) if density_g_cm3 else None

    meta = distribution.metadata
    return GeometricSSAResult(
        source_file=str(distribution.source_file),
        formulation_id=meta.get("formulation_id", meta.get("Product", distribution.source_file.stem)),
        product=meta.get("Product", ""),
        measurement_time=meta.get("Time", ""),
        dispersing_system=meta.get("Dispersing system", ""),
        d32_um=d32_mid,
        d32_bin_lower_um=d32_lower,
        d32_bin_upper_um=d32_upper,
        cumulative_start_percent=float(q[0]),
        cumulative_end_percent=float(q[-1]),
        observed_volume_percent=observed,
        first_bin_volume_percent=first_bin,
        n_intervals=int(np.count_nonzero(increments_percent > 0)),
        density_g_cm3=density_g_cm3,
        density_basis=density_basis if density_g_cm3 else "",
        ssa_geo_m2_g=ssa_mid,
        ssa_geo_bin_low_m2_g=ssa_low,
        ssa_geo_bin_high_m2_g=ssa_high,
        ssa_dose_m2_g=ssa_low,
        qc_flags=";".join(flags),
    )


def analyze_distribution_path(
    path: str | Path,
    *,
    density_g_cm3: float | None = None,
    density_basis: str = "unspecified",
    first_bin_warning_percent: float = 10.0,
) -> pd.DataFrame:
    """Analyze one distribution CSV or every CSV beneath a directory."""

    path = Path(path)
    files = [path] if path.is_file() else sorted(path.rglob("*.csv"))
    if not files:
        raise FileNotFoundError(f"No CSV files found at {path}")

    rows = []
    failures = []
    for file in files:
        try:
            distribution = read_cumulative_q3(file)
            rows.append(calculate_d32(
                distribution,
                density_g_cm3=density_g_cm3,
                density_basis=density_basis,
                first_bin_warning_percent=first_bin_warning_percent,
            ).as_dict())
        except (OSError, ValueError) as exc:
            failures.append(f"{file}: {exc}")
    if not rows:
        raise ValueError("No valid cumulative Q3 files were found:\n" + "\n".join(failures))
    result = pd.DataFrame(rows)
    result.attrs["failures"] = failures
    return result


def summarize_replicates(per_file: pd.DataFrame) -> pd.DataFrame:
    """Summarize file-level D32/SSA results by formulation."""

    metrics = [
        "d32_um", "d32_bin_lower_um", "d32_bin_upper_um",
        "first_bin_volume_percent", "ssa_geo_m2_g", "ssa_dose_m2_g",
    ]
    available = [metric for metric in metrics if per_file[metric].notna().any()]
    aggregations: dict[str, tuple[str, str]] = {
        "n_measurements": ("source_file", "count"),
        "n_flagged": ("qc_flags", lambda values: int((values != "").sum())),
    }
    for metric in available:
        aggregations[f"{metric}_mean"] = (metric, "mean")
        aggregations[f"{metric}_sd"] = (metric, "std")
    return (
        per_file.groupby("formulation_id", dropna=False)
        .agg(**aggregations)
        .reset_index()
    )
