"""Reader for cumulative particle-volume distributions in CSV form."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
import re

import numpy as np


_SUBSCRIPTS = str.maketrans("₀₁₂₃₄₅₆₇₈₉", "0123456789")


@dataclass(frozen=True)
class VolumeDistribution:
    """One cumulative volume distribution and its optional metadata."""

    source_file: Path
    metadata: dict[str, str]
    boundaries_um: np.ndarray
    cumulative_q3_percent: np.ndarray


def _normalise_header(value: str) -> str:
    return value.translate(_SUBSCRIPTS).replace("μ", "µ").strip()


def _find_distribution_header(rows: list[list[str]]) -> int:
    for index, row in enumerate(rows[:40]):
        headers = [_normalise_header(cell) for cell in row]
        has_boundary = any(re.match(r"^xo\s*/", h, re.IGNORECASE) for h in headers)
        has_cumulative_q3 = any(
            h.startswith("Q") and "3" in h and "%" in h for h in headers
        )
        if has_boundary and has_cumulative_q3:
            return index
    raise ValueError("Could not locate an xo/Q3 distribution header")


def read_cumulative_q3(path: str | Path) -> VolumeDistribution:
    """Read cumulative volume Q3 and metadata from a CSV file.

    The parser selects the uppercase cumulative ``Q3 / %`` column rather than
    lowercase differential-density ``q3`` columns.
    """

    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(path)

    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.reader(handle))
    if len(rows) < 4:
        raise ValueError(f"distribution file is too short: {path}")

    metadata: dict[str, str] = {}
    for key, value in zip(rows[0], rows[1]):
        key = key.lstrip("\ufeff").strip()
        if key:
            metadata[key] = value.strip()

    header_index = _find_distribution_header(rows)
    headers = [_normalise_header(cell) for cell in rows[header_index]]
    xo_index = next(
        i for i, header in enumerate(headers)
        if re.match(r"^xo\s*/", header, re.IGNORECASE)
    )
    q3_index = next(
        i for i, header in enumerate(headers)
        if header.startswith("Q") and "3" in header and "%" in header
    )

    boundaries: list[float] = []
    cumulative: list[float] = []
    for row in rows[header_index + 1:]:
        if max(xo_index, q3_index) >= len(row):
            continue
        try:
            boundary = float(row[xo_index])
            q3 = float(row[q3_index])
        except (TypeError, ValueError):
            continue
        if np.isfinite(boundary) and np.isfinite(q3):
            boundaries.append(boundary)
            cumulative.append(q3)

    if len(boundaries) < 2:
        raise ValueError(f"fewer than two valid Q3 boundaries in {path}")

    return VolumeDistribution(
        source_file=path,
        metadata=metadata,
        boundaries_um=np.asarray(boundaries, dtype=float),
        cumulative_q3_percent=np.asarray(cumulative, dtype=float),
    )
