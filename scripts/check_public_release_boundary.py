"""Check tracked and unignored files for repository-content violations.

Ignored local data are allowed. Tracked, force-added, and new unignored files
are checked before a change is released.
"""

from __future__ import annotations

import subprocess
from pathlib import Path, PurePosixPath


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
RESERVED_PREFIXES = (
    PurePosixPath("governed_data"),
    PurePosixPath("local_bundles"),
    PurePosixPath("data/examples"),
)

SYNTHETIC_DATA_PREFIX = PurePosixPath(
    "src/igc_analysis/contracts/igc-neutral-data/0.2.0/examples/synthetic_peak_shape"
)
SYNTHETIC_BET_DATA_PREFIX = PurePosixPath(
    "src/igc_analysis/contracts/igc-neutral-data/0.2.0/examples/synthetic_bet_isotherm"
)
SYNTHETIC_DISPERSIVE_DATA_PREFIX = PurePosixPath(
    "src/igc_analysis/contracts/igc-neutral-data/0.2.0/examples/"
    "synthetic_dispersive_profile"
)
SYNTHETIC_DATA_PREFIXES = (
    SYNTHETIC_DATA_PREFIX,
    SYNTHETIC_BET_DATA_PREFIX,
    SYNTHETIC_DISPERSIVE_DATA_PREFIX,
)
DATA_EXTENSIONS = {
    ".csv",
    ".db",
    ".sqlite",
    ".tsv",
    ".xls",
    ".xlsx",
}

# Assemble path markers in pieces so this guard does not match its own source.
FORBIDDEN_TEXT_MARKERS = (
    ("macOS absolute user path", b"/" + b"Users" + b"/"),
    ("Linux absolute home path", b"/" + b"home" + b"/"),
    ("Windows absolute user path", b":\\" + b"Users" + b"\\"),
    ("institutional Dropbox path", b"Dropbox" + b"-UniversityofMichigan"),
    ("institutional Dropbox path", b"University of Michigan" + b" Dropbox"),
)


def candidate_paths() -> list[PurePosixPath]:
    result = subprocess.run(
        ["git", "ls-files", "-z", "--cached", "--others", "--exclude-standard"],
        cwd=REPOSITORY_ROOT,
        check=True,
        capture_output=True,
    )
    return [
        PurePosixPath(raw.decode("utf-8"))
        for raw in result.stdout.split(b"\0")
        if raw
    ]


def is_reserved(path: PurePosixPath) -> bool:
    return any(path == prefix or prefix in path.parents for prefix in RESERVED_PREFIXES)


def is_allowed_data_file(path: PurePosixPath) -> bool:
    if path.suffix.lower() not in DATA_EXTENSIONS:
        return True
    return any(path == prefix or prefix in path.parents
               for prefix in SYNTHETIC_DATA_PREFIXES)


def forbidden_text_labels(content: bytes) -> list[str]:
    return [label for label, marker in FORBIDDEN_TEXT_MARKERS if marker in content]


def main() -> int:
    paths = candidate_paths()
    reserved_violations = sorted(str(path) for path in paths if is_reserved(path))
    data_violations = sorted(str(path) for path in paths if not is_allowed_data_file(path))
    text_violations: list[tuple[str, list[str]]] = []
    for path in paths:
        content = (REPOSITORY_ROOT / path).read_bytes()
        labels = forbidden_text_labels(content)
        if labels:
            text_violations.append((str(path), labels))

    if reserved_violations:
        print("Repository content check: reserved local-data paths are tracked:")
        for path in reserved_violations:
            print(f"  - {path}")

    if data_violations:
        print("Repository content check: unexpected data files are tracked:")
        for path in data_violations:
            print(f"  - {path}")

    if text_violations:
        print("Repository content check: tracked files contain machine-local paths:")
        for path, labels in text_violations:
            print(f"  - {path}: {', '.join(labels)}")

    if reserved_violations or data_violations or text_violations:
        print("Move study data and machine-specific content outside the repository.")
        return 1

    print("Repository content check: OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
