"""Fail when Git tracks files reserved for local or governed data.

This check intentionally operates on Git's index rather than the working tree:
ignored local data are allowed, while force-added files are rejected.
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


def tracked_paths() -> list[PurePosixPath]:
    result = subprocess.run(
        ["git", "ls-files", "-z"],
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


def main() -> int:
    violations = sorted(str(path) for path in tracked_paths() if is_reserved(path))
    if violations:
        print("Public-release boundary violation: Git tracks reserved local-data paths:")
        for path in violations:
            print(f"  - {path}")
        print("Move these files outside the public repository and remove them from Git history.")
        return 1

    print("Public-release tracked-path boundary: OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
