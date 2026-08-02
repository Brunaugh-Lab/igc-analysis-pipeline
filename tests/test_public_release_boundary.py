import importlib.util
from pathlib import Path, PurePosixPath


BOUNDARY_SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "check_public_release_boundary.py"
SPEC = importlib.util.spec_from_file_location("check_public_release_boundary", BOUNDARY_SCRIPT)
assert SPEC is not None and SPEC.loader is not None
BOUNDARY_MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(BOUNDARY_MODULE)

RESERVED_PREFIXES = BOUNDARY_MODULE.RESERVED_PREFIXES
is_reserved = BOUNDARY_MODULE.is_reserved


def test_reserved_roots_and_descendants_are_rejected():
    for prefix in RESERVED_PREFIXES:
        assert is_reserved(prefix)
        assert is_reserved(prefix / "example.csv")


def test_public_source_and_synthetic_contract_fixture_are_allowed():
    assert not is_reserved(PurePosixPath("src/igc_analysis/analysis/full_peak.py"))
    assert not is_reserved(
        PurePosixPath(
            "src/igc_analysis/contracts/igc-neutral-data/0.2.0/"
            "examples/synthetic_peak_shape/traces.csv"
        )
    )
