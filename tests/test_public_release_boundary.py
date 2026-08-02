from pathlib import PurePosixPath

from scripts.check_public_release_boundary import RESERVED_PREFIXES, is_reserved


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
