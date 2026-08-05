import importlib.util
from pathlib import Path, PurePosixPath


BOUNDARY_SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "check_public_release_boundary.py"
SPEC = importlib.util.spec_from_file_location("check_public_release_boundary", BOUNDARY_SCRIPT)
assert SPEC is not None and SPEC.loader is not None
BOUNDARY_MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(BOUNDARY_MODULE)

RESERVED_PREFIXES = BOUNDARY_MODULE.RESERVED_PREFIXES
SYNTHETIC_DATA_PREFIX = BOUNDARY_MODULE.SYNTHETIC_DATA_PREFIX
SYNTHETIC_BET_DATA_PREFIX = BOUNDARY_MODULE.SYNTHETIC_BET_DATA_PREFIX
SYNTHETIC_DISPERSIVE_DATA_PREFIX = (
    BOUNDARY_MODULE.SYNTHETIC_DISPERSIVE_DATA_PREFIX
)
forbidden_text_labels = BOUNDARY_MODULE.forbidden_text_labels
is_allowed_data_file = BOUNDARY_MODULE.is_allowed_data_file
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
    assert is_allowed_data_file(SYNTHETIC_BET_DATA_PREFIX / "traces.csv")
    assert is_allowed_data_file(SYNTHETIC_DISPERSIVE_DATA_PREFIX / "traces.csv")


def test_only_contract_fixture_data_files_are_allowed():
    assert is_allowed_data_file(SYNTHETIC_DATA_PREFIX / "traces.csv")
    assert is_allowed_data_file(SYNTHETIC_BET_DATA_PREFIX / "traces.csv")
    assert is_allowed_data_file(SYNTHETIC_DISPERSIVE_DATA_PREFIX / "traces.csv")
    assert not is_allowed_data_file(PurePosixPath("tests/fixtures/measurements.csv"))
    assert not is_allowed_data_file(PurePosixPath("study.xlsx"))
    assert is_allowed_data_file(PurePosixPath("docs/architecture.md"))


def test_absolute_user_paths_are_rejected():
    mac_path = b"/" + b"Users" + b"/person/private.csv"
    win_path = b"C:\\" + b"Users" + b"\\person\\private.csv"

    assert forbidden_text_labels(mac_path)
    assert forbidden_text_labels(win_path)
    assert not forbidden_text_labels(b"synthetic neutral chromatography bundle")
