"""Read validated IGC experiment bundles using the neutral-data contract."""

from __future__ import annotations

import importlib.util
import json
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType

import pandas as pd


CONTRACT_NAME = "igc-neutral-data"
SUPPORTED_CONTRACT_VERSION = "0.2.0"
CONTRACT_DIR = (
    Path(__file__).resolve().parent.parent
    / "contracts"
    / CONTRACT_NAME
    / SUPPORTED_CONTRACT_VERSION
)


class NeutralContractError(ValueError):
    """Raised when a neutral bundle does not satisfy the supported contract."""


@dataclass(frozen=True)
class NeutralBundle:
    """A structurally validated neutral bundle and its typed tables."""

    path: Path
    manifest: dict
    tables: dict[str, pd.DataFrame]

    @property
    def dataset_id(self) -> str:
        return str(self.manifest["dataset_id"])

    @property
    def contract_version(self) -> str:
        return str(self.manifest["contract_version"])

    def table(self, filename: str) -> pd.DataFrame:
        try:
            return self.tables[filename]
        except KeyError as exc:
            raise KeyError(f"neutral bundle has no declared table {filename!r}") from exc


def _load_validator() -> ModuleType:
    path = CONTRACT_DIR / "validate_bundle.py"
    spec = importlib.util.spec_from_file_location("igc_neutral_validator_0_2_0", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load bundled neutral validator from {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_neutral_bundle(bundle_dir: str | Path) -> NeutralBundle:
    """Validate and load one ``igc-neutral-data/0.2.0`` bundle.

    Structural validation—including manifest hashes, exact headers, units,
    controlled values, ordering, and cross-table keys—runs before any table is
    exposed to analysis code.
    """

    path = Path(bundle_dir).expanduser().resolve()
    validator = _load_validator()
    try:
        validator.validate_bundle(path)
    except (validator.ValidationError, OSError, ValueError) as exc:
        raise NeutralContractError(str(exc)) from exc

    try:
        manifest = json.loads((path / "manifest.json").read_text(encoding="utf-8-sig"))
    except (OSError, json.JSONDecodeError) as exc:
        raise NeutralContractError(f"cannot read neutral manifest: {exc}") from exc

    if manifest.get("contract_name") != CONTRACT_NAME:
        raise NeutralContractError(f"unsupported contract {manifest.get('contract_name')!r}")
    if manifest.get("contract_version") != SUPPORTED_CONTRACT_VERSION:
        raise NeutralContractError(
            f"unsupported {CONTRACT_NAME} version {manifest.get('contract_version')!r}; "
            f"supported version is {SUPPORTED_CONTRACT_VERSION}"
        )

    tables = {
        filename: pd.read_csv(path / filename, keep_default_na=False)
        for filename in manifest["files"]
    }
    return NeutralBundle(path=path, manifest=manifest, tables=tables)


def bundled_contract_path() -> Path:
    """Return the installed path of the supported source-neutral contract."""

    return CONTRACT_DIR
