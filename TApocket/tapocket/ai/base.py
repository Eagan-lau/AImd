from __future__ import annotations

from pathlib import Path
from typing import Any, Protocol

from tapocket.core.schema import CandidatePocket


class AIPocketPredictor(Protocol):
    name: str

    def predict(self, query_pdb: str | Path, output_dir: str | Path, run_id: str | None = None) -> list[CandidatePocket]:
        ...
