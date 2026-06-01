from __future__ import annotations

from pathlib import Path

from tapocket.core.schema import MCSARecord
from tapocket.utils.fs import resolve_path
from tapocket.utils.jsonl import read_jsonl
from tapocket.retrievers.foldseek import normalize_foldseek_target_id


class MCSADB:
    def __init__(self, records: list[MCSARecord], root: str | Path):
        self.root = Path(root).resolve()
        self.records = records
        self.by_id: dict[str, MCSARecord] = {}
        for record in records:
            aliases = {
                record.mcsa_id,
                normalize_foldseek_target_id(record.mcsa_id),
                normalize_foldseek_target_id(Path(record.reference_protein_path).stem),
            }
            for alias in aliases:
                self.by_id[alias] = record

    @classmethod
    def from_manifest(cls, manifest_path: str | Path, root: str | Path) -> "MCSADB":
        records = [MCSARecord(**record) for record in read_jsonl(manifest_path)]
        return cls(records, root)

    def get(self, mcsa_id: str) -> MCSARecord | None:
        norm = normalize_foldseek_target_id(mcsa_id)
        return self.by_id.get(mcsa_id) or self.by_id.get(norm)

    def resolve(self, path: str | Path) -> Path:
        return resolve_path(path, self.root)

    def reference_path(self, record: MCSARecord) -> Path:
        return self.resolve(record.reference_protein_path)

    def active_site_path(self, record: MCSARecord) -> Path:
        return self.resolve(record.active_site_path)
