from __future__ import annotations

from pathlib import Path

from tapocket.core.schema import TemplateRecord
from tapocket.utils.fs import resolve_path
from tapocket.utils.jsonl import read_jsonl
from tapocket.retrievers.foldseek import normalize_foldseek_target_id


class TemplateDB:
    def __init__(self, records: list[TemplateRecord], root: str | Path):
        self.root = Path(root).resolve()
        self.records = records
        self.by_id: dict[str, TemplateRecord] = {}
        for record in records:
            aliases = {
                record.template_id,
                normalize_foldseek_target_id(record.template_id),
                normalize_foldseek_target_id(Path(record.protein_path).stem),
            }
            for alias in aliases:
                self.by_id[alias] = record

    @classmethod
    def from_manifest(cls, manifest_path: str | Path, root: str | Path) -> "TemplateDB":
        records = [TemplateRecord(**record) for record in read_jsonl(manifest_path)]
        return cls(records, root)

    def get(self, template_id: str) -> TemplateRecord | None:
        norm = normalize_foldseek_target_id(template_id)
        return self.by_id.get(template_id) or self.by_id.get(norm)

    def resolve(self, path: str | Path) -> Path:
        return resolve_path(path, self.root)

    def protein_path(self, record: TemplateRecord) -> Path:
        return self.resolve(record.protein_path)

    def pocket_paths(self, record: TemplateRecord) -> list[Path]:
        return [self.resolve(p) for p in record.pocket_paths]
