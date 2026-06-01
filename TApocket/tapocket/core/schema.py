from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Any, Optional


@dataclass
class TemplateRecord:
    template_id: str
    group_dir: str
    protein_path: str
    pocket_paths: list[str]
    union_pocket_path: Optional[str] = None
    source_sample_ids_path: Optional[str] = None
    source_samples_csv_path: Optional[str] = None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class MCSARecord:
    mcsa_id: str
    group_dir: str
    reference_protein_path: str
    active_site_path: str

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class FoldseekHit:
    query: str
    target: str
    raw: dict[str, Any] = field(default_factory=dict)
    rank: int = 0
    normalized_target_id: Optional[str] = None
    score: float = 0.0

    def get_float(self, name: str, default: float = 0.0) -> float:
        value = self.raw.get(name, default)
        try:
            return float(value)
        except (TypeError, ValueError):
            return default

    @property
    def qtmscore(self) -> float:
        return self.get_float("qtmscore")

    @property
    def ttmscore(self) -> float:
        return self.get_float("ttmscore")

    @property
    def alntmscore(self) -> float:
        return self.get_float("alntmscore")

    @property
    def evalue(self) -> float:
        return self.get_float("evalue", 999.0)

    @property
    def bits(self) -> float:
        return self.get_float("bits")

    def compute_score(self) -> float:
        self.score = max(self.qtmscore, self.ttmscore, self.alntmscore, 0.0)
        return self.score

    def to_dict(self) -> dict[str, Any]:
        return {
            "query": self.query,
            "target": self.target,
            "normalized_target_id": self.normalized_target_id,
            "rank": self.rank,
            "score": self.score,
            "raw": self.raw,
        }


@dataclass(frozen=True)
class ResidueKey:
    chain: str
    resi: str
    icode: str = ""
    resn: str = ""

    def label(self) -> str:
        icode = self.icode if self.icode else ""
        chain = self.chain if self.chain else "_"
        return f"{self.resn}:{chain}:{self.resi}{icode}"

    def to_dict(self) -> dict[str, str]:
        return {"chain": self.chain, "resi": self.resi, "icode": self.icode, "resn": self.resn}


@dataclass
class QueryPocketResidue:
    query_id: str
    template_id: str
    template_pocket_id: str
    chain: str
    resi: str
    icode: str
    resn: str
    min_distance_to_mapped_pocket: float

    def key(self) -> tuple[str, str, str, str]:
        return (self.chain, self.resi, self.icode, self.resn)

    def label(self) -> str:
        chain = self.chain or "_"
        icode = self.icode or ""
        return f"{chain}:{self.resi}{icode}:{self.resn}"

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class CandidatePocket:
    query_id: str
    template_id: str
    template_pocket_id: str
    source: str
    mapped_pocket_path: str = ""
    template_protein_path: str = ""
    template_pocket_path: str = ""
    foldseek_hit: dict[str, Any] = field(default_factory=dict)
    query_pocket_residues: list[QueryPocketResidue] = field(default_factory=list)
    mapping_quality: dict[str, Any] = field(default_factory=dict)
    mcsa_support: dict[str, Any] = field(default_factory=dict)
    ai_support: dict[str, Any] = field(default_factory=dict)
    box: dict[str, Any] = field(default_factory=dict)
    final_score: float = 0.0
    kept: bool = True
    drop_reason: Optional[str] = None
    custom_pocket_id: Optional[str] = None
    cluster_id: Optional[str] = None
    supporting_pockets: list[dict[str, Any]] = field(default_factory=list)
    extra: dict[str, Any] = field(default_factory=dict)

    @property
    def pocket_id(self) -> str:
        if self.custom_pocket_id:
            return self.custom_pocket_id
        return f"{self.template_id}:{self.template_pocket_id}"

    def to_dict(self) -> dict[str, Any]:
        return {
            "query_id": self.query_id,
            "pocket_id": self.pocket_id,
            "cluster_id": self.cluster_id,
            "source": self.source,
            "template_id": self.template_id,
            "template_pocket_id": self.template_pocket_id,
            "mapped_pocket_path": self.mapped_pocket_path,
            "template_protein_path": self.template_protein_path,
            "template_pocket_path": self.template_pocket_path,
            "foldseek_hit": self.foldseek_hit,
            "query_pocket_residues": [r.to_dict() for r in self.query_pocket_residues],
            "mapping_quality": self.mapping_quality,
            "mcsa_support": self.mcsa_support,
            "ai_support": self.ai_support,
            "box": self.box,
            "final_score": self.final_score,
            "kept": self.kept,
            "drop_reason": self.drop_reason,
            "supporting_pockets": self.supporting_pockets,
            "extra": self.extra,
        }


@dataclass
class PipelineRunSummary:
    run_id: str
    query_id: str
    query_path: str
    selected_template_count: int
    candidate_pocket_count: int
    final_pocket_count: int
    mcsa_enabled: bool
    mcsa_action: str
    ai_enabled: bool
    ai_used: bool
    output_mode: str
    run_dir: str

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)
