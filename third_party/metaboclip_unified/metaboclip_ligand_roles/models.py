from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional


@dataclass
class FunctionalGroupMatch:
    group_id: str
    instance_id: str
    roles: Dict[str, int]
    atoms: List[int] = field(default_factory=list)
    role_metadata: Dict[str, Dict[str, Any]] = field(default_factory=dict)
    evidence: Dict[str, Any] = field(default_factory=dict)
    subtype: Optional[str] = None
    confidence: float = 1.0
    priority: int = 0

    def to_dict(self) -> Dict[str, Any]:
        out: Dict[str, Any] = {
            "group_id": self.group_id,
            "instance_id": self.instance_id,
            "roles": self.roles,
            "atoms": self.atoms,
            "role_metadata": self.role_metadata,
            "evidence": self.evidence,
            "confidence": self.confidence,
            "priority": self.priority,
        }
        if self.subtype is not None:
            out["subtype"] = self.subtype
        return out


@dataclass
class LigandAnnotation:
    ligand_id: str
    source_file: str
    index_base: int = 0
    framework_policy: Dict[str, Any] = field(default_factory=dict)
    functional_groups: List[FunctionalGroupMatch] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "ligand_id": self.ligand_id,
            "source_file": self.source_file,
            "index_base": self.index_base,
            "framework_policy": self.framework_policy,
            "functional_groups": [m.to_dict() for m in self.functional_groups],
        }
