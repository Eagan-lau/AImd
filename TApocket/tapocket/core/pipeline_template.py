from __future__ import annotations

import json
import shutil
import time
from pathlib import Path
from typing import Any

from tapocket.ai.script_adapter import ScriptAIPredictor
from tapocket.ai.deeppocket_db_adapter import DeepPocketDBPredictor
from tapocket.core.config import TapocketConfig
from tapocket.core.exporter import (
    export_boxes_tsv,
    export_candidates_json,
    export_final_pdb,
    export_hits_json,
    export_query_pocket_residues,
    export_run_summary,
    export_summary_tsv,
    write_json,
)
from tapocket.core.provenance import write_provenance
from tapocket.core.schema import CandidatePocket, PipelineRunSummary
from tapocket.databases.mcsa_db import MCSADB
from tapocket.databases.template_db import TemplateDB
from tapocket.decision.ai_trigger import should_run_ai
from tapocket.filters.mcsa_filter import MCSAFilter
from tapocket.filters.pocket_dedup import deduplicate_pockets
from tapocket.mapping.pocket_mapper import PocketMapper
from tapocket.retrievers.foldseek import FoldseekRunner, select_hits, write_hits_tsv
from tapocket.utils.fs import ensure_dir


def _query_id_from_path(query_path: str | Path) -> str:
    return Path(query_path).stem


def _write_log(log_path: Path, message: str) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8") as handle:
        handle.write(message.rstrip() + "\n")
    print(message)


def _copy_if_exists(src: Path, dst: Path) -> None:
    if src.exists():
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)


class TApocketPipeline:
    """Template-first TApocket pipeline with optional AI model fallback.

    Main branch:
      Foldseek -> template pocket mapping -> CA-to-CA residue extraction -> QC -> pocket dedup -> M-CSA box filter.

    AI branch:
      Optional script adapter. Default mode is fallback_only when enabled.
    """

    def __init__(self, config: TapocketConfig):
        self.config = config
        self.root = config.root

    def _run_ai_predictors(
        self,
        query_pdb: Path,
        ai_dir: Path,
        query_id: str,
        run_root: Path,
        run_id: str,
        log_path: Path,
    ) -> list[CandidatePocket]:
        predictors = self.config.get("ai_models", "predictors", default=[]) or []
        ai_candidates: list[CandidatePocket] = []
        for predictor_cfg in predictors:
            if not predictor_cfg.get("enabled", True):
                continue
            ptype = str(predictor_cfg.get("type", "script"))
            name = str(predictor_cfg.get("name", "ai_model"))
            model_dir = ensure_dir(ai_dir / name)
            if ptype == "script":
                predictor = ScriptAIPredictor(self.config, predictor_cfg, query_id=query_id, run_root=run_root)
            elif ptype == "deeppocket_db":
                predictor = DeepPocketDBPredictor(self.config, predictor_cfg, query_id=query_id, run_root=run_root)
            else:
                raise RuntimeError(f"Unsupported AI predictor type={ptype!r}. Use type: deeppocket_db or script.")
            _write_log(log_path, f"[AI] running predictor={name} type={ptype}")
            predicted = predictor.predict(query_pdb, model_dir, run_id=run_id)
            _write_log(log_path, f"[AI] predictor={name}, predicted_pockets={len(predicted)}")
            ai_candidates.extend(predicted)
        return ai_candidates

    def run(self, query: str | Path, run_id: str | None = None) -> PipelineRunSummary:
        query = Path(query).resolve()
        if not query.exists():
            raise FileNotFoundError(query)
        query_id = _query_id_from_path(query)
        if run_id is None:
            run_id = f"{query_id}_{time.strftime('%Y%m%d_%H%M%S')}"

        output_mode = str(self.config.get("output", "mode", default="compact"))
        if output_mode not in {"compact", "detailed"}:
            raise ValueError("output.mode must be compact or detailed")

        run_root = self.config.path("paths", "run_root") / run_id
        if output_mode == "compact":
            stage_root = ensure_dir(run_root / "_intermediate")
            input_dir = ensure_dir(stage_root / "input")
            foldseek_dir = ensure_dir(stage_root / "foldseek")
            mapping_dir = ensure_dir(stage_root / "template_mapping")
            mcsa_dir = ensure_dir(stage_root / "mcsa")
            ai_dir = ensure_dir(stage_root / "ai_model")
            final_dir = ensure_dir(run_root)
            log_path = run_root / "run.log"
        else:
            stage_root = ensure_dir(run_root)
            input_dir = ensure_dir(run_root / "input")
            foldseek_dir = ensure_dir(run_root / "foldseek")
            mapping_dir = ensure_dir(run_root / "template_mapping")
            mcsa_dir = ensure_dir(run_root / "mcsa")
            ai_dir = ensure_dir(run_root / "ai_model")
            final_dir = ensure_dir(run_root / "final")
            log_path = run_root / "logs" / "run.log"

        run_query = input_dir / "query.pdb"
        shutil.copy2(query, run_query)
        _write_log(log_path, f"[TApocket] run_id={run_id}")
        _write_log(log_path, f"[TApocket] query={query}")
        _write_log(log_path, f"[TApocket] output_mode={output_mode}")

        template_manifest = self.config.path("paths", "template_manifest")
        if not template_manifest.exists():
            raise FileNotFoundError(
                f"Template manifest not found: {template_manifest}\n"
                "Run: tapocket build-manifest --config configs/tapocket_template_v1.yaml"
            )
        template_db = TemplateDB.from_manifest(template_manifest, self.root)

        foldseek_runner = FoldseekRunner(
            binary=self.config.get("foldseek", "binary", default="foldseek"),
            tmp_dir=self.config.path("foldseek", "tmp_dir"),
            format_fields=self.config.get("foldseek", "format_output"),
            extra_args=self.config.get("foldseek", "extra_args", default=[]),
        )

        raw_hits_tsv = foldseek_dir / "template_hits.raw.tsv"
        _write_log(log_path, "[Step 1] Foldseek search against template DB")
        foldseek_runner.search(
            query_pdb=run_query,
            db_path=self.config.path("paths", "foldseek_template_db"),
            output_tsv=raw_hits_tsv,
            max_seqs=int(self.config.get("foldseek", "max_seqs", default=50)),
            log_json=foldseek_dir / "template_search_command.json",
        )
        raw_hits = foldseek_runner.parse_hits(raw_hits_tsv)
        selected_hits = select_hits(
            raw_hits,
            min_qtmscore=float(self.config.get("template_selection", "min_qtmscore", default=0.45)),
            min_ttmscore=float(self.config.get("template_selection", "min_ttmscore", default=0.45)),
            min_alntmscore=float(self.config.get("template_selection", "min_alntmscore", default=0.45)),
            keep_top_k=int(self.config.get("template_selection", "keep_top_k", default=10)),
            deduplicate_target=bool(self.config.get("template_deduplication", "keep_best_per_template", default=True)),
        )
        write_hits_tsv(selected_hits, foldseek_dir / "selected_template_hits.tsv")
        export_hits_json(foldseek_dir / "selected_template_hits.json", selected_hits)
        _write_log(log_path, f"[Step 1] selected non-redundant template hits={len(selected_hits)}")

        _write_log(log_path, "[Step 2] Map template pockets to query coordinate system")
        mapper = PocketMapper(self.config, template_db, query_id=query_id, run_root=stage_root)
        candidates: list[CandidatePocket] = []
        missing_templates: list[dict[str, Any]] = []
        for hit in selected_hits:
            record = template_db.get(hit.normalized_target_id or hit.target)
            if not record:
                missing_templates.append(hit.to_dict())
                continue
            try:
                candidates.extend(mapper.map_template_hit(run_query, hit, record, mapping_dir))
            except Exception as exc:
                _write_log(log_path, f"[Step 2][ERROR] template={record.template_id}: {exc}")
                missing_templates.append({"hit": hit.to_dict(), "error": str(exc)})

        export_candidates_json(mapping_dir / "mapped_template_pockets.all.json", candidates)
        candidates_before_qc = len(candidates)
        if str(self.config.get("pocket_mapping", "quality_control", "low_quality_action", default="drop")) == "drop":
            candidates = [c for c in candidates if c.kept]
        export_candidates_json(mapping_dir / "mapped_template_pockets.qc.json", candidates)
        export_query_pocket_residues(mapping_dir / "query_pocket_residues.tsv", candidates)
        (mapping_dir / "mapping_errors.json").write_text(json.dumps(missing_templates, indent=2, ensure_ascii=False), encoding="utf-8")
        _write_log(log_path, f"[Step 2] candidate pockets={len(candidates)} after QC, raw={candidates_before_qc}")

        dedup_report: dict[str, Any] = {"enabled": False, "cluster_count": len(candidates)}
        if bool(self.config.get("pocket_deduplication", "enabled", default=True)):
            _write_log(log_path, "[Step 3] Pocket-level deduplication")
            candidates, dedup_report = deduplicate_pockets(
                candidates,
                overlap_threshold=float(self.config.get("pocket_deduplication", "overlap_threshold", default=0.70)),
                representative=str(self.config.get("pocket_deduplication", "representative", default="best_score")),
            )
            write_json(mapping_dir / "pocket_deduplication_report.json", dedup_report)
            export_candidates_json(mapping_dir / "mapped_template_pockets.dedup.json", candidates)
            _write_log(log_path, f"[Step 3] pocket clusters={len(candidates)}")

        final_candidates = candidates
        mcsa_report: dict[str, Any] = {"mcsa_enabled": False, "filter_action": "skip_mcsa"}

        if bool(self.config.get("mcsa", "enabled", default=True)):
            _write_log(log_path, "[Step 4] M-CSA optional catalytic box filtering")
            mcsa_manifest = self.config.path("paths", "mcsa_manifest")
            if not mcsa_manifest.exists():
                raise FileNotFoundError(
                    f"M-CSA manifest not found: {mcsa_manifest}\n"
                    "Run: tapocket build-manifest --config configs/tapocket_template_v1.yaml"
                )
            mcsa_db = MCSADB.from_manifest(mcsa_manifest, self.root)
            mcsa_raw_tsv = mcsa_dir / "mcsa_hits.raw.tsv"
            foldseek_runner.search(
                query_pdb=run_query,
                db_path=self.config.path("paths", "foldseek_mcsa_db"),
                output_tsv=mcsa_raw_tsv,
                max_seqs=int(self.config.get("mcsa", "max_seqs", default=20)),
                log_json=mcsa_dir / "mcsa_search_command.json",
            )
            mcsa_raw_hits = foldseek_runner.parse_hits(mcsa_raw_tsv)
            mcsa_hits = select_hits(
                mcsa_raw_hits,
                min_qtmscore=float(self.config.get("mcsa", "min_qtmscore", default=0.40)),
                min_ttmscore=float(self.config.get("mcsa", "min_ttmscore", default=0.40)),
                min_alntmscore=float(self.config.get("mcsa", "min_alntmscore", default=0.40)),
                keep_top_k=int(self.config.get("mcsa", "keep_top_k", default=10)),
                deduplicate_target=True,
            )
            write_hits_tsv(mcsa_hits, mcsa_dir / "selected_mcsa_hits.tsv")
            export_hits_json(mcsa_dir / "selected_mcsa_hits.json", mcsa_hits)
            mcsa_filter = MCSAFilter(self.config, mcsa_db, query_id=query_id, run_root=stage_root)
            final_candidates, mcsa_report = mcsa_filter.filter_candidates(run_query, candidates, mcsa_hits, mcsa_dir)
            _write_log(log_path, f"[Step 4] mcsa_hits={len(mcsa_hits)}, action={mcsa_report.get('filter_action')}, final={len(final_candidates)}")

        ai_used = False
        ai_report: dict[str, Any] = {"enabled": bool(self.config.get("ai_models", "enabled", default=False)), "used": False}
        run_ai, ai_reason = should_run_ai(self.config, final_candidates, selected_template_count=len(selected_hits))
        ai_report["trigger_reason"] = ai_reason
        if run_ai:
            _write_log(log_path, f"[Step 5] AI model branch triggered: {ai_reason}")
            ai_candidates = self._run_ai_predictors(run_query, ai_dir, query_id, stage_root, run_id, log_path)
            ai_used = True
            ai_report["used"] = True
            ai_report["predicted_pocket_count"] = len(ai_candidates)
            ai_mode = str(self.config.get("ai_models", "mode", default="fallback_only"))
            if ai_mode == "fallback_only":
                if not final_candidates:
                    final_candidates = ai_candidates
            elif ai_mode == "merge":
                final_candidates = final_candidates + ai_candidates
            elif ai_mode in {"always", "rerank"}:
                # Keep template results as primary; store AI outputs for detailed inspection.
                # True reranking can be added after benchmark calibration.
                pass
            export_candidates_json(ai_dir / "ai_candidates.json", ai_candidates)
        else:
            _write_log(log_path, f"[Step 5] AI model branch skipped: {ai_reason}")

        _write_log(log_path, "[Step 6] Export final results")
        export_final_pdb(final_dir / "final_pockets.pdb", final_candidates, stage_root)

        export_candidates = final_candidates
        if output_mode == "compact":
            # Compact mode removes intermediate mapped-pocket files by default. Keep the
            # machine-readable output clean by pointing each pocket to the combined final PDB.
            import copy
            export_candidates = copy.deepcopy(final_candidates)
            for idx, cand in enumerate(export_candidates, start=1):
                cand.mapped_pocket_path = "final_pockets.pdb"
                cand.extra["final_pdb_model_index"] = idx

        export_candidates_json(final_dir / "final_pockets.json", export_candidates)
        export_query_pocket_residues(final_dir / "final_pocket_residues.tsv", export_candidates)
        export_boxes_tsv(final_dir / "final_boxes.tsv", export_candidates)
        export_summary_tsv(final_dir / "summary.tsv", export_candidates)

        summary = PipelineRunSummary(
            run_id=run_id,
            query_id=query_id,
            query_path=str(query),
            selected_template_count=len(selected_hits),
            candidate_pocket_count=len(candidates),
            final_pocket_count=len(final_candidates),
            mcsa_enabled=bool(self.config.get("mcsa", "enabled", default=True)),
            mcsa_action=str(mcsa_report.get("filter_action", "skip_mcsa")),
            ai_enabled=bool(self.config.get("ai_models", "enabled", default=False)),
            ai_used=ai_used,
            output_mode=output_mode,
            run_dir=str(run_root),
        )
        export_run_summary(run_root / "run_summary.json", summary)
        write_provenance(
            run_root / "provenance.json",
            config_path=self.config.config_path,
            project_root=self.root,
            run_id=run_id,
            query_path=query,
            extra={"mcsa_report": mcsa_report, "dedup_report": dedup_report, "ai_report": ai_report},
        )
        _write_log(log_path, "[TApocket] Done")

        if output_mode == "compact" and not bool(self.config.get("output", "compact", "keep_intermediate", default=False)):
            shutil.rmtree(stage_root, ignore_errors=True)

        return summary


# Backward-compatible name used by earlier CLI versions.
TemplateOnlyPipeline = TApocketPipeline
