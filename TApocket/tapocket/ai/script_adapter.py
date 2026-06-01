from __future__ import annotations

import json
import shlex
import subprocess
import time
from pathlib import Path
from typing import Any

from tapocket.ai.parsers import locate_ai_output, parse_tapocket_json
from tapocket.core.schema import CandidatePocket


class ScriptAIPredictor:
    """Generic script adapter for DeepPocket-DB or any other AI pocket model.

    The external script is responsible for inference. TApocket only runs the command
    and parses a standard JSON output file.
    """

    def __init__(self, config: Any, predictor_cfg: dict[str, Any], query_id: str, run_root: str | Path):
        self.config = config
        self.predictor_cfg = predictor_cfg
        self.query_id = query_id
        self.run_root = Path(run_root).resolve()
        self.name = str(predictor_cfg.get("name", "ai_model"))
        self.device = str(predictor_cfg.get("device", config.get("ai_models", "device", default="cuda")))

    def _format_command(self, query_pdb: str | Path, output_dir: str | Path, run_id: str | None) -> list[str]:
        command = self.predictor_cfg.get("command", [])
        if isinstance(command, str):
            command = shlex.split(command)
        values = {
            "query_pdb": str(Path(query_pdb).resolve()),
            "output_dir": str(Path(output_dir).resolve()),
            "device": self.device,
            "run_id": run_id or "",
            "project_root": str(self.config.root),
        }
        return [str(part).format(**values) for part in command]

    def predict(self, query_pdb: str | Path, output_dir: str | Path, run_id: str | None = None) -> list[CandidatePocket]:
        output_dir = Path(output_dir).resolve()
        output_dir.mkdir(parents=True, exist_ok=True)
        command = self._format_command(query_pdb, output_dir, run_id)
        if not command:
            raise RuntimeError(f"AI predictor {self.name} has an empty command")

        started = time.time()
        completed = subprocess.run(command, cwd=str(self.config.root), capture_output=True, text=True)
        elapsed = time.time() - started
        log = {
            "predictor": self.name,
            "command": command,
            "returncode": completed.returncode,
            "stdout": completed.stdout,
            "stderr": completed.stderr,
            "elapsed_seconds": elapsed,
        }
        (output_dir / "ai_command_log.json").write_text(json.dumps(log, indent=2, ensure_ascii=False), encoding="utf-8")
        if completed.returncode != 0:
            raise RuntimeError(
                f"AI predictor {self.name} failed.\nCommand: {' '.join(command)}\nSTDERR:\n{completed.stderr}"
            )

        parser = str(self.predictor_cfg.get("output_parser", "tapocket_json"))
        output_json = locate_ai_output(output_dir, self.predictor_cfg.get("output_json"))
        if parser != "tapocket_json":
            raise RuntimeError(
                f"Unsupported output_parser={parser!r}. Add a parser in tapocket/ai/parsers.py or use tapocket_json."
            )
        return parse_tapocket_json(output_json, self.query_id, self.name, self.run_root, output_dir)
