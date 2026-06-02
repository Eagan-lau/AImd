from __future__ import annotations

import json
import shlex
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

from tapocket.ai.parsers import locate_ai_output, parse_tapocket_json
from tapocket.core.schema import CandidatePocket


VARIANT_SPECS: dict[str, dict[str, Any]] = {
    "basic31": {
        "input_dim": 31,
        "feature_type": "basic31",
        "require_pocketminer": False,
        "pocketminer_mode": None,
    },
    "pocketminer_score": {
        "input_dim": 32,
        "feature_type": "basic31_plus_pocketminer_score",
        "require_pocketminer": True,
        "pocketminer_mode": "score",
    },
    "pocketminer_hidden": {
        "input_dim": 131,
        "feature_type": "basic31_plus_pocketminer_hidden",
        "require_pocketminer": True,
        "pocketminer_mode": "hidden",
        "pocketminer_hidden_dim": 100,
    },
}


class DeepPocketDBPredictor:
    """Native adapter for the user's DeepPocket-DB models.

    This adapter is intentionally still script-based at the execution boundary so that
    TApocket and DeepPocket-DB can keep separate environments. The native part is the
    configuration validation and automatic handling of three DeepPocket-DB variants:

    - basic31: 31-D TApocket residue/geometry features only.
    - pocketminer_score: basic31 + PocketMiner residue score, input dim 32.
    - pocketminer_hidden: basic31 + PocketMiner hidden vector, input dim 131.

    The DeepPocket-DB script should output a TApocket-compatible JSON file, normally
    named tapocket_ai_pockets.json. If your script uses a different raw output, add a
    parser in tapocket/ai/parsers.py and set inference.output_parser accordingly.
    """

    def __init__(self, config: Any, predictor_cfg: dict[str, Any], query_id: str, run_root: str | Path):
        self.config = config
        self.predictor_cfg = predictor_cfg
        self.query_id = query_id
        self.run_root = Path(run_root).resolve()
        self.name = str(predictor_cfg.get("name", "deeppocket_db"))
        self.variant = str(predictor_cfg.get("variant", "basic31"))
        self.device = str(predictor_cfg.get("device", config.get("ai_models", "device", default="cuda")))
        self.root = Path(config.root).resolve()
        self.spec = self._validate_variant_config()

    def _root_path(self, value: str | Path | None, *, required: bool = False, field_name: str = "path") -> Path | None:
        if value is None or str(value) == "":
            if required:
                raise FileNotFoundError(f"Missing required DeepPocket-DB field: {field_name}")
            return None
        path = Path(str(value))
        if not path.is_absolute():
            path = self.root / path
        return path.resolve()

    def _validate_variant_config(self) -> dict[str, Any]:
        if self.variant not in VARIANT_SPECS:
            allowed = ", ".join(sorted(VARIANT_SPECS))
            raise ValueError(f"Unsupported DeepPocket-DB variant={self.variant!r}. Allowed: {allowed}")

        spec = dict(VARIANT_SPECS[self.variant])
        input_dim = int(self.predictor_cfg.get("input_dim", spec["input_dim"]))
        if input_dim != int(spec["input_dim"]):
            raise ValueError(
                f"DeepPocket-DB variant={self.variant} requires input_dim={spec['input_dim']}, "
                f"but config has input_dim={input_dim}."
            )
        spec["input_dim"] = input_dim

        feature_cfg = self.predictor_cfg.get("feature", {}) or {}
        require_pm = bool(feature_cfg.get("require_pocketminer", spec["require_pocketminer"]))
        if require_pm != bool(spec["require_pocketminer"]):
            raise ValueError(
                f"DeepPocket-DB variant={self.variant} has require_pocketminer={spec['require_pocketminer']}, "
                f"but config has require_pocketminer={require_pm}."
            )

        checkpoint = self._root_path(self.predictor_cfg.get("checkpoint"), required=True, field_name="checkpoint")
        if checkpoint and not checkpoint.exists():
            raise FileNotFoundError(f"DeepPocket-DB checkpoint not found: {checkpoint}")

        model_config = self._root_path(self.predictor_cfg.get("model_config"), required=False, field_name="model_config")
        if model_config and not model_config.exists():
            raise FileNotFoundError(f"DeepPocket-DB model_config not found: {model_config}")

        inference_cfg = self.predictor_cfg.get("inference", {}) or {}
        script = self._root_path(inference_cfg.get("script"), required=True, field_name="inference.script")
        if script and not script.exists():
            raise FileNotFoundError(f"DeepPocket-DB inference script not found: {script}")

        cache_dir = self._root_path(feature_cfg.get("cache_dir"), required=False, field_name="feature.cache_dir")
        if cache_dir:
            cache_dir.mkdir(parents=True, exist_ok=True)

        if spec["require_pocketminer"]:
            pm_root = self._root_path(feature_cfg.get("pocketminer_root"), required=True, field_name="feature.pocketminer_root")
            pm_ckpt = self._root_path(feature_cfg.get("pocketminer_checkpoint"), required=True, field_name="feature.pocketminer_checkpoint")
            if pm_root and not pm_root.exists():
                raise FileNotFoundError(f"PocketMiner root not found: {pm_root}")
            if pm_ckpt and not pm_ckpt.exists():
                raise FileNotFoundError(f"PocketMiner checkpoint not found: {pm_ckpt}")
            pm_mode = str(feature_cfg.get("pocketminer_mode", spec["pocketminer_mode"]))
            if pm_mode != spec["pocketminer_mode"]:
                raise ValueError(
                    f"DeepPocket-DB variant={self.variant} requires pocketminer_mode={spec['pocketminer_mode']}, "
                    f"but config has pocketminer_mode={pm_mode}."
                )
            spec["pocketminer_root"] = str(pm_root)
            spec["pocketminer_checkpoint"] = str(pm_ckpt)
            spec["pocketminer_mode"] = pm_mode

        spec["checkpoint"] = str(checkpoint)
        spec["model_config"] = str(model_config) if model_config else ""
        spec["script"] = str(script)
        spec["feature_type"] = str(feature_cfg.get("type", spec["feature_type"]))
        spec["cache_dir"] = str(cache_dir) if cache_dir else ""
        return spec

    def _format_sequence(
        self,
        sequence: str | list[Any] | None,
        values: dict[str, str],
    ) -> list[str]:
        if not sequence:
            return []
        if isinstance(sequence, str):
            parts = shlex.split(sequence)
        else:
            parts = [str(x) for x in sequence]
        return [part.format(**values) for part in parts]

    def _values(self, query_pdb: str | Path, output_dir: str | Path, run_id: str | None) -> dict[str, str]:
        output_dir = Path(output_dir).resolve()
        return {
            "query_pdb": str(Path(query_pdb).resolve()),
            "output_dir": str(output_dir),
            "device": self.device,
            "run_id": run_id or "",
            "project_root": str(self.root),
            "checkpoint": self.spec["checkpoint"],
            "model_config": self.spec.get("model_config", ""),
            "feature_type": self.spec["feature_type"],
            "input_dim": str(self.spec["input_dim"]),
            "variant": self.variant,
            "feature_cache_dir": self.spec.get("cache_dir", ""),
            "pocketminer_root": self.spec.get("pocketminer_root", ""),
            "pocketminer_checkpoint": self.spec.get("pocketminer_checkpoint", ""),
            "pocketminer_mode": str(self.spec.get("pocketminer_mode") or ""),
            "pocketminer_hidden_dim": str(self.spec.get("pocketminer_hidden_dim", "")),
            "python_executable": sys.executable,
        }

    def _default_inference_command(self, values: dict[str, str]) -> list[str]:
        command = [
            values["python_executable"],
            values["script"],
            "--input-pdb",
            values["query_pdb"],
            "--checkpoint",
            values["checkpoint"],
            "--out-dir",
            values["output_dir"],
            "--device",
            values["device"],
            "--feature-type",
            values["feature_type"],
            "--input-dim",
            values["input_dim"],
        ]
        if values.get("model_config"):
            command += ["--model-config", values["model_config"]]
        if values.get("feature_cache_dir"):
            command += ["--feature-cache-dir", values["feature_cache_dir"]]
        if self.spec["require_pocketminer"]:
            command += [
                "--pocketminer-root",
                values["pocketminer_root"],
                "--pocketminer-checkpoint",
                values["pocketminer_checkpoint"],
                "--pocketminer-output",
                values["pocketminer_mode"],
            ]
            if values.get("pocketminer_hidden_dim"):
                command += ["--pocketminer-hidden-dim", values["pocketminer_hidden_dim"]]
        extra_args = self.predictor_cfg.get("inference", {}).get("extra_args", []) or []
        command.extend(self._format_sequence(extra_args, values))
        return command

    @staticmethod
    def _timeout_seconds(stage_cfg: dict[str, Any]) -> int | None:
        value = stage_cfg.get("timeout")
        if value in {None, "", 0, "0"}:
            return None
        return int(value)

    @staticmethod
    def _timeout_output(value: Any) -> str:
        if value is None:
            return ""
        if isinstance(value, bytes):
            return value.decode(errors="replace")
        return str(value)

    def _run_command(self, command: list[str], output_dir: Path, label: str, timeout: int | None = None) -> dict[str, Any]:
        started = time.time()
        try:
            completed = subprocess.run(command, cwd=str(self.root), capture_output=True, text=True, timeout=timeout)
        except subprocess.TimeoutExpired as exc:
            elapsed = time.time() - started
            log = {
                "predictor": self.name,
                "variant": self.variant,
                "stage": label,
                "command": command,
                "returncode": "timeout",
                "stdout": self._timeout_output(exc.stdout),
                "stderr": self._timeout_output(exc.stderr),
                "elapsed_seconds": elapsed,
                "timeout_seconds": timeout,
            }
            (output_dir / f"{label}_command_log.json").write_text(json.dumps(log, indent=2, ensure_ascii=False), encoding="utf-8")
            raise TimeoutError(
                f"DeepPocket-DB {self.name} {label} timed out after {timeout}s. "
                f"Command: {' '.join(command)}"
            ) from exc
        elapsed = time.time() - started
        log = {
            "predictor": self.name,
            "variant": self.variant,
            "stage": label,
            "command": command,
            "returncode": completed.returncode,
            "stdout": completed.stdout,
            "stderr": completed.stderr,
            "elapsed_seconds": elapsed,
            "timeout_seconds": timeout,
        }
        (output_dir / f"{label}_command_log.json").write_text(json.dumps(log, indent=2, ensure_ascii=False), encoding="utf-8")
        if completed.returncode != 0:
            raise RuntimeError(
                f"DeepPocket-DB {self.name} {label} failed.\n"
                f"Command: {' '.join(command)}\nSTDERR:\n{completed.stderr}"
            )
        return log

    def predict(self, query_pdb: str | Path, output_dir: str | Path, run_id: str | None = None) -> list[CandidatePocket]:
        output_dir = Path(output_dir).resolve()
        output_dir.mkdir(parents=True, exist_ok=True)
        values = self._values(query_pdb, output_dir, run_id)
        values["script"] = self.spec["script"]

        feature_cfg = self.predictor_cfg.get("feature", {}) or {}
        precompute_command = self._format_sequence(feature_cfg.get("command"), values)
        command_logs: list[dict[str, Any]] = []
        if precompute_command:
            command_logs.append(self._run_command(precompute_command, output_dir, "feature_precompute", timeout=self._timeout_seconds(feature_cfg)))

        inference_cfg = self.predictor_cfg.get("inference", {}) or {}
        inference_command = self._format_sequence(inference_cfg.get("command"), values)
        if not inference_command:
            inference_command = self._default_inference_command(values)
        command_logs.append(self._run_command(inference_command, output_dir, "inference", timeout=self._timeout_seconds(inference_cfg)))

        meta = {
            "predictor": self.name,
            "type": "deeppocket_db",
            "variant": self.variant,
            "input_dim": self.spec["input_dim"],
            "feature_type": self.spec["feature_type"],
            "require_pocketminer": self.spec["require_pocketminer"],
            "command_logs": command_logs,
        }
        (output_dir / "deeppocket_db_adapter_meta.json").write_text(json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")

        parser = str(inference_cfg.get("output_parser", self.predictor_cfg.get("output_parser", "tapocket_json")))
        output_json_value = inference_cfg.get("output_json", self.predictor_cfg.get("output_json"))
        output_json = locate_ai_output(output_dir, output_json_value)
        if parser != "tapocket_json":
            raise RuntimeError(
                f"Unsupported DeepPocket-DB output_parser={parser!r}. Add a parser in tapocket/ai/parsers.py "
                "or use tapocket_json."
            )
        pockets = parse_tapocket_json(output_json, self.query_id, self.name, self.run_root, output_dir)
        for pocket in pockets:
            pocket.ai_support.setdefault("model_name", self.name)
            pocket.ai_support["adapter_type"] = "deeppocket_db"
            pocket.ai_support["variant"] = self.variant
            pocket.ai_support["feature_type"] = self.spec["feature_type"]
            pocket.extra.setdefault("deeppocket_db", {})
            pocket.extra["deeppocket_db"].update({
                "variant": self.variant,
                "input_dim": self.spec["input_dim"],
                "feature_type": self.spec["feature_type"],
            })
        return pockets
