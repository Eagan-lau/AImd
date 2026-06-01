from __future__ import annotations

import json
import platform
import subprocess
import sys
import time
from pathlib import Path
from typing import Any


def _command_version(command: list[str]) -> str | None:
    try:
        result = subprocess.run(command, capture_output=True, text=True, timeout=10)
    except Exception:
        return None
    text = (result.stdout or result.stderr).strip()
    return text.splitlines()[0] if text else None


def write_provenance(
    path: str | Path,
    *,
    config_path: str | Path,
    project_root: str | Path,
    run_id: str,
    query_path: str | Path,
    extra: dict[str, Any] | None = None,
) -> None:
    data = {
        "run_id": run_id,
        "created_at_epoch": time.time(),
        "project_root": str(Path(project_root).resolve()),
        "config_path": str(Path(config_path).resolve()),
        "query_path": str(Path(query_path).resolve()),
        "python": sys.version,
        "platform": platform.platform(),
        "foldseek_version": _command_version(["foldseek", "version"]),
        "extra": extra or {},
    }
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8")
