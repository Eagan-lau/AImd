#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import annotations
import shlex
import shutil
import subprocess
from pathlib import Path
from typing import Any
try:
    import yaml
except ImportError as exc:
    raise RuntimeError("PyYAML is required. Install with: pip install PyYAML") from exc

def default_tools_config(aimd_root: str | Path | None = None) -> Path:
    return Path(aimd_root or ".").resolve() / "third_party" / "tools.yaml"

def load_tools_config(config_path: str | Path | None = None, aimd_root: str | Path | None = None) -> dict[str, Any]:
    root = Path(aimd_root or ".").resolve()
    path = Path(config_path).expanduser() if config_path else default_tools_config(root)
    if not path.is_absolute():
        path = root / path
    if not path.exists():
        return {"tools": {}}
    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {"tools": {}}

def _resolve_executable_path(executable: str, aimd_root: str | Path | None = None) -> str:
    executable = str(executable).strip()
    if not executable:
        return executable
    if any(ch.isspace() for ch in executable):
        parts = shlex.split(executable)
        if not parts:
            return executable
        first = parts[0]
        resolved_first = _resolve_executable_path(first, aimd_root)
        return executable.replace(first, resolved_first, 1)
    p = Path(executable).expanduser()
    if p.is_absolute() and p.exists():
        return str(p)
    root = Path(aimd_root or ".").resolve()
    if not p.is_absolute():
        candidate = root / p
        if candidate.exists():
            return str(candidate)
        local_bin = root / "third_party" / "bin" / executable
        if local_bin.exists():
            return str(local_bin)
    found = shutil.which(executable)
    return found or executable

def resolve_tool(tool_name: str, configured_executable: str | None = None, aimd_root: str | Path | None = None, tools_config: str | Path | None = None) -> str:
    root = Path(aimd_root or ".").resolve()
    if configured_executable and str(configured_executable).strip().lower() not in {"", "auto", "null", "none"}:
        return _resolve_executable_path(str(configured_executable), root)
    cfg = load_tools_config(tools_config, root)
    tool_cfg = cfg.get("tools", {}).get(tool_name, {}) or {}
    executable = str(tool_cfg.get("executable", tool_name))
    return _resolve_executable_path(executable, root)

def check_tool(tool_name: str, aimd_root: str | Path | None = None, tools_config: str | Path | None = None) -> dict[str, Any]:
    root = Path(aimd_root or ".").resolve()
    cfg = load_tools_config(tools_config, root)
    tool_cfg = cfg.get("tools", {}).get(tool_name, {}) or {}
    executable = resolve_tool(tool_name, tool_cfg.get("executable"), root, tools_config)
    enabled = bool(tool_cfg.get("enabled", True))
    required = bool(tool_cfg.get("required", False))
    if not enabled:
        return {"tool": tool_name, "enabled": enabled, "required": required, "executable": executable, "status": "disabled", "return_code": "", "message": "disabled"}
    first = shlex.split(executable)[0] if executable else tool_name
    exists = Path(first).exists() or shutil.which(first) is not None
    if not exists:
        return {"tool": tool_name, "enabled": enabled, "required": required, "executable": executable, "status": "missing", "return_code": "127", "message": "executable not found"}
    check_command = str(tool_cfg.get("check_command", "")).strip()
    if not check_command:
        return {"tool": tool_name, "enabled": enabled, "required": required, "executable": executable, "status": "ok", "return_code": "", "message": "executable found"}
    command = check_command.replace("{executable}", shlex.quote(executable))
    try:
        proc = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=int(tool_cfg.get("check_timeout", 20)))
        out = (proc.stdout or "").strip().split("\n")[0][:300]
        return {"tool": tool_name, "enabled": enabled, "required": required, "executable": executable, "status": "ok" if proc.returncode == 0 else "failed", "return_code": str(proc.returncode), "message": out}
    except subprocess.TimeoutExpired:
        return {"tool": tool_name, "enabled": enabled, "required": required, "executable": executable, "status": "timeout", "return_code": "124", "message": "check command timed out"}
