from __future__ import annotations

import importlib
from pathlib import Path
from typing import Callable

from core.config import load_config, normalize_family_key, normalize_family_section
from engine.executor import run_family_spec


class FamilyPlugin:
    def __init__(
        self,
        key: str,
        display_name: str,
        version: str = "1.0.0",
        description: str = "",
        default_config_name: str = "family.yaml",
    ):
        self.key = normalize_family_key(key) or key
        self.display_name = display_name
        self.version = version
        self.description = description
        self.default_config_name = default_config_name
        self.module_path: Path | None = None
        self.origin: str = "unknown"

    def bind_origin(self, module_path: str | Path | None, origin: str = "builtin") -> "FamilyPlugin":
        self.module_path = Path(module_path).resolve() if module_path else None
        self.origin = origin
        return self

    def resolve_default_config_path(self) -> Path | None:
        if self.module_path is None:
            return None
        if self.module_path.is_file():
            return self.module_path
        return self.module_path / self.default_config_name

    def load_default_config(self) -> dict:
        path = self.resolve_default_config_path()
        cfg = load_config(path) if path and path.exists() else {}
        return normalize_family_section(cfg, fallback_key=self.key, fallback_name=self.display_name)

    def describe(self) -> dict:
        return {
            "key": self.key,
            "display_name": self.display_name,
            "version": self.version,
            "description": self.description,
            "origin": self.origin,
            "module_path": str(self.module_path) if self.module_path else None,
            "default_config_path": str(self.resolve_default_config_path()) if self.resolve_default_config_path() else None,
        }

    def execute(self, *, cfg: dict, protein_dir: Path, docking_dir: Path, out_dir: Path, registry: dict) -> dict:
        raise NotImplementedError


class YamlFamilyPlugin(FamilyPlugin):
    def __init__(self, spec_path: str | Path, origin: str = "builtin"):
        spec_path = Path(spec_path).resolve()
        cfg = load_config(spec_path)
        family = dict(cfg.get("family") or {})
        key = family.get("key") or family.get("plugin") or family.get("runner") or family.get("name") or spec_path.parent.name
        display_name = family.get("name") or str(key)
        version = str(family.get("version", "1.0.0"))
        description = family.get("description", "")
        super().__init__(key=key, display_name=display_name, version=version, description=description, default_config_name=spec_path.name)
        self.bind_origin(spec_path, origin=origin)

    def execute(self, *, cfg: dict, protein_dir: Path, docking_dir: Path, out_dir: Path, registry: dict) -> dict:
        return run_family_spec(cfg=cfg, protein_dir=protein_dir, docking_dir=docking_dir, out_dir=out_dir, spec_path=self.resolve_default_config_path(), registry=registry)


class RunnerBackedFamilyPlugin(FamilyPlugin):
    def __init__(
        self,
        key: str,
        display_name: str,
        version: str = "1.0.0",
        description: str = "",
        runner_import: str | None = None,
        runner: Callable | None = None,
        default_config_name: str = "default.yaml",
    ):
        super().__init__(key=key, display_name=display_name, version=version, description=description, default_config_name=default_config_name)
        self.runner_import = runner_import
        self.runner = runner

    def _load_runner(self) -> Callable:
        if callable(self.runner):
            return self.runner
        if not self.runner_import:
            raise RuntimeError(f"Family plugin '{self.key}' does not define a runner")
        mod_name, attr = self.runner_import.split(":", 1)
        module = importlib.import_module(mod_name)
        fn = getattr(module, attr)
        if not callable(fn):
            raise RuntimeError(f"Runner import '{self.runner_import}' for family '{self.key}' is not callable")
        return fn

    def execute(self, *, cfg: dict, protein_dir: Path, docking_dir: Path, out_dir: Path, registry: dict) -> dict:
        runner = self._load_runner()
        return runner(cfg=cfg, protein_dir=protein_dir, docking_dir=docking_dir, out_dir=out_dir, registry=registry)
