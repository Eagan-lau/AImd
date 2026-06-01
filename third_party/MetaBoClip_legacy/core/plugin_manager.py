from __future__ import annotations

import sys
from pathlib import Path

from core.plugins import YamlFamilyPlugin, FamilyPlugin


class PluginManager:
    def __init__(self, plugin_dirs=None):
        self.project_root = Path(__file__).resolve().parent.parent
        if str(self.project_root) not in sys.path:
            sys.path.insert(0, str(self.project_root))
        self.plugin_dirs = [Path(p).resolve() for p in (plugin_dirs or [])]
        self.plugins: dict[str, FamilyPlugin] = {}
        self._discover_builtin_plugins()
        self._discover_external_plugins()

    def _register_spec(self, spec_path: Path, origin: str):
        plugin = YamlFamilyPlugin(spec_path=spec_path, origin=origin)
        self.plugins[plugin.key] = plugin

    def _discover_builtin_plugins(self):
        families_root = self.project_root / "families"
        if not families_root.exists():
            return
        for spec_file in sorted(families_root.glob("*/family.yaml")):
            self._register_spec(spec_file, origin="builtin")
        for spec_file in sorted(families_root.glob("*.yaml")):
            self._register_spec(spec_file, origin="builtin")

    def _iter_external_spec_files(self):
        for root in self.plugin_dirs:
            if not root.exists():
                continue
            if root.is_file() and root.suffix.lower() in {".yaml", ".yml"}:
                yield root
                continue
            direct = root / "family.yaml"
            if direct.exists():
                yield direct
            direct2 = root / "default.yaml"
            if direct2.exists():
                yield direct2
            for spec_file in sorted(root.glob("*/family.yaml")):
                yield spec_file
            for spec_file in sorted(root.glob("*.yaml")):
                yield spec_file

    def _discover_external_plugins(self):
        for spec_file in self._iter_external_spec_files():
            self._register_spec(spec_file, origin="external")

    def get(self, key: str) -> FamilyPlugin:
        norm = str(key).strip().lower().replace("-", "_").replace(" ", "_")
        if norm not in self.plugins:
            raise RuntimeError(f"Unknown family plugin '{key}'. Available: {sorted(self.plugins)}")
        return self.plugins[norm]

    def describe_all(self):
        return [self.plugins[k].describe() for k in sorted(self.plugins)]
