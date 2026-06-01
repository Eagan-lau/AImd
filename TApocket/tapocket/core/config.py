from __future__ import annotations

import os
from pathlib import Path
from typing import Any

try:
    import yaml
except ImportError as exc:  # pragma: no cover
    raise RuntimeError("PyYAML is required. Install with: pip install PyYAML") from exc


class TapocketConfig:
    def __init__(self, data: dict[str, Any], config_path: Path):
        self.data = data
        self.config_path = Path(config_path).resolve()
        self.root = self._detect_root()

    def _detect_root(self) -> Path:
        env_root = os.environ.get("TAPOCKET_ROOT")
        if env_root:
            return Path(env_root).resolve()

        root_value = str(self.data.get("project", {}).get("root", "."))
        root_path = Path(root_value)
        if root_path.is_absolute():
            return root_path.resolve()

        if root_value == ".":
            if self.config_path.parent.name == "configs":
                return self.config_path.parent.parent.resolve()
            return Path.cwd().resolve()

        return (self.config_path.parent / root_path).resolve()

    def get(self, *keys: str, default: Any = None) -> Any:
        current: Any = self.data
        for key in keys:
            if not isinstance(current, dict) or key not in current:
                return default
            current = current[key]
        return current

    def path(self, *keys: str, default: str | None = None) -> Path:
        value = self.get(*keys, default=default)
        if value is None:
            raise KeyError(f"Missing path config: {'.'.join(keys)}")
        p = Path(str(value))
        if p.is_absolute():
            return p
        return (self.root / p).resolve()

    def relpath(self, path: str | Path) -> str:
        p = Path(path).resolve()
        try:
            return str(p.relative_to(self.root))
        except ValueError:
            return str(p)


def load_config(config_path: str | Path) -> TapocketConfig:
    config_path = Path(config_path).resolve()
    with config_path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    return TapocketConfig(data, config_path)
