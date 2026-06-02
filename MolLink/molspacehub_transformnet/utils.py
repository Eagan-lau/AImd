from __future__ import annotations

import json
from pathlib import Path
from typing import Any

def ensure_dir(path: str | Path) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p

def write_json(obj: dict[str, Any], path: str | Path) -> Path:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w", encoding="utf-8") as handle:
        json.dump(obj, handle, ensure_ascii=False, indent=2)
    return p

def semicolon_join(values, max_items: int | None = None) -> str:
    seen = []
    for v in values:
        if v is None:
            continue
        s = str(v)
        if s == "" or s.lower() == "nan":
            continue
        if s not in seen:
            seen.append(s)
    if max_items is not None and len(seen) > max_items:
        return ";".join(seen[:max_items]) + f";...(+{len(seen)-max_items})"
    return ";".join(seen)

def bool_to_int(v: bool) -> int:
    return 1 if bool(v) else 0
