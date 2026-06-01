REQUIRED_TOP_LEVEL = [
    "family",
    "input",
    "aggregation",
    "output",
    "globals",
    "pose",
    "candidates",
    "features",
    "gating",
    "scoring",
]


def validate_config(cfg: dict):
    missing = [k for k in REQUIRED_TOP_LEVEL if k not in cfg]
    if missing:
        raise RuntimeError(f"Missing required top-level keys: {missing}")
    fam = cfg.get("family", {})
    if "name" not in fam:
        raise RuntimeError("family.name is required")
    if "key" not in fam and "plugin" not in fam:
        raise RuntimeError("family.key or family.plugin is required")
    return True
