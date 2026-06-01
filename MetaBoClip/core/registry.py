from engine.operators import OBJECT_OPS, FEATURE_OPS
from engine.scoring import SCORING_OPS


def build_default_registry():
    return {
        "object_ops": sorted(OBJECT_OPS),
        "feature_ops": sorted(FEATURE_OPS),
        "scoring_ops": sorted(SCORING_OPS),
    }
