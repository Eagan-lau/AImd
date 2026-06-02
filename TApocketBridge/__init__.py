"""AImd adapter for batch TApocket pocket prediction.

This adapter connects RGPC protein batches (`data/data_output/protein_batches/file_*`) to
TApocket single-query runs and exports standardized AImd pocket manifests.
"""

from .runner import run_tapocket_batch

__all__ = ["run_tapocket_batch"]
