#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Minimal adapter contract for future AImd tool plug-ins.

Current modules are already modular and config-driven. This file is intentionally
small: it documents the expected adapter surface for future replacements such as
P2Rank/FPocket, Gnina/Uni-Dock, Leiden/MCL, or alternative ensemble samplers.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any


class BaseAImdAdapter(ABC):
    """Base interface for future external-tool adapters."""

    name: str = "base"
    task_type: str = "generic"

    @abstractmethod
    def validate_inputs(self, config: dict[str, Any]) -> None:
        """Validate required input files, directories, and executable availability."""

    @abstractmethod
    def run(self, config: dict[str, Any]) -> Path | dict[str, Any] | list[dict[str, Any]]:
        """Execute the tool and return its primary output path or manifest rows."""

    @abstractmethod
    def collect_outputs(self, config: dict[str, Any]) -> Path | dict[str, Any] | list[dict[str, Any]]:
        """Convert native outputs into an AImd-standard manifest."""
