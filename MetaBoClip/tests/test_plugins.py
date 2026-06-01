import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import unittest

from core.plugin_manager import PluginManager
from core.runtime import resolve_plugin_and_config


class TestBuiltins(unittest.TestCase):
    def test_builtin_specs_discovered(self):
        pm = PluginManager()
        self.assertEqual(set(pm.plugins.keys()), {"ugt", "act", "cyp450", "fe2og", "ach"})

    def test_default_spec_loads_and_is_normalized(self):
        pm = PluginManager()
        cfg = pm.get("ugt").load_default_config()
        self.assertEqual(cfg["family"]["plugin"], "ugt")
        self.assertIn("globals", cfg)
        self.assertIn("features", cfg)
        self.assertEqual(cfg["family"]["version"], '2.0.0')

    def test_overlay_merges_into_family_default(self):
        raw_cfg = {
            "family": {"runner": "act"},
            "input": {"affinity_threshold": -4.5},
        }
        _, plugin, merged, _ = resolve_plugin_and_config(raw_cfg)
        self.assertEqual(plugin.key, "act")
        self.assertEqual(merged["family"]["plugin"], "act")
        self.assertEqual(merged["input"]["affinity_threshold"], -4.5)
        self.assertIn("globals", merged)
