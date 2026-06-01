import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import tempfile
import textwrap
import unittest
from pathlib import Path

from core.plugin_manager import PluginManager


class TestExternalYamlSpecs(unittest.TestCase):
    def test_external_family_yaml_is_discovered(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            fam = root / "myfam"
            fam.mkdir(parents=True)
            (fam / "family.yaml").write_text(
                textwrap.dedent(
                    """
                    version: 1
                    family:
                      key: myfam
                      plugin: myfam
                      runner: myfam
                      name: MyFam
                    input:
                      protein_extension: .pdbqt
                      ligand_extension: .pdbqt
                      ligand_pattern: '{ligand_id}@{protein_name}{ligand_extension}'
                      out_pattern: '{ligand_id}_{protein_name}_cavity_1.out'
                      protein_object: protein
                      ligand_object: ligand
                    aggregation:
                      total_confs: 6
                      cover_t: 70.0
                      alpha: 0.3
                    output:
                      save_pse: false
                    resources: {}
                    globals: {}
                    pose:
                      add_hydrogens: false
                      objects: {}
                    candidates: {}
                    features: []
                    gating:
                      logic: all
                      rules: []
                    scoring:
                      transforms: []
                      intermediates: []
                      pose_score:
                        op: weighted_sum
                        scale: 100.0
                        terms: []
                    """
                ).strip()
                + "\n",
                encoding="utf-8",
            )

            pm = PluginManager(plugin_dirs=[root])
            self.assertIn("myfam", pm.plugins)
            cfg = pm.get("myfam").load_default_config()
            self.assertEqual(cfg["family"]["plugin"], "myfam")

    def test_external_spec_overrides_builtin_key(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            fam = root / "override_ugt"
            fam.mkdir(parents=True)
            (fam / "family.yaml").write_text(
                textwrap.dedent(
                    """
                    version: 1
                    family:
                      key: ugt
                      plugin: ugt
                      runner: ugt
                      name: UGT override
                      version: 9.9.9
                    input:
                      protein_extension: .pdbqt
                      ligand_extension: .pdbqt
                      ligand_pattern: '{ligand_id}@{protein_name}{ligand_extension}'
                      out_pattern: '{ligand_id}_{protein_name}_cavity_1.out'
                      protein_object: protein
                      ligand_object: ligand
                    aggregation:
                      total_confs: 6
                      cover_t: 70.0
                      alpha: 0.3
                    output:
                      save_pse: false
                    resources: {}
                    globals: {}
                    pose:
                      add_hydrogens: false
                      objects: {}
                    candidates: {}
                    features: []
                    gating:
                      logic: all
                      rules: []
                    scoring:
                      transforms: []
                      intermediates: []
                      pose_score:
                        op: weighted_sum
                        scale: 100.0
                        terms: []
                    """
                ).strip()
                + "\n",
                encoding="utf-8",
            )

            pm = PluginManager(plugin_dirs=[root])
            self.assertEqual(pm.get("ugt").version, "9.9.9")
            self.assertEqual(pm.get("ugt").origin, "external")
