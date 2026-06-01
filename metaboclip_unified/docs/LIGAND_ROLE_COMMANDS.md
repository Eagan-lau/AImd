# LigandRoleMap Commands

This patch adds package modules for LigandRoleMap command-line tools.

## Direct module usage

```bash
python -m metaboclip_ligand_roles.build_role_table --help
python -m metaboclip_ligand_roles.build_role_tables_batch --help
python -m metaboclip_ligand_roles.recover_atom_map --help
python -m metaboclip_ligand_roles.extract_role_pose_coords --help
python -m metaboclip_ligand_roles.list_atom_labels --help
```

## Console script usage

Run the pyproject patch helper once from the project root:

```bash
python scripts/patch_pyproject_ligand_role_scripts.py
pip install -e .
hash -r
```

Then these commands should be available:

```bash
metaboclip-build-role-table --help
metaboclip-build-role-tables-batch --help
metaboclip-recover-atom-map --help
metaboclip-extract-role-coords --help
metaboclip-list-atom-labels --help
```

## Single-ligand role table

```bash
metaboclip-build-role-table \
  --ligand-source data_input/ligand_raw/1.mol2 \
  --prepared-pdbqt data_input/ligand_pdbqt/1.pdbqt \
  --ligand-id 1 \
  --rules rules/functional_groups.yaml \
  --out-role-table data_output/ligand_roles/1.role_table.csv \
  --out-annotation data_output/ligand_roles/1.annotation.json \
  --out-atom-map data_output/ligand_roles/1.atom_map.json
```
