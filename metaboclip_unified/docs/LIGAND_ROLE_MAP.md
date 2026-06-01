# LigandRoleMap

LigandRoleMap converts original ligand chemistry into stable cross-molecule atom labels and maps these labels to ligand PDBQT atom order.

## Purpose

LigandRoleMap answers:

```text
Which ligand atoms can mechanism.yaml refer to?
What is the stable atom_label of each relevant atom?
What is the PDBQT atom order used to retrieve pose coordinates?
```

## Key rule

Use original SDF/MOL/MOL2/SMILES for chemistry. Use PDBQT only for coordinate mapping.

Hydrogen atoms may be used for functional group recognition, but generic catalytic geometry uses heavy atoms only.

Examples:

```text
hydroxyl.o
  H confirms hydroxyl identity.
  O is used for geometry.

c_h_site.c
  H confirms the carbon has hydrogen.
  C is used for geometry.
```

## Main output: role_table.csv

Required columns:

```text
ligand_id
group_id
instance_id
atom_label
atom_class
atom_role
source_atom_index
element
pdbqt_order
subtype
confidence
```

Example rows:

```text
1,hydroxyl,hydroxyl_1,hydroxyl.o,hydroxyl_oxygen,o,22,O,23,aliphatic_hydroxyl,1.0
1,protic_nucleophile,protic_nucleophile_1,protic_nucleophile.atom,protic_heteroatom,atom,22,O,23,,1.0
1,acetyl,acetyl_1,acetyl.acyl_c,acyl_carbonyl_c,acyl_c,31,C,19,,1.0
1,glycosyl_pyranose,glycosyl_pyranose_1,glycosyl_pyranose.c1,glycosyl_anomeric_c,c1,34,C,22,hexopyranosyl_like,1.0
```

## Generate a single role table

```bash
metaboclip-build-role-table \
  --ligand-source /path/to/original_ligands/1.sdf \
  --prepared-pdbqt /path/to/prepared_ligands/1.pdbqt \
  --ligand-id 1 \
  --rules rules/functional_groups.yaml \
  --out-role-table /path/to/ligand_roles/1.role_table.csv \
  --out-annotation /path/to/ligand_roles/1.annotation.json \
  --out-atom-map /path/to/ligand_roles/1.atom_map.json
```

## Generate role tables in batch

```bash
metaboclip-build-role-tables-batch \
  --ligand-source-dir /path/to/original_ligands \
  --prepared-pdbqt-dir /path/to/prepared_ligands \
  --rules rules/functional_groups.yaml \
  --out-dir /path/to/ligand_roles
```

## Inspect atom labels

```bash
cut -d, -f4,5 /path/to/ligand_roles/1.role_table.csv | sort | uniq
```

## Ligand site matching in mechanism.yaml

Precise selection:

```yaml
ligand_sites:
  acceptor_o:
    atom_labels:
      - hydroxyl.o
    required: true
```

Class selection:

```yaml
ligand_sites:
  carbonyl_c:
    atom_classes:
      - acyl_carbonyl_c
      - ester_carbonyl_c
    required: true
```

Linked atom selection:

```yaml
ligand_sites:
  carbonyl_c:
    atom_classes:
      - acyl_carbonyl_c
    required: true

  carbonyl_o:
    linked_to: carbonyl_c
    atom_classes:
      - carbonyl_oxygen
    required: true
```

`linked_to` requires the same `group_id` and `instance_id` as the target ligand site.

## PDBQT order rule

`pdbqt_order` is the atom line order inside each PDBQT pose. It is not necessarily the ATOM serial number printed in the PDBQT file.
