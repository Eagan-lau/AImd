# Hardening Notes

This version keeps the same MetaboClip logic and strengthens runtime reliability.

## Core hardening changes

1. Ligand PDBQT pose atom order is now based on atom line order inside each pose.
   The original PDBQT file serial is stored separately as `file_serial`.

2. Linked ligand sites now require the same `group_id` and `instance_id`.
   This prevents mixing atoms from different functional group instances.

3. Optional missing geometry features are ignored during scoring.
   Missing optional support roles no longer add zero-valued score levels.

4. Protein collection roles support:
   - `min_count`
   - `preferred_count`
   - `same_residue`

5. The mechanism validator is available through:

```bash
metaboclip validate-mechanism --mechanism path/to/mechanism.yaml
```

6. `site_set_id` values now include pose IDs.

## Recommended mechanism checks

Before running a new mechanism file:

```bash
metaboclip validate-mechanism --mechanism examples/new_family/mechanism.yaml
```

## PDBQT note

For ligand docking pose PDBQT files, `pdbqt_order` means atom line order inside each pose, not the `ATOM` serial value printed in the file.

## v2.2 conformation-level score aggregation

Protein score is now computed only after all conformations belonging to the same base protein have been processed.

A protein conformation name ending in an underscore followed by digits is interpreted as an alternative conformation. The unsuffixed structure is interpreted as conformation 0.

Examples:

```text
Tca01g00352.t1_pLDDT_91.1_pTM_0.91   -> protein_id=Tca01g00352.t1_pLDDT_91.1_pTM_0.91, conformation_id=0
Tca01g00352.t1_pLDDT_91.1_pTM_0.91_1 -> protein_id=Tca01g00352.t1_pLDDT_91.1_pTM_0.91, conformation_id=1
Tca01g00352.t1_pLDDT_91.1_pTM_0.91_5 -> protein_id=Tca01g00352.t1_pLDDT_91.1_pTM_0.91, conformation_id=5
```

Per-conformation output folders contain candidate, pose, and conformation scores only. A merged conformation table and the final protein score table are written at the run output root:

```text
merged_conformation_scores.csv
protein_scores.csv
```

This prevents a single PDBQT conformation from being scored as a full protein.
