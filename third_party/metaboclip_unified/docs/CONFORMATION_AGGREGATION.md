# Conformation and Protein-Level Aggregation

## Important distinction

One protein PDBQT file is one conformation, not one protein.

Example:

```text
Tca01g00352.t1_pLDDT=91.1_pTM=0.91.pdbqt    -> conformation_id 0
Tca01g00352.t1_pLDDT=91.1_pTM=0.91_1.pdbqt  -> conformation_id 1
Tca01g00352.t1_pLDDT=91.1_pTM=0.91_2.pdbqt  -> conformation_id 2
Tca01g00352.t1_pLDDT=91.1_pTM=0.91_3.pdbqt  -> conformation_id 3
Tca01g00352.t1_pLDDT=91.1_pTM=0.91_4.pdbqt  -> conformation_id 4
Tca01g00352.t1_pLDDT=91.1_pTM=0.91_5.pdbqt  -> conformation_id 5
```

The first file without suffix is treated as conformation 0.

## Output levels

Each conformation folder contains:

```text
candidate_scores.csv
pose_scores.csv
conformation_scores.csv
resolved_ligand_sites.csv
resolved_protein_roles.csv
geometry_features.csv
passing_candidates.csv
```

The top-level output directory contains protein-level files:

```text
merged_conformation_scores.csv
protein_scores.csv
summary.json
```

## ProteinScore

ProteinScore is calculated only after all conformations are collected for the same base protein ID.

Default scoring profile:

```text
active_conformation = ConformationScore >= coverage_threshold
coverage = active_conformation_count / total_conformations
quality_score = max ConformationScore
coverage_factor = (1 - alpha) + alpha * min(coverage / coverage_target, 1.0)
ProteinScore = quality_score * coverage_factor
```

Default parameters:

```yaml
coverage_threshold: 70.0
coverage_target: 0.70
alpha: 0.30
total_conformations: 6
```

## Common mistake

Do not interpret `protein_scores.csv` inside a conformation folder as a valid protein-level score. Protein-level scores must be generated after grouping all conformations.
