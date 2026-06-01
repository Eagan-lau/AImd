# Module Architecture

## Core modules

```text
LigandRoleMap
  Original ligand chemistry -> stable ligand atom labels.

LigandPoseParser
  Multi-pose PDBQT -> pose-specific coordinates by pdbqt_order.

MechanismConfig
  profile.yaml + mechanism.yaml + CLI paths.

ProteinSiteMap
  explicit anchor -> radius search -> protein role atoms.

TemplateAligner
  structure-first template alignment and local target atom search.

GeometryEngine
  distance, angle_3pt, best_fit_plane, oriented_axis, axis_deviation.

CatalyticGate
  required feature gate filtering.

CatalyticScorer
  affinity and hierarchical geometry scoring.

PyMOLExporter
  optional selected candidate visualization.

PaperLockedBackend
  exact original-script reproduction.
```

## Data flow

```text
role_table.csv
  + docked ligand PDBQT
  + protein PDBQT
  + mechanism.yaml
  + profile.yaml
      -> resolved_ligand_sites.csv
      -> resolved_protein_roles.csv
      -> geometry_features.csv
      -> passing_candidates.csv
      -> candidate_scores.csv
      -> pose_scores.csv
      -> conformation_scores.csv
      -> protein_scores.csv
```

## PyMOL-free core

The core calculation does not use PyMOL.

PyMOL is only used for:

```text
visual inspection
PML generation
PSE generation
figure preparation
```

## No global pocket

ProteinSiteMap uses explicit anchors only.

Allowed anchors:

```text
ligand.<site>
protein.<role>
template.<template_name>
```
