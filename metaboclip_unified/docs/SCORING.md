# Catalytic Gate and Scoring

## Core rule

MetaboClip first applies hard geometry filters, then scores only passing candidates.

```text
CatalyticGate
  -> pass / fail

CatalyticScorer
  -> ranked scores for passing candidates
```

## CandidateScore

Default formula:

```text
CandidateScore = 100 * (0.30 * AffinityScore + 0.70 * HierarchicalGeometryScore)
```

This keeps the original family-script idea that affinity is important but geometry is dominant.

## AffinityScore

Default transform:

```text
affinity <= -7.0  -> 1.0
affinity >= -3.0  -> 0.0
otherwise          -> linear interpolation
```

Profile settings:

```yaml
affinity:
  scoring:
    full: -7.0
    zero: -3.0
    weight: 0.30
```

## HierarchicalGeometryScore

Geometry is organized by ProteinSiteMap levels.

```text
Level 1:
  ligand atom label -> first protein role

Level 2:
  first protein role -> second protein role

Level 3:
  second protein role -> third protein role
```

Default level weights:

```yaml
level_weights:
  1: 0.60
  2: 0.25
  3: 0.15
```

If only Level 1 is declared, Level 1 is normalized to 1.0.
If Level 1 and Level 2 are declared, weights are normalized over those two levels.
Missing undeclared levels do not penalize the score.

## Gate parameters

Example:

```yaml
features:
  attack_distance:
    type: distance
    atoms:
      - protein.nucleophile
      - ligand.carbonyl_c
    gate:
      min: 2.0
      max: 5.5
      required: true
```

Rules:

```text
required feature missing -> candidate fails
required feature out of range -> candidate fails
optional feature missing -> ignored
optional feature out of range -> ignored unless configured otherwise
```

## Feature transforms

Distance:

```yaml
score:
  level: 1
  transform: distance_piecewise
  best: 3.2
  cutoff: 5.5
  weight: 0.65
```

Angle:

```yaml
score:
  level: 1
  transform: angle_gaussian
  target: 107.0
  sigma: 20.0
  weight: 0.35
```

Axis deviation:

```yaml
score:
  level: 1
  transform: angle_deviation
  best: 0.0
  cutoff: 45.0
  weight: 0.30
```

## Aggregation

```text
CandidateScore -> PoseScore -> ConformationScore -> ProteinScore
```

Default:

```text
PoseScore = max CandidateScore within a pose
ConformationScore = max PoseScore within a protein conformation
ProteinScore = coverage-weighted aggregation over conformations
```
