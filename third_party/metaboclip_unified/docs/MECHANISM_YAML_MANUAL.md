# Mechanism YAML Manual

`mechanism.yaml` is the only file users usually edit when adding a new enzyme family.

Keep it focused on computation. Do not put family metadata, output paths, or global score defaults here.

## Minimal sections

```text
ligand_sites
protein_templates
protein_roles
geometry_refs
features
```

`protein_templates` and `geometry_refs` are optional.

## Minimal template

```yaml
ligand_sites:
  target:
    atom_labels: []
    atom_classes:
      - acyl_carbonyl_c
    required: true

  angle_atom:
    linked_to: target
    atom_labels: []
    atom_classes:
      - carbonyl_oxygen
    required: false

protein_templates: {}

protein_roles:
  - role: atom1
    from: ligand.target
    radius: 6.0
    required: true
    residues:
      SER: [OG]
      CYS: [SG]
      ASP: [OD1, OD2]
      GLU: [OE1, OE2]

  - role: atom2
    from: protein.atom1
    radius: 5.0
    required: false
    residues:
      HIS: [ND1, NE2]

features:
  primary_distance:
    type: distance
    atoms:
      - protein.atom1
      - ligand.target
    gate:
      min: 2.0
      max: 6.0
      required: true
    score:
      level: 1
      transform: distance_piecewise
      best: 3.2
      cutoff: 6.0
      weight: 1.0
```

## ligand_sites

Use `atom_labels` for exact labels:

```yaml
ligand_sites:
  acceptor_o:
    atom_labels:
      - hydroxyl.o
    required: true
```

Use `atom_classes` for broader selection:

```yaml
ligand_sites:
  carbonyl_c:
    atom_classes:
      - acyl_carbonyl_c
      - ester_carbonyl_c
    required: true
```

Use `linked_to` for atoms from the same functional group instance:

```yaml
ligand_sites:
  carbonyl_o:
    linked_to: carbonyl_c
    atom_classes:
      - carbonyl_oxygen
    required: true
```

## protein_roles

Role anchors can be:

```text
ligand.<site>
protein.<role>
template.<template_name>
```

Ligand-seeded role:

```yaml
protein_roles:
  - role: catalytic_base
    from: ligand.acceptor_o
    radius: 6.0
    required: true
    residues:
      HIS: [ND1, NE2]
```

Protein-seeded role:

```yaml
protein_roles:
  - role: acid
    from: protein.base
    radius: 5.5
    required: false
    residues:
      ASP: [OD1, OD2]
      GLU: [OE1, OE2]
```

Template-seeded role:

```yaml
protein_roles:
  - role: nucleophile
    from: template.ach_template
    required: true
    template_atom:
      residue_name: SER
      residue_id: 171
      atom_name: OG
    target_search:
      residues:
        SER: [OG]
        CYS: [SG]
      search_radius: 4.0
      choose: nearest
    validate:
      from: ligand.carbonyl_c
      radius: 6.0
```

## geometry_refs

Use this only when a feature needs a plane, axis, or normal vector.

P450-like example:

```yaml
geometry_refs:
  heme_plane:
    type: best_fit_plane
    atoms:
      - protein.heme_n_atoms

  distal_axis:
    type: oriented_axis
    origin: protein.catalytic_fe
    plane: geometry.heme_plane
    orient_away_from: protein.proximal_cys
```

## features

Distance feature:

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
    score:
      level: 1
      transform: distance_piecewise
      best: 3.2
      cutoff: 5.5
      weight: 0.65
```

Three-point angle feature:

```yaml
features:
  attack_angle:
    type: angle_3pt
    enabled: true
    a: protein.nucleophile
    vertex: ligand.carbonyl_c
    c: ligand.carbonyl_o
    gate:
      min: 90.0
      max: 140.0
      required: true
    score:
      level: 1
      transform: angle_gaussian
      target: 107.0
      sigma: 20.0
      weight: 0.35
```

Axis deviation feature:

```yaml
features:
  heme_axis_deviation:
    type: axis_deviation
    axis: geometry.distal_axis
    vector:
      from: protein.catalytic_fe
      to: ligand.reactive_c
    gate:
      min: 0.0
      max: 45.0
      required: true
    score:
      level: 1
      transform: angle_deviation
      best: 0.0
      cutoff: 45.0
      weight: 0.30
```

## Required fields to edit for a new family

```text
ligand_sites
protein_roles
features
```

Edit `protein_templates` only if template mapping is needed.
Edit `geometry_refs` only if plane or axis features are needed.
