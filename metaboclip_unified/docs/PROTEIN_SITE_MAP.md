# ProteinSiteMap

ProteinSiteMap resolves protein-side catalytic atoms and residues from explicit anchors.

## Core rule

There is no global pocket in the default ProteinSiteMap logic.

Every protein role is resolved from one explicit anchor:

```text
ligand.<site>
protein.<previous_role>
template.<template_name>
```

Each role search does:

```text
anchor coordinate
  -> radius search
  -> residue / atom / element filtering
  -> candidate role atoms
```

## Ligand-seeded role

```yaml
protein_roles:
  - role: nucleophile
    from: ligand.carbonyl_c
    radius: 6.0
    required: true
    residues:
      SER: [OG]
      CYS: [SG]
      ASP: [OD1, OD2]
      GLU: [OE1, OE2]
```

This means:

```text
For each ligand carbonyl_c coordinate, find allowed protein atoms within 6 A.
```

## Protein-seeded support role

```yaml
protein_roles:
  - role: base
    from: protein.nucleophile
    radius: 5.0
    required: true
    residues:
      HIS: [ND1, NE2]
```

This means:

```text
For each resolved nucleophile, find HIS ND1 or HIS NE2 within 5 A.
```

## Template-seeded role

```yaml
protein_templates:
  ach_template:
    path: resources/ach_template.pdb
    align_to_protein: true
    align_method: structure_global
    align_backend: tmalign
    executable: TMalign
    required: true

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
      require_residue_type_match: false
    validate:
      from: ligand.carbonyl_c
      radius: 6.0
```

This means:

```text
Align template to target.
Project template SER171 OG into target coordinates.
Search nearby target SER OG or CYS SG.
Validate the mapped atom against ligand.carbonyl_c.
```

## Branching behavior

If a role has multiple candidates, every candidate is retained as a branch.

Example:

```text
ligand.carbonyl_c
  -> nucleophile candidate 1
      -> base candidate 1
      -> base candidate 2
  -> nucleophile candidate 2
      -> base candidate 3
```

Every complete branch becomes one `site_set_id`.

## Collection roles

Some roles are atom collections, such as heme nitrogen atoms.

```yaml
protein_roles:
  - role: heme_n_atoms
    from: protein.catalytic_fe
    radius: 3.0
    required: true
    collection: true
    min_count: 3
    preferred_count: 4
    same_residue: true
    residues:
      HEM: [NA, NB, NC, ND]
      HEME: [NA, NB, NC, ND]
      HEC: [NA, NB, NC, ND]
```

Use `collection: true` for groups used in plane or axis construction.

## Metal or cofactor atoms

```yaml
protein_roles:
  - role: catalytic_fe
    from: ligand.reactive_c
    radius: 6.0
    required: true
    atoms:
      - elem Fe
      - name FE
```

## Output

ProteinSiteMap writes:

```text
resolved_protein_roles.csv
```

Important columns:

```text
protein_id
conformation_id
ligand_id
pose_id
site_set_id
role
from_anchor
chain
resi
resn
atom_name
element
x
y
z
distance_to_anchor
```
