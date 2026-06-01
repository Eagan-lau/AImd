# Template Alignment

## When to use templates

Use templates when protein sequence identity is low but the 3D fold and catalytic site are conserved.

In this mode, the template residue number is used only to locate an atom in the template structure. The target protein residue number does not need to match.

## Recommended logic

```text
Template structure
  -> structure-only alignment to target protein
  -> transform template catalytic atom into target coordinates
  -> search nearby target atoms allowed by mechanism.yaml
  -> validate mapped atom against ligand atom label or previous protein role
```

## TMalign backend

Example:

```yaml
protein_templates:
  ach_template:
    path: resources/ach_template.pdb
    align_to_protein: true
    align_method: structure_global
    align_backend: tmalign
    executable: TMalign
    required: true
    quality:
      min_tm_score: 0.45
      max_rmsd: 5.0
      min_aligned_length: 80
```

## Template-seeded role

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
      require_residue_type_match: false
    validate:
      from: ligand.carbonyl_c
      radius: 6.0
```

This means:

```text
1. Align the template structure to the target structure.
2. Transform template SER171 OG into target coordinates.
3. Search target SER OG or CYS SG near that transformed coordinate.
4. Validate that the mapped atom is close to ligand.carbonyl_c.
```

## Templates can be used at any role level

First level:

```text
template -> protein.nucleophile
```

Second level:

```text
ligand -> protein.nucleophile
template -> protein.base
validate base against nucleophile
```

Third level:

```text
ligand -> protein.nucleophile
protein.nucleophile -> protein.base
template -> protein.acid
validate acid against base
```

## PyMOL is not required

Template alignment does not need PyMOL. Core alignment should use TMalign, US-align, internal Kabsch alignment, or a manual mapping table.

PyMOL remains optional for visualization.
