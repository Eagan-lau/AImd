# Curated Family YAML Files

MetaboClip includes curated mechanism YAML files for five families:

```text
ugt
act
cyp450
fe2og
ach
```

Paths:

```text
metaboclip/config/families/ugt/mechanism.yaml
metaboclip/config/families/act/mechanism.yaml
metaboclip/config/families/cyp450/mechanism.yaml
metaboclip/config/families/fe2og/mechanism.yaml
metaboclip/config/families/ach/mechanism.yaml
```

## UGT

Typical ligand site:

```text
hydroxyl.o
```

Typical protein roles:

```text
catalytic_base
sugar_donor_c1
sugar_donor_o1
```

Typical features:

```text
donor C1 -> hydroxyl O distance
His -> hydroxyl O distance
O1-C1-O angle
```

## ACT

Typical ligand sites:

```text
hydroxyl.o
protic_nucleophile.atom
```

Typical protein roles:

```text
catalytic_base
acyl_donor_c
acyl_donor_o
```

If strict ACT mode produces no candidates, check whether donor/cofactor atoms exist in the protein PDBQT.

## CYP450

Typical ligand site:

```text
c_h_site.c
oxidizable_carbon
```

Typical protein roles:

```text
catalytic_fe
heme_n_atoms
proximal_cys
```

Typical features:

```text
Fe -> reactive C distance
heme axis deviation
```

## Fe2OG

Typical ligand site:

```text
c_h_site.c
oxidizable_carbon
```

Typical protein roles:

```text
catalytic_fe
cofactor atoms or metal context atoms
```

Typical features:

```text
Fe -> reactive C distance
optional cofactor axis deviation
```

## ACH

Typical ligand sites:

```text
acyl_carbonyl_c
ester_carbonyl_c
carbonyl_oxygen
```

Typical protein roles:

```text
nucleophile
base
acid
```

ACH may use template-seeded mode for the nucleophile role.
