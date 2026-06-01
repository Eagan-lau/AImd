# Installation

## Requirements

Recommended environment:

```text
Python >= 3.10
numpy
scipy
PyYAML
rdkit
```

Optional tools:

```text
PyMOL    Only for visualization and PSE export.
TMalign  Optional structure-first template alignment backend.
US-align Optional alternative template alignment backend.
```

## Install in editable mode

```bash
cd /path/to/metaboclip-unified-v2.1-hardened
pip install -e .
```

## Verify commands

```bash
metaboclip --help
metaboclip validate-mechanism --help
metaboclip export-pymol --help
metaboclip export-pymol-batch --help
```

## If editable install fails

If you see:

```text
Multiple top-level packages discovered in a flat-layout
```

see:

```text
docs/PACKAGE_DISCOVERY_FIX.md
```

The usual fix is to make package discovery explicit in `pyproject.toml` and keep runtime output folders outside the source tree.

## Runtime output directory recommendation

Do not place large runtime folders inside the source tree when possible:

```text
results/
ligand_roles/
paper_locked_originals/
```

Recommended layout:

```text
/project/metaboclip_source/
/runtime/metaboclip_results/
/runtime/ligand_roles/
/runtime/docking_results/
/runtime/protein_pdbqt/
```
