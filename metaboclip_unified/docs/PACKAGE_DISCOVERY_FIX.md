# Editable Install Package Discovery Fix

## Symptom

Editable installation fails with:

```text
Multiple top-level packages discovered in a flat-layout
```

This happens when runtime folders such as `results`, `ligand_roles`, or `rules` are inside the source tree and setuptools tries to treat them as Python packages.

## Fix

Edit `pyproject.toml` and ensure package discovery is explicit:

```toml
[tool.setuptools.packages.find]
where = ["."]
include = ["metaboclip*", "metaboclip_ligand_roles*"]
exclude = ["results*", "ligand_roles*", "rules*", "examples*", "docs*", "tests*", "paper_locked_originals*"]
```

Then run:

```bash
pip install -e .
```

## Better long-term practice

Move runtime output folders outside the source tree:

```bash
mkdir -p /path/to/runtime/results
mkdir -p /path/to/runtime/ligand_roles
```

Use absolute paths in all commands.
