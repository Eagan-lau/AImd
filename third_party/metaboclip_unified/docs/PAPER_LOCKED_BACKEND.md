# Paper-Locked Backend

The paper-locked backend is for exact manuscript-result reproduction.

It does not use:

```text
LigandRoleMap
ProteinSiteMap
CatalyticGate
CatalyticScorer
mechanism.yaml
```

It directly runs external original scripts after checksum verification.

## Prepare original scripts

```bash
metaboclip prepare-paper-locked \
  --archive /path/to/original_scoring_modules.zip \
  --dest /path/to/paper_locked_originals \
  --overwrite \
  --report paper_locked_prepare_report.json
```

## Verify original scripts

```bash
metaboclip verify-paper-locked \
  --original-root /path/to/paper_locked_originals \
  --report verify_all.json
```

## Run original scripts

```bash
metaboclip run-paper-locked \
  --family act \
  --original-root /path/to/paper_locked_originals \
  --stage both \
  --file-range 1-100 \
  --out-dir /path/to/paper_locked_logs/act \
  --report paper_locked_act_run.json
```

## Purpose

Use paper-locked mode when you need byte-for-byte legacy behavior.

Use generic or curated YAML mode when you need the new MetaboClip modular framework.
