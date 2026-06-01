# AImd third_party

This directory centralizes external software, executable paths, models, templates, and databases used by AImd modules.

## Check tools

```bash
cd AImd
python third_party/check_tools.py --config third_party/tools.yaml --root .
```

## Use a local symlink

```bash
python third_party/setup_links.py hipmcl /absolute/path/to/hipmcl --root .
```

Then set in `third_party/tools.yaml`:

```yaml
tools:
  hipmcl:
    executable: third_party/bin/hipmcl
```

AImd modules should resolve tools from this registry before falling back to system PATH.
