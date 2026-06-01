from __future__ import annotations

import shutil
from pathlib import Path


def ensure_dir(path: str | Path) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def clean_dir(path: str | Path) -> Path:
    p = Path(path)
    if p.exists():
        shutil.rmtree(p)
    p.mkdir(parents=True, exist_ok=True)
    return p


def rel_or_abs(path: str | Path, root: str | Path) -> str:
    p = Path(path).resolve()
    r = Path(root).resolve()
    try:
        return str(p.relative_to(r))
    except ValueError:
        return str(p)


def resolve_path(path: str | Path, root: str | Path) -> Path:
    p = Path(path)
    if p.is_absolute():
        return p
    return (Path(root) / p).resolve()


def symlink_or_copy(src: str | Path, dst: str | Path, overwrite: bool = True) -> None:
    src = Path(src).resolve()
    dst = Path(dst).resolve()
    if overwrite and (dst.exists() or dst.is_symlink()):
        dst.unlink()
    dst.parent.mkdir(parents=True, exist_ok=True)
    try:
        dst.symlink_to(src)
    except OSError:
        shutil.copy2(src, dst)
