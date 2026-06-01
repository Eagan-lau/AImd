from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
import sys
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import yaml


PACKAGE_DIR = Path(__file__).resolve().parent
MANIFEST_PATH = PACKAGE_DIR / "manifest.yaml"
FAMILY_ALIASES = {
    "ugt": "ugt",
    "act": "act",
    "cyp450": "cyp450",
    "p450": "cyp450",
    "fe2og": "fe2og",
    "2-odd": "fe2og",
    "2odd": "fe2og",
    "ach": "ach",
    "dac": "ach",
}


@dataclass(frozen=True)
class CommandResult:
    stage: str
    command: list[str]
    returncode: int
    log_file: str | None


def load_manifest(path: Path | None = None) -> dict[str, Any]:
    manifest_path = path or MANIFEST_PATH
    with open(manifest_path, "r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle)
    if not isinstance(data, dict) or "families" not in data:
        raise ValueError(f"Invalid paper-locked manifest: {manifest_path}")
    return data


def normalize_family(family: str) -> str:
    key = family.strip().lower()
    if key not in FAMILY_ALIASES:
        allowed = ", ".join(sorted(set(FAMILY_ALIASES.values())))
        raise ValueError(f"Unknown paper-locked family: {family}. Allowed: {allowed}")
    return FAMILY_ALIASES[key]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _family_dir_candidates(root: Path, source_dir: str, family_key: str) -> list[Path]:
    names = [source_dir, family_key, family_key.upper()]
    candidates: list[Path] = []
    for name in names:
        p = root / name
        if p.is_dir():
            candidates.append(p)
    for p in root.rglob(source_dir):
        if p.is_dir():
            candidates.append(p)
    seen: set[Path] = set()
    unique: list[Path] = []
    for p in candidates:
        rp = p.resolve()
        if rp not in seen:
            seen.add(rp)
            unique.append(p)
    return unique


def find_family_dir(original_root: Path, family_key: str, manifest: dict[str, Any] | None = None) -> Path:
    data = manifest or load_manifest()
    family = data["families"][normalize_family(family_key)]
    source_dir = family["source_dir"]
    candidates = _family_dir_candidates(original_root, source_dir, normalize_family(family_key))
    if not candidates:
        raise FileNotFoundError(f"Could not find source directory {source_dir} under {original_root}")
    return candidates[0]


def verify_original_scripts(original_root: Path, family: str | None = None, manifest_path: Path | None = None) -> dict[str, Any]:
    data = load_manifest(manifest_path)
    families = [normalize_family(family)] if family else sorted(data["families"].keys())
    report: dict[str, Any] = {"root": str(original_root), "families": {}}
    for fam in families:
        spec = data["families"][fam]
        fam_dir = find_family_dir(original_root, fam, data)
        fam_report: dict[str, Any] = {"family_dir": str(fam_dir), "files": {}, "ok": True}
        for filename, expected in spec.get("expected_sha256", {}).items():
            path = fam_dir / filename
            if not path.exists():
                fam_report["files"][filename] = {"exists": False, "ok": False, "sha256": None, "expected": expected}
                fam_report["ok"] = False
                continue
            got = sha256_file(path)
            ok = got == expected
            fam_report["files"][filename] = {"exists": True, "ok": ok, "sha256": got, "expected": expected}
            fam_report["ok"] = fam_report["ok"] and ok
        report["families"][fam] = fam_report
    report["ok"] = all(item["ok"] for item in report["families"].values())
    return report


def _extract_archive(archive: Path, scratch_dir: Path) -> Path:
    scratch_dir.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(archive, "r") as zf:
        zf.extractall(scratch_dir)
    return scratch_dir


def prepare_original_scripts(archive: Path, dest: Path, overwrite: bool = False, manifest_path: Path | None = None) -> dict[str, Any]:
    data = load_manifest(manifest_path)
    dest.mkdir(parents=True, exist_ok=True)
    scratch = dest / ".extract_tmp"
    if scratch.exists():
        shutil.rmtree(scratch)
    _extract_archive(archive, scratch)
    report: dict[str, Any] = {"archive": str(archive), "dest": str(dest), "families": {}}
    for fam, spec in data["families"].items():
        src_dir = find_family_dir(scratch, fam, data)
        out_dir = dest / fam
        if out_dir.exists():
            if overwrite:
                shutil.rmtree(out_dir)
            else:
                raise FileExistsError(f"Destination family directory already exists: {out_dir}")
        out_dir.mkdir(parents=True)
        copied: list[str] = []
        for filename in spec.get("expected_sha256", {}):
            src = src_dir / filename
            if src.exists():
                shutil.copy2(src, out_dir / filename)
                copied.append(filename)
        report["families"][fam] = {"source_dir": str(src_dir), "prepared_dir": str(out_dir), "files": copied}
    shutil.rmtree(scratch)
    report["verification"] = verify_original_scripts(dest, manifest_path=manifest_path)
    return report


class PaperLockedRunner:
    def __init__(
        self,
        original_root: Path,
        family: str,
        python_executable: str | None = None,
        manifest_path: Path | None = None,
        skip_checksums: bool = False,
    ) -> None:
        self.original_root = original_root
        self.family_key = normalize_family(family)
        self.python_executable = python_executable or sys.executable
        self.manifest = load_manifest(manifest_path)
        self.spec = self.manifest["families"][self.family_key]
        self.family_dir = find_family_dir(original_root, self.family_key, self.manifest)
        if not skip_checksums:
            report = verify_original_scripts(original_root, self.family_key, manifest_path)
            if not report.get("ok"):
                raise RuntimeError("Original script checksum verification failed. Use verify-paper-locked for details.")

    def build_command(
        self,
        stage: str,
        file_range: str | None = None,
        docking_base_dir: Path | None = None,
        extra_args: Iterable[str] | None = None,
    ) -> list[str]:
        if stage not in {"gate", "score"}:
            raise ValueError(f"Invalid stage: {stage}")
        script_key = "gating_script" if stage == "gate" else "scoring_script"
        script = self.family_dir / self.spec[script_key]
        if not script.exists():
            raise FileNotFoundError(script)
        cmd = [self.python_executable, str(script)]
        if stage == "gate":
            if docking_base_dir is not None and self.spec.get("gating_docking_base_arg"):
                cmd.extend([self.spec["gating_docking_base_arg"], str(docking_base_dir)])
            if file_range:
                arg = self.spec.get("gating_range_arg")
                if arg:
                    cmd.extend([arg, file_range])
        else:
            cmd.extend(self.spec.get("scoring_args", []))
        if extra_args:
            cmd.extend(list(extra_args))
        return cmd

    def run_stage(
        self,
        stage: str,
        file_range: str | None = None,
        docking_base_dir: Path | None = None,
        extra_args: Iterable[str] | None = None,
        out_dir: Path | None = None,
        dry_run: bool = False,
    ) -> CommandResult:
        cmd = self.build_command(stage, file_range=file_range, docking_base_dir=docking_base_dir, extra_args=extra_args)
        log_file: Path | None = None
        if out_dir:
            out_dir.mkdir(parents=True, exist_ok=True)
            log_file = out_dir / f"paper_locked_{self.family_key}_{stage}.log"
        if dry_run:
            if log_file:
                log_file.write_text(" ".join(cmd) + "\n", encoding="utf-8")
            return CommandResult(stage=stage, command=cmd, returncode=0, log_file=str(log_file) if log_file else None)
        env = os.environ.copy()
        env.setdefault("PYTHONUNBUFFERED", "1")
        if log_file:
            with open(log_file, "w", encoding="utf-8") as handle:
                proc = subprocess.run(cmd, cwd=str(self.family_dir), env=env, stdout=handle, stderr=subprocess.STDOUT, text=True)
        else:
            proc = subprocess.run(cmd, cwd=str(self.family_dir), env=env, text=True)
        if proc.returncode != 0:
            raise RuntimeError(f"Paper-locked stage failed: family={self.family_key} stage={stage} returncode={proc.returncode}")
        return CommandResult(stage=stage, command=cmd, returncode=proc.returncode, log_file=str(log_file) if log_file else None)

    def run(
        self,
        stage: str = "both",
        file_range: str | None = None,
        docking_base_dir: Path | None = None,
        out_dir: Path | None = None,
        dry_run: bool = False,
        extra_gate_args: Iterable[str] | None = None,
        extra_score_args: Iterable[str] | None = None,
    ) -> dict[str, Any]:
        stages = ["gate", "score"] if stage == "both" else [stage]
        results: list[CommandResult] = []
        for item in stages:
            extra = extra_gate_args if item == "gate" else extra_score_args
            results.append(
                self.run_stage(
                    item,
                    file_range=file_range,
                    docking_base_dir=docking_base_dir,
                    extra_args=extra,
                    out_dir=out_dir,
                    dry_run=dry_run,
                )
            )
        return {
            "family": self.family_key,
            "family_dir": str(self.family_dir),
            "stage": stage,
            "dry_run": dry_run,
            "results": [r.__dict__ for r in results],
        }


def write_json_report(report: dict[str, Any], path: Path | None) -> None:
    if path is None:
        print(json.dumps(report, indent=2))
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2)
