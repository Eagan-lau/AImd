from pathlib import Path

from metaboclip.paper_locked.runner import load_manifest, normalize_family, PaperLockedRunner


def test_paper_locked_manifest_loads():
    manifest = load_manifest()
    assert "families" in manifest
    assert set(manifest["families"]) == {"ugt", "act", "cyp450", "fe2og", "ach"}
    assert manifest["families"]["ugt"]["gating_script"] == "6_20250810_score.py"


def test_family_aliases():
    assert normalize_family("P450") == "cyp450"
    assert normalize_family("2-ODD") == "fe2og"
    assert normalize_family("dAC") == "ach"


def test_dry_command_building(tmp_path: Path):
    root = tmp_path / "originals"
    fam = root / "ugt"
    fam.mkdir(parents=True)
    (fam / "6_20250810_score.py").write_text("print('gate')\n", encoding="utf-8")
    (fam / "7_score_UGT-20250901.py").write_text("print('score')\n", encoding="utf-8")
    runner = PaperLockedRunner(root, "ugt", skip_checksums=True)
    cmd = runner.build_command("gate", file_range="1-3")
    assert cmd[-2:] == ["--file-range", "1-3"]
    assert cmd[-3].endswith("6_20250810_score.py")
