import json

def summary_placeholder_runner(cfg, protein_dir, docking_dir, out_dir, registry):
    protein_ext = cfg.get("input", {}).get("protein_extension", ".pdbqt")
    proteins = sorted(protein_dir.glob(f"*{protein_ext}"))
    ligand_id = docking_dir.name.split("_")[-1] if docking_dir.name.startswith("file_") else None
    summary = {
        "family": cfg["family"]["name"],
        "runner": cfg["family"]["runner"],
        "protein_count": len(proteins),
        "ligand_id": ligand_id,
        "status": "placeholder_runner",
    }
    with open(out_dir / "execution_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
    print(f"Execution package loaded for family: {cfg['family']['name']}")
    print(f"Protein files detected: {len(proteins)}")
    print(f"Summary written to: {out_dir / 'execution_summary.json'}")
    return summary
