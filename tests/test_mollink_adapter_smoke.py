from __future__ import annotations

import pandas as pd

from MolLink.runner import run_mollink


def test_mollink_csv_only_smoke(tmp_path):
    molecule_table = tmp_path / "taxane_molecules.csv"
    molecule_table.write_text(
        "source_id,smiles,molecule_name\n"
        "taxane_1,CCO,Example alcohol\n"
        "taxane_2,CC=O,Example aldehyde\n",
        encoding="utf-8",
    )
    output_root = tmp_path / "ligand_transformation"
    result = run_mollink({
        "project": {"root": str(tmp_path)},
        "inputs": {"molecule_table": str(molecule_table), "template_file": None},
        "columns": {"smiles": "smiles", "id": "source_id", "name": "molecule_name"},
        "compute": {"mode": "csv_only", "prefix": "transformnet"},
        "network": {"prefix": "transformnet_network", "write_excel": False, "write_graphs": False, "write_html": True},
        "outputs": {"output_root": str(output_root)},
    })

    manifest_path = result.paths["ligand_source_manifest"]
    manifest = pd.read_csv(manifest_path)
    assert result.summary["mode"] == "csv_only_no_template"
    assert result.summary["valid_molecules"] == 2
    assert result.summary["directed_edges"] == 0
    assert manifest["ligand_id"].tolist() == ["taxane_1", "taxane_2"]
    assert (output_root / "network" / "transformnet_network.nodes.csv").is_file()
