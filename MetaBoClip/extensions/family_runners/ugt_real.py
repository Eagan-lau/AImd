from pathlib import Path
import pandas as pd

from core.pymol_utils import (
    get_pymol,
    parse_vina_out,
    split_ligand_states,
    detect_clashes_pairs,
    dist,
    calculate_angle,
    save_filtered_pse,
)
from core.scoring import (
    energy_linear,
    distance_piecewise,
    gaussian_score,
    coupled_mean_positive,
    aggregate_standard,
)
from core.output import write_rows_csv, write_score_tables

def ugt_real_runner(cfg, protein_dir, docking_dir, out_dir, registry):
    protein_ext = cfg.get("input", {}).get("protein_extension", ".pdbqt")
    proteins = sorted(protein_dir.glob(f"*{protein_ext}"))
    ligand_id = docking_dir.name.split("_")[-1]
    rows = []

    cmd, pm = get_pymol()
    try:
        for protein_path in proteins:
            rows.extend(_process_single_protein(cmd, cfg, protein_path, docking_dir, ligand_id, out_dir))
            try:
                cmd.delete("all")
            except Exception:
                pass
    finally:
        try:
            if pm is not None:
                pm.stop()
        except Exception:
            pass

    gating_csv = out_dir / f"file_{ligand_id}.gating.csv"
    write_rows_csv(rows, gating_csv)
    score_tables = _score_rows(rows, cfg)
    write_score_tables(score_tables, out_dir, f"file_{ligand_id}")
    print(f"UGT real runner finished. Gating rows: {len(rows)}")
    print(f"Gating CSV: {gating_csv}")
    return {"rows": len(rows)}

def _process_single_protein(cmd, cfg, protein_path, docking_dir, ligand_id, out_dir):
    rows = []
    protein_name = protein_path.stem
    input_cfg = cfg["input"]
    ligand_ext = input_cfg.get("ligand_extension", ".pdbqt")
    ligand_pattern = input_cfg.get("ligand_pattern", "{ligand_id}@{protein_name}{ligand_extension}")
    out_pattern = input_cfg.get("out_pattern", "{ligand_id}_{protein_name}_cavity_1.out")
    ligand_file = docking_dir / ligand_pattern.format(ligand_id=ligand_id, protein_name=protein_name, ligand_extension=ligand_ext)
    out_file = docking_dir / out_pattern.format(ligand_id=ligand_id, protein_name=protein_name)

    if not protein_path.exists() or not ligand_file.exists() or not out_file.exists():
        return rows

    affinity_threshold = float(input_cfg.get("affinity_threshold", -3.0))
    list1, affinity_dict = parse_vina_out(str(out_file), affinity_threshold)
    if not list1:
        return rows

    cmd.reinitialize()
    cmd.load(str(protein_path), "protein1")

    donor_c1_sel = cfg["cofactor"]["atoms"]["donor_c1"]["selector"].format(protein_obj="protein1")
    donor_o1_sel = cfg["cofactor"]["atoms"]["donor_o1"]["selector"].format(protein_obj="protein1")

    cmd.select("C1", donor_c1_sel)
    if cmd.count_atoms("C1") == 0:
        return rows
    zb1 = cmd.get_coords("C1", 1)[0]
    stored = {"c1": {}}
    cmd.iterate("C1", "c1['id'] = ID", space=stored)

    cmd.select("near_His", "protein1 and resn HIS within 10 of C1")
    if cmd.count_atoms("near_His") == 0:
        return rows

    cmd.select("NE2_atoms", "near_His and name NE2")
    zb2_list = cmd.get_coords("NE2_atoms", 1)
    if zb2_list is None or len(zb2_list) == 0:
        return rows

    stored["ne2_list"] = []
    cmd.iterate("NE2_atoms", "ne2_list.append({'id': ID})", space=stored)

    cmd.load(str(ligand_file), "ligand")
    all_models = split_ligand_states(cmd, "ligand", "ligand_")
    retained_models = [f"ligand_{str(mode).zfill(4)}" for mode in list1 if f"ligand_{str(mode).zfill(4)}" in cmd.get_names("objects")]

    for m in all_models:
        if m not in retained_models:
            cmd.delete(m)

    passing_models = set()

    for model in retained_models:
        clash_cutoff = float(input_cfg.get("clash_cutoff", 2.0))
        if detect_clashes_pairs(cmd, "protein1", model, clash_cutoff):
            continue

        cmd.select(f"oh_oxygens_{model}", f"({model} and name O) and (neighbor name H)")
        zb3_list = cmd.get_coords(f"oh_oxygens_{model}", 1)
        if zb3_list is None or len(zb3_list) == 0:
            continue

        stored["o_list"] = []
        cmd.iterate(f"oh_oxygens_{model}", "o_list.append({'id': ID})", space=stored)

        mode_num = int(model.split("_")[1])
        affinity = affinity_dict.get(mode_num, float("nan"))
        results_for_model = []

        for o_idx, zb3 in enumerate(zb3_list):
            for ne2_idx, zb2 in enumerate(zb2_list):
                d_ne2_o = dist(zb2, zb3)
                if not (d_ne2_o < 5.0):
                    continue

                d_c1_o = dist(zb1, zb3)
                if not (d_c1_o < 5.0):
                    continue

                cmd.select("O1", donor_o1_sel)
                if cmd.count_atoms("O1") == 0:
                    continue
                zbO1 = cmd.get_coords("O1", 1)[0]
                stored["o1"] = {}
                cmd.iterate("O1", "o1['id'] = ID", space=stored)

                vec1 = [zbO1[0] - zb1[0], zbO1[1] - zb1[1], zbO1[2] - zb1[2]]
                vec2 = [zb3[0] - zb1[0], zb3[1] - zb1[1], zb3[2] - zb1[2]]
                angle = calculate_angle(vec1, vec2)
                if angle is None or not (100.0 < angle < 180.1):
                    continue

                rows.append({
                    "Ligand_id": str(ligand_id),
                    "protein_id": protein_name.split("_")[0],
                    "conformation": protein_name,
                    "mode": int(mode_num),
                    "affinity_kcal": float(affinity),
                    "dist_NE2_O": float(d_ne2_o),
                    "dist_C1_O": float(d_c1_o),
                    "angle_deg": float(angle),
                    "ne2_atom_id": int(stored["ne2_list"][ne2_idx]["id"]),
                    "o_atom_id": int(stored["o_list"][o_idx]["id"]),
                    "c1_atom_id": int(stored["c1"]["id"]),
                    "o1_atom_id": int(stored["o1"]["id"]),
                    "clash_flag": False,
                })
                results_for_model.append(1)

        if results_for_model:
            passing_models.add(model)

    if passing_models and cfg.get("output", {}).get("save_pse", True):
        sessions_dir = out_dir / "sessions"
        sessions_dir.mkdir(parents=True, exist_ok=True)
        save_filtered_pse(cmd, passing_models, sessions_dir / f"{protein_name}.pse")

    return rows

def _score_rows(rows, cfg):
    df = pd.DataFrame(rows)
    if df.empty:
        return {
            "pose_scores": pd.DataFrame(),
            "conformation_scores": pd.DataFrame(),
            "protein_scores": pd.DataFrame(),
        }

    df["affinity_score"] = df["affinity_kcal"].apply(lambda x: energy_linear(x, -7.0, -3.0))
    df["d1_score"] = df["dist_NE2_O"].apply(lambda x: distance_piecewise(x, 3.2, 5.0))
    df["d2_score"] = df["dist_C1_O"].apply(lambda x: distance_piecewise(x, 3.2, 5.0))
    df["D_score"] = [coupled_mean_positive(a, b) for a, b in zip(df["d1_score"], df["d2_score"])]
    df["angle_score"] = df["angle_deg"].apply(lambda x: gaussian_score(x, 180.0, 30.0, flat=25.0, fold180=True))
    df["s_pose"] = 100.0 * (0.30 * df["affinity_score"] + 0.40 * df["D_score"] + 0.30 * df["angle_score"])

    agg = cfg.get("aggregation", {})
    conf_df, prot_df = aggregate_standard(
        df,
        total_confs=int(agg.get("total_confs", 6)),
        cover_t=float(agg.get("cover_t", 70.0)),
        alpha=float(agg.get("alpha", 0.30)),
    )
    return {
        "pose_scores": df,
        "conformation_scores": conf_df,
        "protein_scores": prot_df,
    }
