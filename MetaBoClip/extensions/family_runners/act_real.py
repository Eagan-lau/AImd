from pathlib import Path
import pandas as pd

from core.pymol_utils import (
    get_pymol,
    parse_vina_out,
    split_ligand_states,
    has_clash_around,
    dist,
    angle_three_points,
    save_filtered_pse,
)
from core.scoring import energy_linear, distance_piecewise, aggregate_standard
from core.output import write_rows_csv, write_score_tables


def act_real_runner(cfg, protein_dir, docking_dir, out_dir, registry):
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
    print(f"ACT real runner finished. Gating rows: {len(rows)}")
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
    coa_path = Path(cfg["cofactor"]["reference_file"])

    if not protein_path.exists() or not ligand_file.exists() or not out_file.exists() or not coa_path.exists():
        return rows

    affinity_threshold = float(input_cfg.get("affinity_threshold", -3.0))
    modes, mode2aff = parse_vina_out(str(out_file), affinity_threshold)
    if not modes:
        return rows

    cmd.reinitialize()
    cmd.load(str(protein_path), "protein1")
    cmd.load(str(ligand_file), "ligand")
    cmd.load(str(coa_path), "coa")

    donor_c_sel_list = cfg["cofactor"]["atoms"]["donor_c"]["selector_any"]
    donor_o_sel_list = cfg["cofactor"]["atoms"]["donor_o"]["selector_any"]
    donor_s_sel_list = cfg["cofactor"]["atoms"]["donor_s"]["selector_any"]

    def _first_nonempty(selectors):
        for s in selectors:
            sel = s.format(protein_obj="protein1")
            cmd.select("tmp_sel", sel)
            if cmd.count_atoms("tmp_sel") > 0:
                return sel
        return None

    donor_c_sel = _first_nonempty(donor_c_sel_list)
    donor_o_sel = _first_nonempty(donor_o_sel_list)
    donor_s_sel = _first_nonempty(donor_s_sel_list)
    if donor_c_sel is None or donor_o_sel is None or donor_s_sel is None:
        return rows

    cmd.select("DONOR_C", donor_c_sel)
    cmd.select("DONOR_O", donor_o_sel)
    cmd.select("DONOR_S", donor_s_sel)
    if cmd.count_atoms("DONOR_C") == 0 or cmd.count_atoms("DONOR_O") == 0 or cmd.count_atoms("DONOR_S") == 0:
        return rows

    donor_c = cmd.get_coords("DONOR_C", 1)[0]
    donor_o = cmd.get_coords("DONOR_O", 1)[0]
    donor_s = cmd.get_coords("DONOR_S", 1)[0]

    cmd.select("near_His", "protein1 and resn HIS within 10 of DONOR_C")
    if cmd.count_atoms("near_His") == 0:
        return rows
    cmd.select("NE2_atoms", "near_His and name NE2")
    ne2_list = cmd.get_coords("NE2_atoms", 1)
    if ne2_list is None or len(ne2_list) == 0:
        return rows
    stored = {"ne2_list": []}
    cmd.iterate("NE2_atoms", "ne2_list.append({'id': ID})", space=stored)

    all_models = split_ligand_states(cmd, "ligand", "ligand_")
    retained_models = [f"ligand_{str(mode).zfill(4)}" for mode in modes if f"ligand_{str(mode).zfill(4)}" in cmd.get_names("objects")]
    for m in all_models:
        if m not in retained_models:
            cmd.delete(m)

    passing_models = set()
    for model in retained_models:
        clash_cutoff = float(input_cfg.get("clash_cutoff", 2.0))
        if has_clash_around(cmd, model, clash_cutoff, protein_obj="protein1"):
            continue

        nu_atoms = []
        for elem in ["O", "N", "S"]:
            sel = f"({model} and elem {elem}) and (neighbor ({model} and elem H))"
            if cmd.count_atoms(sel) == 0:
                continue
            mm = cmd.get_model(sel, 1)
            for a in mm.atom:
                nu_atoms.append({"id": int(a.id), "coord": [a.coord[0], a.coord[1], a.coord[2]], "elem": elem})
        if not nu_atoms:
            continue

        mode_num = int(model.split("_")[1])
        affinity = mode2aff.get(mode_num, float("nan"))
        any_pass = False
        for nu in nu_atoms:
            for ne2_idx, ne2 in enumerate(ne2_list):
                d_ne2_nu = dist(ne2, nu["coord"])
                d_nu_c = dist(nu["coord"], donor_c)
                d_nu_s = dist(nu["coord"], donor_s)
                if d_ne2_nu >= 5.0 or d_nu_c >= 5.0 or d_nu_s >= 6.0:
                    continue
                bd = angle_three_points(nu["coord"], donor_c, donor_o)
                if not (85.0 <= bd <= 125.0):
                    continue
                rows.append({
                    "Ligand_id": str(ligand_id),
                    "protein_id": protein_name.split("_")[0],
                    "conformation": protein_name,
                    "mode": int(mode_num),
                    "affinity_kcal": float(affinity),
                    "nu_type": nu["elem"],
                    "nu_atom_id": int(nu["id"]),
                    "dist_NE2_Nu": float(d_ne2_nu),
                    "dist_Nu_C": float(d_nu_c),
                    "dist_Nu_S": float(d_nu_s),
                    "BD_angle_deg": float(bd),
                    "ne2_atom_id": int(stored["ne2_list"][ne2_idx]["id"]),
                    "clash_flag": False,
                })
                any_pass = True
        if any_pass:
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

    def act_angle_score(theta):
        if pd.isna(theta):
            return 0.0
        x = float(theta)
        if 100.0 <= x <= 110.0:
            return 1.0
        if 85.0 <= x < 100.0:
            return (x - 85.0) / 15.0
        if 110.0 < x <= 125.0:
            return (125.0 - x) / 15.0
        return 0.0

    df["affinity_score"] = df["affinity_kcal"].apply(lambda x: energy_linear(x, -7.0, -3.0))
    df["d1_score"] = df["dist_NE2_Nu"].apply(lambda x: distance_piecewise(x, 3.2, 5.0))
    df["d2_score"] = df["dist_Nu_C"].apply(lambda x: distance_piecewise(x, 3.2, 5.0))
    df["angle_score"] = df["BD_angle_deg"].apply(act_angle_score)
    df["s_pose"] = 100.0 * (0.30 * df["affinity_score"] + 0.20 * df["d1_score"] + 0.20 * df["d2_score"] + 0.30 * df["angle_score"])

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
