from pathlib import Path
import pandas as pd

from core.pymol_utils import get_pymol, parse_vina_out, split_ligand_states, angle_three_points, dist, save_filtered_pse
from core.scoring import energy_linear, distance_piecewise, gaussian_score, aggregate_standard
from core.output import write_rows_csv, write_score_tables

PROJECT_ROOT = Path(__file__).resolve().parents[2]


def _resolve_project_path(raw_path):
    path = Path(raw_path)
    if path.is_absolute():
        return path
    candidate = PROJECT_ROOT / path
    return candidate if candidate.exists() else path


def _format_selector(template, **selector_vars):
    if not isinstance(template, str):
        return template
    return template.format(**selector_vars)


def _align_template_to_protein(cmd, template_cfg):
    method = str(template_cfg.get("align_method", "super")).lower()
    mobile = _format_selector(
        template_cfg.get("mobile_selector", "{template_obj} and name CA"),
        template_obj="A1",
        protein_obj="prot",
    )
    target = _format_selector(
        template_cfg.get("target_selector", "{protein_obj} and name CA"),
        template_obj="A1",
        protein_obj="prot",
    )
    if cmd.count_atoms(mobile) <= 0 or cmd.count_atoms(target) <= 0:
        return False
    if method == "align":
        cmd.align(mobile, target)
    elif method == "cealign":
        cmd.cealign(target, mobile)
    else:
        cmd.super(mobile, target)
    return True

def _sel_nonempty(cmd, selname, expr):
    cmd.select(selname, expr)
    return cmd.count_atoms(selname) > 0

def _pick_unique_ser_and_his(cmd):
    if not _sel_nonempty(cmd, 'ser171', 'A1 and resn SER and resi 171'):
        return None
    if not _sel_nonempty(cmd, 'ser171_ca', 'ser171 and name CA'):
        return None
    ser171_ca_xyz = cmd.get_coords('ser171_ca', 1)[0]
    if not _sel_nonempty(cmd, 'near_res', 'prot and byres (ser171 around 2)'):
        return None
    if not _sel_nonempty(cmd, 'ser_near', 'near_res and resn SER'):
        return None
    ser_list = []
    cmd.iterate('ser_near and name CA', 'ser_list.append((chain, resi))', space={'ser_list': ser_list})
    if not ser_list:
        return None
    best_ser = None
    for ser_chain, ser_resi in ser_list:
        expr = f"(prot and {'chain '+ser_chain+' and ' if ser_chain and ser_chain.strip() else ''}resi {ser_resi} and resn SER and name CA)"
        if not _sel_nonempty(cmd, 'ser_ca_tmp', expr):
            continue
        d = dist(cmd.get_coords('ser_ca_tmp', 1)[0], ser171_ca_xyz)
        if best_ser is None or d < best_ser[0]:
            best_ser = (d, ser_chain, ser_resi)
    if best_ser is None:
        return None
    _, ser_chain, ser_resi = best_ser
    current_ser_expr = f"(prot and {'chain '+ser_chain+' and ' if ser_chain and ser_chain.strip() else ''}resi {ser_resi} and resn SER)"
    if not _sel_nonempty(cmd, 'current_ser', current_ser_expr):
        return None
    if not _sel_nonempty(cmd, 'ser_og', 'current_ser and name OG'):
        return None
    ser_og_xyz = cmd.get_coords('ser_og', 1)[0]
    if not _sel_nonempty(cmd, 'around_ser', f'prot and byres ({current_ser_expr} around 5)'):
        return None
    if not _sel_nonempty(cmd, 'his_near', 'around_ser and resn HIS'):
        return None
    his_list = []
    cmd.iterate('his_near and name CA', 'his_list.append((chain, resi))', space={'his_list': his_list})
    if not his_list:
        return None
    best_his = None
    for his_chain, his_resi in his_list:
        expr = f"(prot and {'chain '+his_chain+' and ' if his_chain and his_chain.strip() else ''}resi {his_resi} and name NE2)"
        if not _sel_nonempty(cmd, 'his_ne2', expr):
            continue
        try:
            d_sh = cmd.get_distance('ser_og', 'his_ne2')
        except Exception:
            continue
        if d_sh <= 5.0 and (best_his is None or d_sh < best_his[0]):
            best_his = (d_sh, his_chain, his_resi)
    if best_his is None:
        return None
    d_sh, his_chain, his_resi = best_his
    return (ser_chain, ser_resi, ser_og_xyz, his_chain, his_resi, d_sh)

def ach_real_runner(cfg, protein_dir, docking_dir, out_dir, registry):
    ligand_id = docking_dir.name.split("_")[-1]
    protein_ext = cfg.get("input", {}).get("protein_extension", ".pdbqt")
    proteins = sorted(protein_dir.glob(f"*{protein_ext}"))
    rows = []
    cmd, pm = get_pymol()
    try:
        for protein_path in proteins:
            rows.extend(_process_one(cmd, cfg, protein_path, docking_dir, ligand_id, out_dir))
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
    return {"rows": len(rows)}

def _process_one(cmd, cfg, protein_path, docking_dir, ligand_id, out_dir):
    rows = []
    protein_name = protein_path.stem
    inp = cfg["input"]
    ligand_file = docking_dir / inp["ligand_pattern"].format(ligand_id=ligand_id, protein_name=protein_name, ligand_extension=inp.get("ligand_extension",".pdbqt"))
    out_file = docking_dir / inp["out_pattern"].format(ligand_id=ligand_id, protein_name=protein_name)
    template_cfg = (
        cfg.get("resources", {}).get("template")
        or cfg.get("protein_side", {}).get("template")
        or {"path": "examples/templates/dAC-T.pdb"}
    )
    template_path = _resolve_project_path(template_cfg["path"])
    if not protein_path.exists() or not ligand_file.exists() or not out_file.exists() or not template_path.exists():
        return rows
    modes, mode2aff = parse_vina_out(str(out_file), float(inp.get("affinity_threshold", -3.0)))
    if not modes:
        return rows
    cmd.reinitialize()
    cmd.load(str(protein_path), 'prot')
    cmd.load(str(template_path), 'A1')
    if template_cfg.get("align_to_protein", True):
        if not _align_template_to_protein(cmd, template_cfg):
            return rows
    picked = _pick_unique_ser_and_his(cmd)
    if picked is None:
        return rows
    ser_chain, ser_resi, ser_og_xyz, his_chain, his_resi, dist_ser_his = picked
    his_ne2_expr = f"(prot and {'chain '+his_chain+' and ' if his_chain and his_chain.strip() else ''}resi {his_resi} and name NE2)"
    if not _sel_nonempty(cmd, 'his_ne2', his_ne2_expr):
        return rows
    cmd.load(str(ligand_file), 'ligand')
    all_models = split_ligand_states(cmd, 'ligand', 'ligand_')
    retained = [f"ligand_{str(int(m)).zfill(4)}" for m in modes if f"ligand_{str(int(m)).zfill(4)}" in cmd.get_names('objects')]
    for m in all_models:
        if m not in retained:
            cmd.delete(m)
    passing = set()
    for obj in retained:
        mode_num = int(obj.replace('ligand_', '').lstrip('0') or '1')
        aff = mode2aff.get(mode_num, float('nan'))
        if not _sel_nonempty(cmd, 'carbon_with_o', f"{obj} and elem C and neighbor ({obj} and elem O)"):
            continue
        model = cmd.get_model(obj)
        oxygens = [b for b in model.atom if b.symbol == 'O']
        ligand_atom_ids = ",".join(str(a.index) for a in model.atom)
        carbons = []
        for idx, a in enumerate(cmd.get_model('carbon_with_o').atom, start=1):
            carbon_index = int(a.index)
            carbon_id = int(a.id)
            cxyz = [a.coord[0], a.coord[1], a.coord[2]]
            o_neighbors = []
            for b in oxygens:
                bxyz = [b.coord[0], b.coord[1], b.coord[2]]
                if dist(cxyz, bxyz) <= 2.0:
                    o_neighbors.append((int(b.id), bxyz))
            if len(o_neighbors) != 2:
                continue
            nearest_o_id, nearest_o_xyz = min(o_neighbors, key=lambda item: dist(cxyz, item[1]))
            carbons.append((idx, carbon_index, carbon_id, cxyz, nearest_o_id, nearest_o_xyz))
        pose_any = False
        for idx, carbon_index, cid, cxyz, oid, oxy in carbons:
            d_c_ser = dist(cxyz, ser_og_xyz)
            if d_c_ser > 5.0:
                continue
            sele_name = f"{obj}_Csel_{idx}"
            _sel_nonempty(cmd, sele_name, f"{obj} and index {carbon_index}")
            bd = angle_three_points(ser_og_xyz, cxyz, oxy)
            rows.append({
                "Ligand_id": str(ligand_id),
                "protein_id": protein_name.split("_")[0],
                "conformation": protein_name,
                "mode": int(mode_num),
                "affinity_kcal": float(aff),
                "carbon_index": int(carbon_index),
                "carbon_id": int(cid),
                "selection_name": sele_name,
                "dist_to_SerOG_A": float(d_c_ser),
                "dist_SerOG_HisNE2_A": float(dist_ser_his),
                "bd_angle_deg": float(bd),
                "ligand_atom_ids": ligand_atom_ids,
            })
            pose_any = True
        if pose_any:
            passing.add(obj)
    if passing and cfg.get("output", {}).get("save_pse", True):
        sessions_dir = out_dir / "sessions"
        sessions_dir.mkdir(parents=True, exist_ok=True)
        save_filtered_pse(cmd, passing, sessions_dir / f"{protein_name}_filtered.pse")
    return rows

def _score_rows(rows, cfg):
    df = pd.DataFrame(rows)
    if df.empty:
        return {"pose_scores": pd.DataFrame(), "conformation_scores": pd.DataFrame(), "protein_scores": pd.DataFrame()}
    df["affinity_score"] = df["affinity_kcal"].apply(lambda x: energy_linear(x, -7.0, -3.0))
    df["d1_score"] = df["dist_to_SerOG_A"].apply(lambda x: distance_piecewise(x, 3.2, 5.0))
    df["d2_score"] = df["dist_SerOG_HisNE2_A"].apply(lambda x: distance_piecewise(x, 3.2, 5.0))
    df["D_score"] = 0.90 * df["d1_score"] + 0.10 * df["d2_score"]
    df["angle_score"] = df["bd_angle_deg"].apply(lambda x: gaussian_score(x, 107.0, 20.0, flat=0.0, fold180=False))
    df["s_pose"] = 100.0 * (0.30 * df["affinity_score"] + 0.50 * df["D_score"] + 0.20 * df["angle_score"])
    agg = cfg.get("aggregation", {})
    conf_df, prot_df = aggregate_standard(
        df,
        total_confs=int(agg.get("total_confs", 6)),
        cover_t=float(agg.get("cover_t", 70.0)),
        alpha=float(agg.get("alpha", 0.30)),
    )
    return {"pose_scores": df, "conformation_scores": conf_df, "protein_scores": prot_df}
