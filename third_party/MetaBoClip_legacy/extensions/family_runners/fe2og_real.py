from pathlib import Path
import pandas as pd
from core.pymol_utils import get_pymol, parse_vina_out, split_ligand_states, dist, save_filtered_pse
from core.scoring import energy_linear, aggregate_standard
from core.output import write_rows_csv, write_score_tables
import math

def _get_atom(cmd, sel):
    if cmd.count_atoms(sel) == 0:
        return None, None
    a = cmd.get_model(sel).atom[0]
    return [a.coord[0], a.coord[1], a.coord[2]], int(a.id)

def _get_anchorA(cmd):
    for name in ["O1","O2","C1","C2"]:
        sel = f"protein and resn AKG and name {name}"
        if cmd.count_atoms(sel) > 0:
            a = cmd.get_model(sel).atom[0]
            return [a.coord[0], a.coord[1], a.coord[2]], int(a.id), name
    return None, None, None

def _plane_normal(fe, o5, a):
    v1 = [o5[0]-fe[0], o5[1]-fe[1], o5[2]-fe[2]]
    v2 = [a[0]-fe[0], a[1]-fe[1], a[2]-fe[2]]
    n = [
        v1[1]*v2[2]-v1[2]*v2[1],
        v1[2]*v2[0]-v1[0]*v2[2],
        v1[0]*v2[1]-v1[1]*v2[0],
    ]
    nn = math.sqrt(n[0]**2+n[1]**2+n[2]**2) + 1e-12
    return [n[0]/nn, n[1]/nn, n[2]/nn]

def _build_ch_map(cmd, pose_obj, cutoff=1.25):
    model = cmd.get_model(pose_obj)
    carbons = [a for a in model.atom if a.symbol == 'C']
    hydrogens = [a for a in model.atom if a.symbol == 'H']
    mp = {}
    for h in hydrogens:
        hxyz = [h.coord[0], h.coord[1], h.coord[2]]
        best = None
        for c in carbons:
            cxyz = [c.coord[0], c.coord[1], c.coord[2]]
            d = dist(hxyz, cxyz)
            if d <= cutoff and (best is None or d < best[0]):
                best = (d, int(c.id))
        if best is not None:
            mp[int(h.id)] = best[1]
    return mp

def fe2og_real_runner(cfg, protein_dir, docking_dir, out_dir, registry):
    ligand_id = docking_dir.name.split("_")[-1]
    protein_ext = cfg.get("input", {}).get("protein_extension", ".pdb")
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
    if not protein_path.exists() or not ligand_file.exists() or not out_file.exists():
        return rows
    modes, mode2aff = parse_vina_out(str(out_file), float(inp.get("affinity_threshold", -3.0)))
    if not modes:
        return rows
    cmd.reinitialize()
    cmd.load(str(protein_path), "protein")
    cmd.load(str(ligand_file), "ligand")
    all_models = split_ligand_states(cmd, "ligand", "ligand_")
    retained = [f"ligand_{str(m).zfill(4)}" for m in modes if f"ligand_{str(m).zfill(4)}" in cmd.get_names('objects')]
    for m in all_models:
        if m not in retained:
            cmd.delete(m)
    fe, fe_id = _get_atom(cmd, "protein and name FE")
    if fe is None:
        fe, fe_id = _get_atom(cmd, "protein and elem Fe")
    o5, o5_id = _get_atom(cmd, "protein and resn AKG and name O5")
    a, a_id, a_tag = _get_anchorA(cmd)
    if fe is None or o5 is None or a is None:
        return rows
    normal = _plane_normal(fe, o5, a)
    passing = set()
    for model in retained:
        cmd.h_add(model)
        ch_map = _build_ch_map(cmd, model, 1.25)
        Hs = cmd.get_model(f"{model} and elem H").atom
        if not Hs:
            continue
        pose_rows = []
        for h in Hs:
            hid = int(h.id); cid = ch_map.get(hid)
            if cid is None:
                continue
            hxyz = [h.coord[0], h.coord[1], h.coord[2]]
            feh = dist(fe, hxyz)
            if not (2.5 <= feh <= 5.0):
                continue
            v = [hxyz[0]-fe[0], hxyz[1]-fe[1], hxyz[2]-fe[2]]
            vn = math.sqrt(v[0]**2 + v[1]**2 + v[2]**2) + 1e-12
            cosv = max(-1.0, min(1.0, (v[0]*normal[0] + v[1]*normal[1] + v[2]*normal[2]) / vn))
            angle = math.degrees(math.acos(cosv))
            axis_dev = min(angle, 180-angle)
            if axis_dev > 60.0:
                continue
            pose_rows.append((hid, cid, feh, axis_dev))
        if not pose_rows:
            continue
        clash_sel = f"(polymer and not hydro) within 2.0 of ({model} and not hydro)"
        try:
            if cmd.count_atoms(clash_sel) > 0:
                continue
        except Exception:
            pass
        mode_num = int(model.split("_")[1]); aff = mode2aff.get(mode_num, float("nan"))
        for hid, cid, feh, axis_dev in pose_rows:
            rows.append({
                "Ligand_id": str(ligand_id),
                "protein_id": protein_name.split("_")[0],
                "conformation": protein_name,
                "mode": int(mode_num),
                "affinity_kcal": float(aff),
                "dist_FE_H": float(feh),
                "angle_deg": float(axis_dev),
                "lig_atom_id": int(hid),
                "parent_atom_id": int(cid),
                "fe_atom_id": int(fe_id),
                "o_ref_atom_id": int(o5_id),
                "anchor_atom_id": int(a_id),
                "anchor_tag": a_tag,
                "clash_flag": False,
            })
        passing.add(model)
    if passing and cfg.get("output", {}).get("save_pse", True):
        sessions_dir = out_dir / "sessions"
        sessions_dir.mkdir(parents=True, exist_ok=True)
        save_filtered_pse(cmd, passing, sessions_dir / f"{protein_name}.pse")
    return rows

def _score_rows(rows, cfg):
    df = pd.DataFrame(rows)
    if df.empty:
        return {"pose_scores": pd.DataFrame(), "conformation_scores": pd.DataFrame(), "protein_scores": pd.DataFrame()}
    def dist_score(d):
        if pd.isna(d): return 0.0
        x = float(d)
        if x <= 3.5: return 1.0
        if x >= 5.0: return 0.0
        return (5.0 - x) / 1.5
    def angle_score(a):
        if pd.isna(a): return 0.0
        x = float(a)
        if x <= 30.0: return 1.0
        if x >= 60.0: return 0.0
        return (60.0 - x) / 30.0
    df["affinity_score"] = df["affinity_kcal"].apply(lambda x: energy_linear(x, -7.0, -3.0))
    df["dist_score"] = df["dist_FE_H"].apply(dist_score)
    df["angle_score"] = df["angle_deg"].apply(angle_score)
    df["s_pose"] = 100.0 * (0.30*df["affinity_score"] + 0.50*df["dist_score"] + 0.20*df["angle_score"])
    agg = cfg.get("aggregation", {})
    conf_df, prot_df = aggregate_standard(
        df,
        total_confs=int(agg.get("total_confs", 6)),
        cover_t=float(agg.get("cover_t", 70.0)),
        alpha=float(agg.get("alpha", 0.30)),
    )
    return {"pose_scores": df, "conformation_scores": conf_df, "protein_scores": prot_df}
