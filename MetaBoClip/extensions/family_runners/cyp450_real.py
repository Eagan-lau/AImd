from pathlib import Path
import math
import pandas as pd

from core.pymol_utils import get_pymol, parse_vina_out, split_ligand_states, dist, save_filtered_pse
from core.scoring import energy_linear, aggregate_standard
from core.output import write_rows_csv, write_score_tables


def _get_first_atom(cmd, selectors):
    for sel in selectors:
        if cmd.count_atoms(sel) > 0:
            a = cmd.get_model(sel).atom[0]
            return [a.coord[0], a.coord[1], a.coord[2]], int(a.id), sel
    return None, None, None


def _collect_heme_n(cmd):
    names = ["NA", "NB", "NC", "ND", "N1", "N2", "N3", "N4"]
    atoms = []
    for nm in names:
        sel = f"protein and resn HEM+HEME+HEC and name {nm}"
        if cmd.count_atoms(sel) > 0:
            a = cmd.get_model(sel).atom[0]
            atoms.append(([a.coord[0], a.coord[1], a.coord[2]], int(a.id), nm))
    if len(atoms) >= 3:
        return atoms[:3]
    sel = "protein and resn HEM+HEME+HEC and name N*"
    if cmd.count_atoms(sel) >= 3:
        model = cmd.get_model(sel)
        out = []
        for a in model.atom[:3]:
            out.append(([a.coord[0], a.coord[1], a.coord[2]], int(a.id), a.name))
        return out
    return None


def _plane_normal(a, b, c):
    v1 = [b[0]-a[0], b[1]-a[1], b[2]-a[2]]
    v2 = [c[0]-a[0], c[1]-a[1], c[2]-a[2]]
    n = [
        v1[1]*v2[2]-v1[2]*v2[1],
        v1[2]*v2[0]-v1[0]*v2[2],
        v1[0]*v2[1]-v1[1]*v2[0],
    ]
    nn = math.sqrt(sum(x*x for x in n)) + 1e-12
    return [x/nn for x in n]


def _oriented_normal(fe_xyz, normal, ref_xyz):
    ref_vec = [ref_xyz[0]-fe_xyz[0], ref_xyz[1]-fe_xyz[1], ref_xyz[2]-fe_xyz[2]]
    dot = sum(normal[i]*ref_vec[i] for i in range(3))
    if dot > 0:
        normal = [-x for x in normal]
    return normal


def _build_ch_pairs(cmd, pose_obj, h_cutoff=1.25):
    model = cmd.get_model(pose_obj)
    carbons = [a for a in model.atom if a.symbol == "C"]
    hydrogens = [a for a in model.atom if a.symbol == "H"]
    pairs = []
    for h in hydrogens:
        hxyz = [h.coord[0], h.coord[1], h.coord[2]]
        best = None
        for c in carbons:
            cxyz = [c.coord[0], c.coord[1], c.coord[2]]
            d = dist(hxyz, cxyz)
            if d <= h_cutoff and (best is None or d < best[0]):
                best = (d, c)
        if best is not None:
            c = best[1]
            pairs.append((int(h.id), [h.coord[0], h.coord[1], h.coord[2]], int(c.id), [c.coord[0], c.coord[1], c.coord[2]]))
    return pairs


def _axis_deviation(axis, vec):
    nv = math.sqrt(sum(x*x for x in vec)) + 1e-12
    cosv = max(-1.0, min(1.0, sum(axis[i]*vec[i] for i in range(3)) / nv))
    angle = math.degrees(math.acos(cosv))
    return min(angle, 180.0-angle)


def cyp450_real_runner(cfg, protein_dir, docking_dir, out_dir, registry):
    ligand_id = docking_dir.name.split("_")[-1]
    inp = cfg.get("input", {})
    protein_ext = inp.get("protein_extension", ".pdbqt")
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
    print(f"CYP450 real runner finished. Gating rows: {len(rows)}")
    print(f"Gating CSV: {gating_csv}")
    return {"rows": len(rows)}


def _process_one(cmd, cfg, protein_path, docking_dir, ligand_id, out_dir):
    rows = []
    protein_name = protein_path.stem
    inp = cfg["input"]
    ligand_file = docking_dir / inp["ligand_pattern"].format(
        ligand_id=ligand_id,
        protein_name=protein_name,
        ligand_extension=inp.get("ligand_extension", ".pdbqt"),
    )
    out_file = docking_dir / inp["out_pattern"].format(ligand_id=ligand_id, protein_name=protein_name)
    if not protein_path.exists() or not ligand_file.exists() or not out_file.exists():
        return rows

    modes, mode2aff = parse_vina_out(str(out_file), float(inp.get("affinity_threshold", -3.0)))
    if not modes:
        return rows

    cmd.reinitialize()
    cmd.load(str(protein_path), "protein")
    cmd.load(str(ligand_file), "ligand")

    fe_xyz, fe_id, _ = _get_first_atom(cmd, ["protein and elem Fe", "protein and name FE"])
    heme_n = _collect_heme_n(cmd)
    prox_sg_xyz, prox_sg_id, _ = _get_first_atom(cmd, [
        "protein and name SG within 4 of (protein and elem Fe)",
        "protein and name SG within 6 of (protein and elem Fe)",
    ])
    if fe_xyz is None or heme_n is None or prox_sg_xyz is None:
        return rows

    normal = _plane_normal(heme_n[0][0], heme_n[1][0], heme_n[2][0])
    axis = _oriented_normal(fe_xyz, normal, prox_sg_xyz)

    all_models = split_ligand_states(cmd, "ligand", "ligand_")
    retained = [f"ligand_{str(m).zfill(4)}" for m in modes if f"ligand_{str(m).zfill(4)}" in cmd.get_names("objects")]
    for m in all_models:
        if m not in retained:
            cmd.delete(m)

    passing = set()
    for model in retained:
        try:
            cmd.h_add(model)
        except Exception:
            pass
        pairs = _build_ch_pairs(cmd, model, 1.25)
        if not pairs:
            continue

        try:
            clash_sel = f"(polymer and not hydro) within {float(inp.get('clash_cutoff', 2.0))} of ({model} and not hydro)"
            if cmd.count_atoms(clash_sel) > 0:
                continue
        except Exception:
            pass

        mode_num = int(model.split("_")[1])
        affinity = mode2aff.get(mode_num, float("nan"))
        pose_any = False
        for hid, hxyz, cid, cxyz in pairs:
            d_fe_h = dist(fe_xyz, hxyz)
            d_fe_c = dist(fe_xyz, cxyz)
            axis_dev = _axis_deviation(axis, [cxyz[0]-fe_xyz[0], cxyz[1]-fe_xyz[1], cxyz[2]-fe_xyz[2]])
            if not (3.0 <= d_fe_h <= 5.0):
                continue
            if axis_dev > 45.0:
                continue
            rows.append({
                "Ligand_id": str(ligand_id),
                "protein_id": protein_name.split("_")[0],
                "conformation": protein_name,
                "mode": int(mode_num),
                "affinity_kcal": float(affinity),
                "dist_FE_H": float(d_fe_h),
                "dist_FE_C": float(d_fe_c),
                "axis_dev_deg": float(axis_dev),
                "lig_h_atom_id": int(hid),
                "lig_c_atom_id": int(cid),
                "fe_atom_id": int(fe_id),
                "prox_sg_atom_id": int(prox_sg_id),
                "heme_n1_id": int(heme_n[0][1]),
                "heme_n2_id": int(heme_n[1][1]),
                "heme_n3_id": int(heme_n[2][1]),
                "clash_flag": False,
            })
            pose_any = True
        if pose_any:
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

    def dist_score_h(d):
        if pd.isna(d):
            return 0.0
        x = float(d)
        if 3.0 <= x <= 3.8:
            return 1.0
        if x < 3.0:
            return max(0.0, (x - 2.5) / 0.5)
        if x > 5.0:
            return 0.0
        return (5.0 - x) / 1.2

    def angle_score(a):
        if pd.isna(a):
            return 0.0
        x = float(a)
        if x <= 15.0:
            return 1.0
        if x >= 45.0:
            return 0.0
        return (45.0 - x) / 30.0

    df["affinity_score"] = df["affinity_kcal"].apply(lambda x: energy_linear(x, -7.0, -3.0))
    df["dist_score"] = df["dist_FE_H"].apply(dist_score_h)
    df["angle_score"] = df["axis_dev_deg"].apply(angle_score)
    df["s_pose"] = 100.0 * (0.30 * df["affinity_score"] + 0.40 * df["dist_score"] + 0.30 * df["angle_score"])

    agg = cfg.get("aggregation", {})
    conf_df, prot_df = aggregate_standard(
        df,
        total_confs=int(agg.get("total_confs", 6)),
        cover_t=float(agg.get("cover_t", 70.0)),
        alpha=float(agg.get("alpha", 0.30)),
    )
    return {"pose_scores": df, "conformation_scores": conf_df, "protein_scores": prot_df}
