#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import math
import pymol2
import glob
import numpy as np
import csv
import argparse

# -----------------------------
# Tunable thresholds (gating)
# -----------------------------
AFFINITY_THRESHOLD = -3.0
ANGLE_THRESHOLD = 30.0          # angle to heme plane normal; we use min(theta, 180-theta)
FEH_DIST_MIN = 3.0              # Fe–H distance lower bound (Å)
FEH_DIST_MAX = 5.0              # Fe–H distance upper bound (Å)
CLASH_CUTOFF  = 2.0             # Å
CH_NEAR_CUTOFF = 1.25           # Å, judge if H belongs to a carbon
SG_SEARCH_RADIUS = 5.0          # Å, proximal SG search radius

# -----------------------------
# Object / residue naming
# -----------------------------
PROT_OBJ = "protein1"           # 使用 protein1 作为读入蛋白的对象名
HEM_RESN = "HEM"                # 若有 HEC/HEA，可改为 "HEM+HEC+HEA"

# -----------------------------
# Paths
# -----------------------------
protein_dir = '/public/agis/yanjianbin_group/liuyugeng/taxus_alphaflow_docking/CYP450/data_input/data_protein/all'
docking_base_dir = '/public/agis/yanjianbin_group/liuyugeng/taxus_alphaflow_docking/CYP450/data_output/4_docking_results'
output_base_dir  = '/public/agis/yanjianbin_group/liuyugeng/taxus_alphaflow_docking/CYP450/data_output/6_conformational_analysis_20250811'

# -----------------------------
# CSV output (per-ligand file)
# -----------------------------
CSV_DIR = "/public/agis/yanjianbin_group/liuyugeng/taxus_alphaflow_docking/CYP450/data_output/6_CSV_OUTPUT/gating_csv"
os.makedirs(CSV_DIR, exist_ok=True)

CSV_HEADER = [
    "Ligand_id",        # file_num
    "protein_id",       # protein family ID (prefix before first underscore)
    "conformation",     # protein conformation (protein file name without .pdb)
    "mode",             # Vina pose index
    "affinity_kcal",    # docking affinity from .out
    "dist_FE_H",        # Å (Fe to selected H)
    "angle_deg",        # ° (effective angle to heme normal: min(theta, 180-theta))
    "H_atom_id",        # chosen H atom id (PyMOL)
    "C_atom_id",        # the parent carbon atom id (PyMOL)
    "fe_atom_id",       # Fe atom id (PyMOL)
    "dist_fe_sg_A",     # Fe–SG distance in Å (or "NA" if not found)
    "prox_sg_found",    # True/False
    "clash_flag"        # False if passed clash gating
]

# -----------------------------
# Utilities
# -----------------------------
def parse_docking_out(out_file, affinity_threshold=AFFINITY_THRESHOLD):
    """
    Parse Vina .out table. Read rows until affinity fails the threshold.
    Return: (list1, mode2aff)
      - list1: modes passing the energy threshold (in order)
      - mode2aff: dict mapping mode -> affinity
    """
    list1 = []
    mode2aff = {}
    with open(out_file, 'r') as f:
        lines = f.readlines()

    table_started = False
    for line in lines:
        s = line.strip()
        if s.startswith("mode"):
            table_started = True
            continue
        if not table_started:
            continue
        if s.startswith("-----") or s == "":
            continue
        parts = s.split()
        if len(parts) < 2:
            continue
        try:
            mode_number = int(parts[0])
            affinity = float(parts[1])
        except ValueError:
            continue

        mode2aff[mode_number] = affinity
        if affinity < affinity_threshold:
            list1.append(mode_number)
        else:
            break
    return list1, mode2aff


def compute_heme_normal(cmd, fe_coords, sg_coord=None):
    """
    Fit the heme plane using HEM N* atoms (SVD), take the normal vector,
    and orient it to point away from proximal Cys SG if sg_coord provided.
    Return normalized [x,y,z] or None if cannot be computed.
    """
    heme_model = cmd.get_model(f"{PROT_OBJ} and resn {HEM_RESN} and name N*")
    if len(heme_model.atom) < 3:
        return None

    coords = np.array([atom.coord for atom in heme_model.atom])
    centroid = np.mean(coords, axis=0)
    centered = coords - centroid
    _, _, vt = np.linalg.svd(centered)
    normal = vt[2, :]
    n = np.linalg.norm(normal)
    if n == 0:
        return None
    normal = normal / n

    if sg_coord is not None:
        fe_np = np.array(fe_coords)
        vec_to_s = np.array(sg_coord) - fe_np
        vn = np.linalg.norm(vec_to_s)
        if vn > 0:
            vec_to_s = vec_to_s / vn
            if np.dot(normal, vec_to_s) > 0:
                normal = -normal
            print("Heme plane normal oriented using proximal Cys SG (5 Å search).")
    else:
        print("Proximal Cys SG not found within 5 Å; normal orientation not adjusted.")

    return normal.tolist()


def find_fe_and_sg(cmd):
    """
    Find Fe in heme and proximal SG within SG_SEARCH_RADIUS.
    Return: (fe_atom_id, fe_coords(numpy array), sg_found(bool), fe_sg_dist(float or 'NA'), sg_coord or None)
    """
    fe_model = cmd.get_model(f"{PROT_OBJ} and name FE")
    if not fe_model.atom:
        return None, None, False, "NA", None
    fe_atom = fe_model.atom[0]
    fe_coords = np.array(fe_atom.coord)

    # SG search within SG_SEARCH_RADIUS
    tmp = "_fe_tmp_for_sg"
    cmd.pseudoatom(tmp, pos=fe_coords.tolist())
    try:
        s_model = cmd.get_model(f"({PROT_OBJ} and resn CYS and name SG) within {SG_SEARCH_RADIUS} of {tmp}")
    finally:
        cmd.delete(tmp)

    if not s_model.atom:
        return int(fe_atom.id), fe_coords, False, "NA", None

    # choose nearest SG
    best_at = min(s_model.atom, key=lambda at: np.linalg.norm(np.array(at.coord) - fe_coords))
    fe_sg_dist = float(np.linalg.norm(np.array(best_at.coord) - fe_coords))
    return int(fe_atom.id), fe_coords, True, round(fe_sg_dist, 3), np.array(best_at.coord)


def nearest_parent_carbon_id(cmd, obj_name, h_id, cutoff=CH_NEAR_CUTOFF):
    """
    Return the nearest carbon atom id that owns this H (within cutoff), or None.
    """
    sel = f"({obj_name} and elem C) within {cutoff} of id {h_id}"
    cands = cmd.get_model(sel)
    if not cands.atom:
        return None
    if len(cands.atom) == 1:
        return int(cands.atom[0].id)
    # choose nearest C if more than one (rare)
    h_coord = np.array(cmd.get_model(f"id {h_id}").atom[0].coord)
    best = min(cands.atom, key=lambda at: np.linalg.norm(np.array(at.coord) - h_coord))
    return int(best.id)


# ========= NEW: fast H->C parent mapping (minimal change; logic-equivalent) =========
def build_CH_parent_map_fast(cmd, obj_name, cutoff=CH_NEAR_CUTOFF):
    """
    用一次 get_model 拿到该 pose 的所有 H 与 C 坐标，在 NumPy 中做成对距离，
    对每个 H 选出 <= cutoff Å 的最近 C；若没有满足者 -> None。
    返回: dict {H_id -> C_id 或 None}
    （功能与 nearest_parent_carbon_id 等价，但避免对每个 H 触发一次选择器）
    """
    m = cmd.get_model(obj_name)
    if not m.atom:
        return {}

    H_ids, H_xyz = [], []
    C_ids, C_xyz = [], []
    for a in m.atom:
        if a.symbol == 'H':
            H_ids.append(int(a.id)); H_xyz.append(a.coord)
        elif a.symbol == 'C':
            C_ids.append(int(a.id)); C_xyz.append(a.coord)

    if not H_ids or not C_ids:
        return {hid: None for hid in H_ids}

    H_xyz = np.array(H_xyz, dtype=float)   # (nH,3)
    C_xyz = np.array(C_xyz, dtype=float)   # (nC,3)

    # pairwise distances (nH, nC)
    diff = H_xyz[:, None, :] - C_xyz[None, :, :]
    dist = np.sqrt(np.sum(diff * diff, axis=2))

    near_mask = dist <= float(cutoff)
    mapping = {}
    for i, hid in enumerate(H_ids):
        idxs = np.nonzero(near_mask[i])[0]
        if idxs.size == 0:
            mapping[hid] = None
        else:
            j = idxs[np.argmin(dist[i, idxs])]
            mapping[hid] = int(C_ids[j])
    return mapping
# ================================================================================


def expand_file_range_str(range_str):
    """
    Parse a range string like:
      '1-10' -> [1,2,...,10]
      '3,5,8' -> [3,5,8]
      '1-3,7,10-12' -> [1,2,3,7,10,11,12]
    Return a sorted list of unique ints. Raise ValueError on bad input.
    """
    nums = set()
    for part in range_str.split(','):
        part = part.strip()
        if not part:
            continue
        if '-' in part:
            a, b = part.split('-', 1)
            a = int(a.strip()); b = int(b.strip())
            if a > b:
                a, b = b, a
            nums.update(range(a, b + 1))
        else:
            nums.add(int(part))
    return sorted(nums)

# -----------------------------
# Core processing
# -----------------------------
def process_protein(cmd, protein_name, protein_dir, docking_dir, output_dir, file_num,
                    csv_writer,
                    angle_threshold=ANGLE_THRESHOLD, affinity_threshold=AFFINITY_THRESHOLD,
                    feh_min=FEH_DIST_MIN, feh_max=FEH_DIST_MAX, clash_cutoff=CLASH_CUTOFF):
    """
    One protein conformation:
      1) Parse {file_num}_{protein_name}_cavity_1.out for passing modes
      2) Load protein (object name = PROT_OBJ) and ligand; split ligand poses
      3) Keep only passing modes
      4) For each kept mode:
         - add hydrogens
         - find C–H pairs (H within 1.25 Å of a carbon)
         - gate by Fe–H distance [feh_min, feh_max] and angle to heme normal
         - if a carbon has multiple H passing gate, keep only the H closest to Fe
         - clash gate vs polymer heavy atoms
         - write CSV rows with (H_id, C_id); save a sele including all passing H for that pose
      5) Save .pse if any valid mode remains; write angles.txt summary
    """
    out_file = os.path.join(docking_dir, f'{file_num}_' + protein_name + '_cavity_1.out')
    if not os.path.exists(out_file):
        print("Docking .out not found:", out_file)
        return

    list1, mode2aff = parse_docking_out(out_file, affinity_threshold=affinity_threshold)
    if not list1:
        print(f"{protein_name}: no modes passed affinity threshold.")
        return

    # Load protein (object name = PROT_OBJ)
    protein_path = os.path.join(protein_dir, protein_name + '.pdb')
    if not os.path.exists(protein_path):
        print("Protein PDB not found:", protein_path)
        return
    cmd.load(protein_path, PROT_OBJ)

    # Load ligand
    ligand_path = os.path.join(docking_dir, f'{file_num}@' + protein_name + '.pdbqt')
    if not os.path.exists(ligand_path):
        print("Ligand file not found:", ligand_path)
        cmd.delete("all")
        return
    cmd.load(ligand_path, "ligand")

    # Split ligand states -> pose_1, pose_2, ...
    cmd.split_states("ligand", prefix="object_ligand")
    cmd.delete("ligand")
    objs = cmd.get_object_list("object_ligand*")
    objs.sort(key=lambda x: int(x.replace('object_ligand', '')))
    for i, obj in enumerate(objs, start=1):
        cmd.set_name(obj, f"pose_{i}")

    # Energy filter: delete modes not in list1
    total_states = len(objs)
    for mode in range(1, total_states + 1):
        name = f"pose_{mode}"
        if (mode not in list1) and (name in cmd.get_object_list("all")):
            cmd.delete(name)

    # Find Fe and proximal SG (within 5 Å), compute heme normal (flip by SG if found)
    fe_atom_id, fe_coords, sg_found, fe_sg_dist, sg_coord = find_fe_and_sg(cmd)
    if fe_atom_id is None:
        print(f"{protein_name}: Fe atom not found.")
        return

    heme_normal = compute_heme_normal(cmd, fe_coords, sg_coord=sg_coord)
    if heme_normal is None:
        print(f"{protein_name}: cannot compute heme plane normal; angle gating will be skipped.")

    valid_ligand_objects = []
    protein_angles_info = []  # [(mode_name, [(feh_dist, eff_angle, H_id, C_id), ...], affinity), ...]

    # Iterate passing modes
    for mode in list1:
        mode_name = f"pose_{mode}"
        if mode_name not in cmd.get_object_list("all"):
            continue

        cmd.h_add(mode_name)

        # --- build H->C parent mapping once per pose (fast, logic-equivalent) ---
        CH_map = build_CH_parent_map_fast(cmd, mode_name, cutoff=CH_NEAR_CUTOFF)

        # --- find candidate H atoms belonging to carbons (within 1.25 Å) ---
        Hs = cmd.get_model(f"{mode_name} and elem H").atom
        if not Hs:
            # no hydrogens in this pose after h_add?
            cmd.delete(mode_name)
            continue

        # Per-carbon best H (closest Fe) among those passing gate
        best_H_for_C = {}  # C_id -> (H_id, feh_dist, eff_angle)

        for h_atom in Hs:
            H_id = int(h_atom.id)
            # use fast mapping (minimal change: drop per-H selector)
            C_id = CH_map.get(H_id, None)
            if C_id is None:
                continue

            hx, hy, hz = h_atom.coord
            # Fe–H distance
            dx = fe_coords[0] - hx
            dy = fe_coords[1] - hy
            dz = fe_coords[2] - hz
            feh_dist = math.sqrt(dx*dx + dy*dy + dz*dz)

            if not (FEH_DIST_MIN <= feh_dist <= FEH_DIST_MAX):
                continue

            # angle between Fe->H and heme normal
            if heme_normal is not None:
                v = [hx - fe_coords[0], hy - fe_coords[1], hz - fe_coords[2]]
                v_norm = math.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])
                if v_norm == 0:
                    continue
                dot = (v[0]*heme_normal[0] + v[1]*heme_normal[1] + v[2]*heme_normal[2]) / v_norm
                dot = max(min(dot, 1.0), -1.0)
                angle = math.degrees(math.acos(dot))
                eff_angle = min(angle, 180 - angle)
            else:
                eff_angle = None

            if eff_angle is None or eff_angle > ANGLE_THRESHOLD:
                continue

            # keep only the closest H to Fe for this carbon
            prev = best_H_for_C.get(C_id)
            if (prev is None) or (feh_dist < prev[1]):
                best_H_for_C[C_id] = (H_id, feh_dist, eff_angle)

        # --- collect chosen Hs (one per carbon) ---
        chosen = []
        for C_id, (H_id, feh_dist, eff_angle) in best_H_for_C.items():
            chosen.append((H_id, C_id, feh_dist, eff_angle))

        # clash check (only if we have at least one candidate)
        within_range = len(chosen) > 0
        if within_range and heme_normal is not None:
            clash_sel = f"(polymer and not hydro) within {CLASH_CUTOFF} of ({mode_name} and not hydro)"
            clash_count = cmd.count_atoms(clash_sel)
            if clash_count > 0:
                within_range = False
                print(f"{mode_name}: clash detected (count={clash_count}) -> invalid.")

        # save info & CSV rows
        if within_range and chosen:
            aff = mode2aff.get(mode, float('nan'))
            # record for angles report
            info_list = []
            # create a PyMOL selection for these H atoms (saved in .pse)
            h_ids = [hid for (hid, cid, d, ang) in chosen]
            if h_ids:
                sel_expr = " or ".join([f"id {hid}" for hid in h_ids])
                cmd.select(f"CH_hits_{mode}", f"{mode_name} and ({sel_expr})")

            for (H_id, C_id, feh_dist, eff_angle) in chosen:
                info_list.append((feh_dist, eff_angle, H_id, C_id))
                # write CSV (one row per C–H kept)
                csv_writer.writerow([
                    file_num,
                    protein_name.split('_')[0],
                    protein_name,
                    mode,
                    float(aff),
                    float(feh_dist),
                    float(eff_angle),
                    int(H_id),
                    int(C_id),
                    int(fe_atom_id),
                    (f"{fe_sg_dist:.3f}" if sg_found else "NA"),
                    bool(sg_found),
                    False
                ])

            protein_angles_info.append((mode_name, info_list, aff))
            valid_ligand_objects.append(mode_name)
        else:
            # delete the pose if no passing C–H
            if mode_name in cmd.get_object_list("all"):
                cmd.delete(mode_name)

    # Write angles summary
    angles_filename = os.path.join(output_dir, protein_name + "_angles.txt")
    try:
        with open(angles_filename, "w") as f_angles:
            f_angles.write("Parameters:\n")
            f_angles.write(f"Affinity Threshold: {AFFINITY_THRESHOLD}\n")
            f_angles.write(f"Angle Threshold: {ANGLE_THRESHOLD}°\n")
            f_angles.write(f"Fe–H Distance Range: {FEH_DIST_MIN} - {FEH_DIST_MAX} Å\n")
            f_angles.write(f"C–H ownership cutoff: {CH_NEAR_CUTOFF} Å\n")
            f_angles.write(f"Clash Cutoff: {CLASH_CUTOFF} Å\n")
            f_angles.write(f"Proximal SG found (≤ {SG_SEARCH_RADIUS} Å): {sg_found}\n")
            f_angles.write(f"Fe–SG distance: {(f'{fe_sg_dist:.3f} Å' if sg_found else 'NA')}\n\n")

            for mode_name, info, aff in protein_angles_info:
                f_angles.write(f"Mode: {mode_name} | affinity: {aff:.3f} kcal/mol\n")
                if info:
                    for (feh_dist, eff_angle, H_id, C_id) in info:
                        f_angles.write(
                            f"  Fe–H: {feh_dist:.2f} Å, Effective angle: {eff_angle:.2f}°, "
                            f"H_id: {H_id}, C_id: {C_id}\n"
                        )
                else:
                    f_angles.write("  No candidate C–H met the criteria.\n")
                f_angles.write("\n")
        print("Angles info written:", angles_filename)
    except Exception as e:
        print("Failed to write angles info:", e)

    # Save .pse if any valid mode remains
    if valid_ligand_objects:
        out_session = os.path.join(output_dir, protein_name + ".pse")
        cmd.save(out_session)
        print(f"{protein_name}: session saved -> {out_session}")
    else:
        print(f"{protein_name}: no ligand modes passed gating.")


def main(cmd, file_range_str=None):
    # 确定要处理的 file_N 列表
    if file_range_str:
        try:
            file_nums = expand_file_range_str(file_range_str)
        except Exception as e:
            raise SystemExit(f"Invalid --file-range: {file_range_str} ({e})")
        # 仅保留实际存在的子目录
        file_nums = [n for n in file_nums if os.path.isdir(os.path.join(docking_base_dir, f"file_{n}"))]
        if not file_nums:
            print(f"No matching subfolders under {docking_base_dir} for range '{file_range_str}'. Exit.")
            return
    else:
        # 默认：处理所有子目录
        docking_files = glob.glob(os.path.join(docking_base_dir, "file_*"))
        file_nums = sorted([
            int(os.path.basename(f).split('_')[-1])
            for f in docking_files
            if os.path.isdir(f) and os.path.basename(f).split('_')[-1].isdigit()
        ])

    print(f"Processing file_nums: {file_nums}")

    for file_num in file_nums:
        print(f"\nProcessing file_num: {file_num}")
        docking_dir = f"{docking_base_dir}/file_{file_num}"
        output_dir  = f"{output_base_dir}/file_{file_num}"
        os.makedirs(output_dir, exist_ok=True)

        # per-ligand CSV
        csv_path_this = os.path.join(CSV_DIR, f"file_{file_num}.csv")
        with open(csv_path_this, "w", newline="") as csv_fh:
            csv_writer = csv.writer(csv_fh)
            csv_writer.writerow(CSV_HEADER)

            for file in os.listdir(protein_dir):
                if file.endswith('.pdb'):
                    protein_name = file[:-4]
                    print("Protein:", protein_name)
                    process_protein(cmd, protein_name, protein_dir, docking_dir, output_dir, file_num, csv_writer)
                    cmd.delete("all")

        print(f"CSV written: {csv_path_this}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="CYP450 C–H gating with Fe–SG info and selectable file_* range.")
    parser.add_argument("--file-range", type=str, default=None,
                        help="e.g., '1-10' or '3,5,8' or '1-3,7,10-12'. If omitted, process all file_* under docking_base_dir.")
    args = parser.parse_args()

    with pymol2.PyMOL() as pymol:
        cmd = pymol.cmd
        main(cmd, file_range_str=args.file_range)
