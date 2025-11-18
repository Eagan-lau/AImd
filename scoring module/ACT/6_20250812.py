#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import glob
import argparse
import pymol
from pymol import cmd
import numpy as np
import csv

# -------------------------
# Base directories
# -------------------------
docking_base_dir = "/public/agis/yanjianbin_group/liuyugeng/taxus_alphaflow_docking/ACT/yixianji/data_output/4_docking_results"

# Protein directory
protein_dir = "/public/agis/yanjianbin_group/liuyugeng/taxus_alphaflow_docking/ACT/yixianji/data_output/1_PDBQT/protein/file_1"

# Acetyl-CoA PDB
acetyl_coa_pdb = "/public/agis/yanjianbin_group/liuyugeng/taxus_alphaflow_docking/ACT/yixianji/data_output/acetyl-CoA.pdb"

# User parameters
affinity_cut = -3.0
CLASH_CUTOFF = 2.0
max_d_Nu_S1P = 6.0
max_d_Nu_C9 = 5.0
BD_target = 105.0
BD_window = 20.0
max_dist_NE2_Nu = 5.0
ALLOWED_NU_ELEMS = {"O", "N", "S"}

# CSV output
CSV_DIR = "/public/agis/yanjianbin_group/liuyugeng/taxus_alphaflow_docking/ACT/yixianji/data_output/6_CSV_OUTPUT/gating_csv"
os.makedirs(CSV_DIR, exist_ok=True)

CSV_HEADER = [
    "Ligand_id",
    "protein_id",
    "conformation",
    "mode",
    "affinity_kcal",
    "nu_type",
    "nu_atom_id",
    "dist_NE2_Nu",
    "dist_Nu_S1P",
    "dist_Nu_C9",
    "BD_angle_deg",
    "sele_name",
    "clash_flag"
]

# -------------------------
# Argparse: optional --file-range
# -------------------------
def parse_file_range_arg(arg_str):
    """
    Accepts:
      - single integer: "7" -> [7]
      - inclusive range: "1-100" -> [1,2,...,100]
    """
    if arg_str is None:
        return None
    s = arg_str.strip()
    if re.fullmatch(r"\d+", s):
        return [int(s)]
    m = re.fullmatch(r"(\d+)-(\d+)", s)
    if m:
        start = int(m.group(1))
        end = int(m.group(2))
        if start > end:
            raise ValueError(f"--file-range start ({start}) is greater than end ({end}).")
        return list(range(start, end + 1))
    raise ValueError(f"Invalid --file-range format: '{arg_str}'. Expected 'N' or 'N-M'.")

def discover_all_file_nums(base_dir):
    """
    Find all subfolders named 'file_<num>' and return the list of <num> strings, sorted numerically.
    """
    subdirs = [d for d in glob.glob(os.path.join(base_dir, "file_*")) if os.path.isdir(d)]
    nums = []
    for d in subdirs:
        base = os.path.basename(d)
        m = re.fullmatch(r"file_(\d+)", base)
        if m:
            nums.append(m.group(1))
    return sorted(nums, key=lambda x: int(x))

parser = argparse.ArgumentParser(description="ACT gating script with optional subfolder range control.")
parser.add_argument(
    "--file-range",
    type=str,
    default=None,
    help="Restrict processed subfolders under docking_base_dir. Examples: '7' or '1-100'. "
         "Targets directories named file_<num>."
)
args = parser.parse_args()

# -------------------------
# Resolve target file_nums
# -------------------------
if args.file_range is None:
    file_nums = discover_all_file_nums(docking_base_dir)
else:
    nums = parse_file_range_arg(args.file_range)
    file_nums = [str(n) for n in nums]

if not file_nums:
    print(f"No target subfolders found under {docking_base_dir}. Exiting.")
    exit()
print(f"Resolved file_nums: {file_nums}")

# -------------------------
# Discover protein files
# -------------------------
protein_files = glob.glob(os.path.join(protein_dir, "*.pdbqt"))
if not protein_files:
    print(f"No .pdbqt files found in {protein_dir}. Exiting.")
    exit()
print(f"Found {len(protein_files)} protein files to process.")

# -------------------------
# Initialize PyMOL (headless)
# -------------------------
pymol.finish_launching(['pymol', '-c'])

# -------------------------
# Geometry helpers
# -------------------------
def distance(coord1, coord2):
    return np.linalg.norm(np.array(coord1) - np.array(coord2))

def angle_between_points(a, b, c):
    ba = np.array(a) - np.array(b)
    bc = np.array(c) - np.array(b)
    na = np.linalg.norm(ba)
    nc = np.linalg.norm(bc)
    if na < 1e-6 or nc < 1e-6:
        return float('nan')
    cosine_angle = np.clip(np.dot(ba, bc) / (na * nc), -1.0, 1.0)
    angle = np.arccos(cosine_angle)
    return np.degrees(angle)

def has_clash(pose_obj: str, cutoff: float) -> bool:
    """
    Clash check aligned with UGT/CYP style: heavy-atom proximity within cutoff.
    """
    if pose_obj not in cmd.get_names('objects'):
        return False
    try:
        if cmd.count_atoms(pose_obj) == 0:
            return False
    except pymol.CmdException:
        return False
    sel = f"protein1 and not hydro and (({pose_obj} and not hydro) around {cutoff})"
    try:
        return cmd.count_atoms(sel) > 0
    except pymol.CmdException:
        return False

def find_strict_nucleophiles(lig_obj):
    """
    Return list of dicts: [{"id": int, "coord": [x,y,z], "elem": "O/N/S"}]
    Requires X–H explicitly present (X in ALLOWED_NU_ELEMS).
    """
    nu_atoms = []
    for X in ALLOWED_NU_ELEMS:
        sel = f"({lig_obj} and elem {X}) and (neighbor ({lig_obj} and elem H))"
        try:
            if cmd.count_atoms(sel) == 0:
                continue
            m = cmd.get_model(sel, 1)
            for a in m.atom:
                nu_atoms.append({
                    "id": int(a.id),
                    "coord": [a.coord[0], a.coord[1], a.coord[2]],
                    "elem": X
                })
        except pymol.CmdException:
            continue
    return nu_atoms

# -------------------------
# Main loop
# -------------------------
last_csv_path = None

for file_num in file_nums:
    print(f"\nProcessing file_num: {file_num}")

    csv_path_this = os.path.join(CSV_DIR, f"file_{file_num}.csv")
    with open(csv_path_this, "w", newline="") as csv_fh:
        csv_writer = csv.writer(csv_fh)
        csv_writer.writerow(CSV_HEADER)

        docking_dir = os.path.join(docking_base_dir, f"file_{file_num}")
        output_dir  = f"/public/agis/yanjianbin_group/liuyugeng/taxus_alphaflow_docking/ACT/yixianji/data_output/6_conformational_analysis/file_{file_num}"
        os.makedirs(output_dir, exist_ok=True)

        for protein_file in protein_files:
            protein1 = os.path.basename(protein_file).replace(".pdbqt", "")
            print(f"Processing {protein1}...")

            try:
                result_file = os.path.join(docking_dir, f"{file_num}_{protein1}_cavity_1.out")
                ligand_file = os.path.join(docking_dir, f"{file_num}@{protein1}.pdbqt")

                if not os.path.exists(result_file) or not os.path.exists(ligand_file):
                    print(f"Skip {protein1}: missing result or ligand file.")
                    continue

                list1 = []
                mode2aff = {}
                with open(result_file, "r") as f:
                    lines = f.readlines()
                for line in lines[2:]:
                    parts = line.split()
                    if len(parts) >= 2:
                        try:
                            mode, affinity = int(parts[0]), float(parts[1])
                            mode2aff[mode] = affinity
                            if affinity < affinity_cut:
                                list1.append(mode)
                        except ValueError:
                            continue
                if not list1:
                    print(f"Skip {protein1}: no mode with affinity < {affinity_cut}.")
                    continue

                cmd.reinitialize()
                cmd.load(protein_file, "protein1")
                cmd.load(ligand_file, "ligand")
                cmd.load(acetyl_coa_pdb, "coa")

                residue_type = None
                C_name = None
                O_name = None
                S_name = "S1P"

                cmd.select("C9_H5L", "coa and resn H5L and name C9")
                if cmd.count_atoms("C9_H5L") > 0:
                    residue_type = "H5L"
                    C_name = "C9"
                    O_name = "O10"
                else:
                    cmd.select("CM1_MLC", "coa and resn MLC and name CM1")
                    if cmd.count_atoms("CM1_MLC") > 0:
                        residue_type = "MLC"
                        C_name = "CM1"
                        O_name = "OM2"

                if residue_type is None:
                    print(f"Skip {protein1}: neither H5L(C9) nor MLC(CM1) found.")
                    continue

                zb1 = cmd.get_coords(f"coa and resn {residue_type} and name {C_name}", 1)[0]
                cmd.select("O_atom", f"coa and resn {residue_type} and name {O_name}")
                if cmd.count_atoms("O_atom") == 0:
                    print(f"Skip {protein1}: {residue_type} {O_name} not found.")
                    continue
                zb4 = cmd.get_coords("O_atom", 1)[0]

                cmd.select("S1P_atom", f"coa and resn {residue_type} and name {S_name}")
                if cmd.count_atoms("S1P_atom") == 0:
                    print(f"Skip {protein1}: {residue_type} S1P not found.")
                    continue
                S1P_coord = cmd.get_coords("S1P_atom", 1)[0]

                cmd.select("near_His", f"protein1 and resn HIS within 10 of (coa and resn {residue_type} and name {C_name})")
                if cmd.count_atoms("near_His") == 0:
                    print(f"Skip {protein1}: no His within 10 Å of anchor carbon.")
                    continue

                cmd.select("NE2_atoms", "near_His and name NE2")
                zb2_list = cmd.get_coords("NE2_atoms", 1)
                if zb2_list is None or len(zb2_list) == 0:
                    print(f"Skip {protein1}: NE2 coordinates not found.")
                    continue

                cmd.split_states("ligand")
                total_states = cmd.count_states("ligand")
                model_names = [f"ligand_{str(i).zfill(4)}" for i in range(1, total_states + 1)]
                for i in range(1, total_states + 1):
                    if i not in list1:
                        cmd.delete(model_names[i - 1])
                retained_models = [f"ligand_{str(mode).zfill(4)}" for mode in list1]

                results = {}
                for model in retained_models[:]:
                    if model not in cmd.get_names('objects'):
                        continue

                    if has_clash(model, CLASH_CUTOFF):
                        cmd.delete(model)
                        continue

                    nu_list = find_strict_nucleophiles(model)
                    if not nu_list:
                        cmd.delete(model)
                        continue

                    kept = False
                    hit_idx = 0
                    mode_num = int(model.split('_')[1])
                    aff = mode2aff.get(mode_num, "")

                    for nu in nu_list:
                        Nu_xyz = nu["coord"]
                        X = nu["elem"]
                        nu_id = nu["id"]

                        for zb2 in zb2_list:
                            dist_NE2_Nu = distance(zb2, Nu_xyz)
                            if dist_NE2_Nu < max_dist_NE2_Nu:
                                dist_Nu_C9  = distance(zb1, Nu_xyz)
                                dist_Nu_S1P = distance(S1P_coord, Nu_xyz)
                                if dist_Nu_C9 < max_d_Nu_C9 and dist_Nu_S1P < max_d_Nu_S1P:
                                    ang = angle_between_points(Nu_xyz, zb1, zb4)
                                    if not np.isnan(ang) and abs(ang - BD_target) <= BD_window:
                                        sel_name = f"nu_atom_{mode_num:04d}_{hit_idx:02d}"
                                        try:
                                            cmd.select(sel_name, f"{model} and id {nu_id}")
                                        except:
                                            pass
                                        hit_idx += 1

                                        if model not in results:
                                            results[model] = []
                                        results[model].append({
                                            "nu_type": X,
                                            "nu_id": nu_id,
                                            "dist_NE2_Nu": dist_NE2_Nu,
                                            "dist_Nu_S1P": dist_Nu_S1P,
                                            "dist_Nu_C9": dist_Nu_C9,
                                            "angle": ang,
                                            "sele_name": sel_name,
                                            "affinity": aff
                                        })

                                        protein_id = protein1.split('_')[0]
                                        csv_writer.writerow([
                                            file_num,
                                            protein_id,
                                            protein1,
                                            mode_num,
                                            float(aff) if aff != "" else "",
                                            X,
                                            int(nu_id),
                                            float(dist_NE2_Nu),
                                            float(dist_Nu_S1P),
                                            float(dist_Nu_C9),
                                            float(ang),
                                            sel_name,
                                            False
                                        ])
                                        kept = True
                                        # keep collecting multiple hits within the same pose

                    if not kept:
                        cmd.delete(model)

                if results:
                    pse_file = os.path.join(output_dir, f"{protein1}.pse")
                    cmd.save(pse_file)

                    txt_file = os.path.join(output_dir, f"{protein1}.txt")
                    with open(txt_file, "w") as f:
                        for model, infos in results.items():
                            mode_num = int(model.split('_')[1])
                            f.write(f"Model: {mode_num}\n")
                            for info in infos:
                                f.write(
                                    f"  nu_type: {info['nu_type']}, "
                                    f"nu_id: {info['nu_id']}, "
                                    f"dist_NE2_Nu: {info['dist_NE2_Nu']:.2f}, "
                                    f"dist_Nu_S1P: {info['dist_Nu_S1P']:.2f}, "
                                    f"dist_Nu_C9: {info['dist_Nu_C9']:.2f}, "
                                    f"angle: {info['angle']:.2f}, "
                                    f"affinity: {info['affinity']}, "
                                    f"sele: {info['sele_name']}\n"
                                )
                    print(f"[Saved] {pse_file} and {txt_file}")
                else:
                    print(f"Skip {protein1}: no poses satisfy distance/angle criteria.")

            except Exception as e:
                print(f"Error while processing {protein1}: {str(e)}")
                continue

    last_csv_path = csv_path_this
    print(f"[CSV written] {csv_path_this}")

# Cleanup
cmd.quit()
print("All protein files processed.")
