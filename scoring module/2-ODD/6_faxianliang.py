#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import math
import glob
import argparse
import pymol2
import numpy as np
import csv

# -------------------------
# Thresholds
# -------------------------
AFFINITY_THRESHOLD = -3.0
# Angle threshold now applies to axis deviation (0..90): pass if axis_dev <= ANGLE_THRESHOLD
ANGLE_THRESHOLD    = 60.0
DIST_MIN           = 2.5
DIST_MAX           = 5.0
CLASH_CUTOFF       = 2.0  # hard clash cutoff (Å)

# -------------------------
# CSV output (one file per ligand)
# -------------------------
CSV_DIR = "/public/agis/yanjianbin_group/liuyugeng/taxus_alphaflow_docking/Fe2OG/data_output/6_CSV_OUTPUT/gating_csv"
os.makedirs(CSV_DIR, exist_ok=True)

# angle_deg now stores axis deviation (0..90) between Fe->H and the plane normal of {Fe,O5,A}
CSV_HEADER = [
    "Ligand_id",            # file_num
    "protein_id",           # prefix before '_' from protein_name
    "conformation",         # protein_name (pdb basename)
    "mode",                 # Vina pose index
    "affinity_kcal",        # from .out
    "dist_FE_H",            # Å (Fe···H)
    "angle_deg",            # axis deviation (deg): min(angle(n̂, Fe->H), 180-angle)
    "lig_atom_symbol",      # 'H'
    "lig_atom_name",        # e.g., H22
    "lig_atom_id",          # PyMOL atom id
    "parent_atom_element",  # 'C'
    "parent_atom_name",     # e.g., C12
    "parent_atom_id",       # PyMOL atom id of C
    "fe_atom_id",           # Fe atom id
    "o_ref_atom_id",        # AKG O5 atom id (reference atom ID stored for bookkeeping)
    "clash_flag"            # False (only rows that passed gating are written)
]

# -------------------------
# Vina .out parsing
# -------------------------
def parse_docking_out(out_file, affinity_threshold=AFFINITY_THRESHOLD):
    """
    Parse a Vina .out file.

    Returns:
      list1: list of modes whose affinity < threshold (stops at first failure)
      mode2aff: dict mapping mode -> affinity
    """
    list1 = []
    mode2aff = {}
    if not os.path.exists(out_file):
        return list1, mode2aff
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

# -------------------------
# Atom locators
# -------------------------
def get_atom_once(cmd, sel):
    """
    Return (coord_tuple, atom_id) for a selection expected to have at least one atom; otherwise (None, None).
    """
    model = cmd.get_model(sel)
    if not model.atom:
        return None, None
    a = model.atom[0]
    return a.coord, a.id

def get_fe_coords(cmd):
    return get_atom_once(cmd, "protein and name FE")

def get_o5_coords(cmd):
    return get_atom_once(cmd, "protein and resn AKG and name O5")

def get_anchorA_coords(cmd):
    """
    Prefer O1, then fallback to O2 -> C1 -> C2 (all within resn AKG).
    Returns (coord, atom_id, name_tried) or (None, None, None).
    """
    candidates = [
        ("protein and resn AKG and name O1", "O1"),
        ("protein and resn AKG and name O2", "O2"),
        ("protein and resn AKG and name C1", "C1"),
        ("protein and resn AKG and name C2", "C2"),
    ]
    for sel, tag in candidates:
        coord, aid = get_atom_once(cmd, sel)
        if coord:
            return coord, aid, tag
    return None, None, None

# -------------------------
# Geometry helpers
# -------------------------
def v_sub(a, b):
    return [a[0]-b[0], a[1]-b[1], a[2]-b[2]]

def v_len(v):
    return math.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])

def v_dot(a, b):
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

def v_cross(a, b):
    return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]

def angle_between(v1, v2):
    """
    Return angle between v1 and v2 in degrees, range 0..180.
    """
    n1 = v_len(v1)
    n2 = v_len(v2)
    if n1 == 0 or n2 == 0:
        return None
    dot = v_dot(v1, v2) / (n1 * n2)
    dot = max(min(dot, 1.0), -1.0)
    return math.degrees(math.acos(dot))

def axis_deviation(fe, o5, a, h):
    """
    Compute axis deviation (0..90) between plane normal n̂ (from Fe->O5 and Fe->A) and Fe->H.
    Returns None if degenerate.
    """
    vO5 = v_sub(o5, fe)
    vA  = v_sub(a, fe)
    n   = v_cross(vO5, vA)
    nlen = v_len(n)
    if nlen < 1e-8:
        return None
    n = [n[0]/nlen, n[1]/nlen, n[2]/nlen]
    vH = v_sub(h, fe)
    theta = angle_between(n, vH)  # 0..180
    if theta is None:
        return None
    return min(theta, 180.0 - theta)  # 0..90

def find_parent_carbon(ligand_model, h_atom, cutoff=1.25):
    """
    Find the nearest carbon to a given hydrogen within cutoff (Å) in the ligand model.
    Returns (C_atom, dist) if found, otherwise (None, None).
    """
    hx, hy, hz = h_atom.coord
    best = (None, None)
    dmin2 = cutoff * cutoff
    for a in ligand_model.atom:
        if a.symbol != 'C':
            continue
        dx = a.coord[0] - hx
        dy = a.coord[1] - hy
        dz = a.coord[2] - hz
        d2 = dx*dx + dy*dy + dz*dz
        if d2 <= dmin2:
            dmin2 = d2
            best = (a, math.sqrt(d2))
    return best

# -------------------------
# Main per-protein worker
# -------------------------
def process_protein(cmd, protein_name, protein_dir, docking_dir, output_dir, file_num,
                    csv_writer,
                    angle_threshold=ANGLE_THRESHOLD, affinity_threshold=AFFINITY_THRESHOLD,
                    dist_min=DIST_MIN, dist_max=DIST_MAX, clash_cutoff=CLASH_CUTOFF):
    """
    Apply gating for a given protein and write per-passing-H rows to CSV.
    Angle gating uses axis deviation relative to the equatorial/chelation-plane normal.
    """
    out_file = os.path.join(docking_dir, f'{file_num}_' + protein_name + '_cavity_1.out')
    if not os.path.exists(out_file):
        print("Docking .out file not found:", out_file)
        return

    list1, mode2aff = parse_docking_out(out_file, affinity_threshold=affinity_threshold)
    if not list1:
        print(f"No modes passed affinity threshold for protein {protein_name}")
        return

    # load protein
    protein_path = os.path.join(protein_dir, protein_name + '.pdb')
    if not os.path.exists(protein_path):
        print("Protein PDB not found:", protein_path)
        return
    cmd.load(protein_path, "protein")

    # load ligand
    ligand_file = f'{file_num}@' + protein_name + '.pdbqt'
    ligand_path = os.path.join(docking_dir, ligand_file)
    if not os.path.exists(ligand_path):
        print("Ligand file not found:", ligand_path)
        cmd.delete("all")
        return
    cmd.load(ligand_path, "ligand")

    # split ligand states to objects named "1","2",...
    cmd.split_states("ligand", prefix="object_ligand")
    objs = cmd.get_object_list("object_ligand*")
    objs.sort(key=lambda x: int(x.replace('object_ligand', '')))
    for i, obj in enumerate(objs, start=1):
        cmd.set_name(obj, str(i))

    # Fetch Fe, O5, and anchor A (O1 preferred; fallback O2 -> C1 -> C2)
    fe_coords, fe_atom_id = get_fe_coords(cmd)
    if not fe_coords:
        cmd.delete("all")
        return
    o5_coords, o5_atom_id = get_o5_coords(cmd)
    if not o5_coords:
        cmd.delete("all")
        return
    a_coords, a_atom_id, a_tag = get_anchorA_coords(cmd)
    if not a_coords:
        print("Failed to locate AKG anchor atom (O1/O2/C1/C2).")
        cmd.delete("all")
        return

    # Precompute plane normal quality (for logging)
    _axis_dev_probe = axis_deviation(fe_coords, o5_coords, a_coords, (fe_coords[0]+1.0, fe_coords[1], fe_coords[2]))
    # Not used in gating; ensures the normal is not degenerate. If None, we'll try to fail fast later.

    valid_ligand_objects = []
    protein_angles_info = []

    for mode in list1:
        ligand_obj_name = str(mode)
        if ligand_obj_name not in cmd.get_names("all"):
            print("Pose object missing:", ligand_obj_name)
            continue

        # add hydrogens for this pose
        cmd.h_add(ligand_obj_name)
        ligand_model = cmd.get_model(ligand_obj_name)

        within_range = False
        candidate_info = []  # (dist, axis_dev_deg, H_name, H_id, C_name, C_id)

        # iterate over H atoms only
        for atom in ligand_model.atom:
            if atom.symbol != 'H':
                continue

            # keep only C–H hydrogens (nearest C within 1.25 Å)
            c_atom, ch_dist = find_parent_carbon(ligand_model, atom, cutoff=1.25)
            if c_atom is None:
                continue  # skip O–H / N–H / others

            # distance Fe···H
            x, y, z = atom.coord
            dx = fe_coords[0] - x
            dy = fe_coords[1] - y
            dz = fe_coords[2] - z
            dist_fe_h = math.sqrt(dx*dx + dy*dy + dz*dz)

            if dist_min <= dist_fe_h <= dist_max:
                # axis deviation against plane normal n̂(Fe,O5,A)
                axis_dev_deg = axis_deviation(fe_coords, o5_coords, a_coords, (x, y, z))
                if axis_dev_deg is None:
                    continue

                if axis_dev_deg <= angle_threshold:
                    candidate_info.append((
                        dist_fe_h,
                        axis_dev_deg,
                        atom.name, atom.id,
                        c_atom.name, c_atom.id
                    ))
                    within_range = True

        # hard clash at 2.0 Å (protein heavy atoms)
        if within_range:
            clash_sel = f"protein and not hydro and (({ligand_obj_name} and not hydro) around {clash_cutoff})"
            clash_count = cmd.count_atoms(clash_sel)
            if clash_count > 0:
                within_range = False
                print(f"Pose {ligand_obj_name} has clashes (cutoff={clash_cutoff} Å, count={clash_count}); marked invalid.")

        # for angles.txt
        if within_range and candidate_info:
            protein_angles_info.append((ligand_obj_name, candidate_info))

        # CSV rows per passing H atom
        if within_range and candidate_info:
            parsed_protein_id = protein_name.split('_')[0]
            conformation = protein_name
            aff = mode2aff.get(mode, float('nan'))
            for (dist, axis_dev_deg, hname, hid, cname, cid) in candidate_info:
                csv_writer.writerow([
                    file_num,
                    parsed_protein_id,
                    conformation,
                    mode,
                    float(aff),
                    float(dist),
                    float(axis_dev_deg),
                    "H",
                    hname,
                    int(hid),
                    "C",
                    cname,
                    int(cid),
                    int(fe_atom_id),
                    int(o5_atom_id),
                    False
                ])

        if within_range:
            valid_ligand_objects.append(ligand_obj_name)
            # visualization: select qualified hydrogens
            qualifying_atoms = [hid for (_, _, _, hid, _, _) in candidate_info]
            if qualifying_atoms:
                sel_str = " or ".join([f"id {atom_id}" for atom_id in qualifying_atoms])
                cmd.select(f"candidates_{ligand_obj_name}", f"{ligand_obj_name} and ({sel_str})")
        else:
            cmd.delete(ligand_obj_name)

    # write angles.txt (describe the new angle definition)
    angles_filename = os.path.join(output_dir, protein_name + "_angles.txt")
    try:
        with open(angles_filename, "w") as f_angles:
            f_angles.write("Parameters:\n")
            f_angles.write(f"Affinity Threshold: {affinity_threshold}\n")
            f_angles.write(f"Angle Threshold (axis deviation): {angle_threshold} deg (pass if <= threshold)\n")
            f_angles.write("Axis definition: normal to the equatorial/chelation plane spanned by Fe->O5 and Fe->A, A in {O1,O2,C1,C2}\n")
            f_angles.write(f"Distance Range (Fe···H): {dist_min} - {dist_max} Å\n")
            f_angles.write(f"Clash Cutoff (hard): {clash_cutoff} Å\n\n")
            for ligand_obj, info in protein_angles_info:
                f_angles.write(f"Ligand conformation: {ligand_obj}\n")
                for entry in info:
                    dist, axis_dev_deg, hname, hid, cname, cid = entry
                    f_angles.write(
                        " Fe-H: {:.2f} Å, Axis deviation: {:.2f} deg, H(name:id) {}:{}, C(name:id) {}:{}\n"
                        .format(dist, axis_dev_deg, hname, hid, cname, cid)
                    )
                f_angles.write("\n")
        print("Fe-Ligand angle info saved to file:", angles_filename)
    except Exception as e:
        print("Failed to write angle info file:", e)

    if valid_ligand_objects:
        out_session = os.path.join(output_dir, protein_name + ".pse")
        cmd.save(out_session)
        print(f"Saved session for protein {protein_name}: {out_session}")
    else:
        print(f"All ligand poses failed gating for protein {protein_name}")

# -------------------------
# CLI helpers for folder selection
# -------------------------

DEFAULT_DOCKING_BASE_DIR = "/public/agis/yanjianbin_group/liuyugeng/taxus_alphaflow_docking/Fe2OG/data_output/4_docking_results"

def parse_file_range(spec):
    """
    Parse a range spec like '1-100' or a single number '7'.
    Returns a sorted list of strings: ['1','2',...].
    """
    if spec is None:
        return None
    spec = spec.strip()
    m = re.fullmatch(r"(\d+)\s*-\s*(\d+)", spec)
    if m:
        a, b = int(m.group(1)), int(m.group(2))
        if a > b:
            a, b = b, a
        return [str(i) for i in range(a, b + 1)]
    m2 = re.fullmatch(r"\d+", spec)
    if m2:
        return [spec]
    raise ValueError(f"Invalid --file-range spec: {spec}. Use 'START-END' or a single integer, e.g., '1-100' or '7'.")

def discover_file_nums(docking_base_dir, file_range=None):
    """
    Determine which file_* folders to process.
    - If file_range is None: scan all file_* subfolders.
    - If file_range is a list of numbers as strings: use those (filtering to existing folders).
    Returns a sorted list of numbers as strings.
    """
    if file_range is None:
        docking_dirs = glob.glob(os.path.join(docking_base_dir, "file_*"))
        file_nums = sorted(
            [os.path.basename(f).split('_')[-1]
             for f in docking_dirs
             if os.path.isdir(f) and '_' in os.path.basename(f)],
            key=int
        )
        return file_nums
    else:
        file_nums = []
        for n in file_range:
            path = os.path.join(docking_base_dir, f"file_{n}")
            if os.path.isdir(path):
                file_nums.append(str(n))
            else:
                print(f"Warning: directory not found, skipping {path}")
        file_nums = sorted(set(file_nums), key=int)
        return file_nums

def main(cmd, docking_base_dir, file_nums):
    protein_dir = '/public/agis/yanjianbin_group/liuyugeng/taxus_alphaflow_docking/Fe2OG/data_input/data_protein/partition_1'
    output_base_dir  = '/public/agis/yanjianbin_group/liuyugeng/taxus_alphaflow_docking/Fe2OG/data_output/6_conformational_analysis_20250811'

    if not file_nums:
        print(f"No eligible folders found in {docking_base_dir}. Exiting.")
        exit()

    print(f"Discovered file_nums: {file_nums}")

    for file_num in file_nums:
        print(f"\nProcessing file_num: {file_num}")
        docking_dir = os.path.join(docking_base_dir, f"file_{file_num}")
        output_dir  = os.path.join(output_base_dir, f"file_{file_num}")
        os.makedirs(output_dir, exist_ok=True)

        # per-ligand CSV
        csv_path_this = os.path.join(CSV_DIR, f"file_{file_num}.csv")
        with open(csv_path_this, "w", newline="") as csv_fh:
            writer = csv.writer(csv_fh)
            writer.writerow(CSV_HEADER)

            for file in os.listdir(protein_dir):
                if file.endswith('.pdb'):
                    protein_name = file[:-4]
                    print("Processing protein:", protein_name)
                    process_protein(cmd, protein_name, protein_dir, docking_dir, output_dir, file_num, writer)
                    cmd.delete("all")

        print(f"Wrote CSV: {csv_path_this}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Fe/2OG gating pipeline (geometry-based prefilter).")
    parser.add_argument(
        "--docking-base-dir",
        default=DEFAULT_DOCKING_BASE_DIR,
        help=f"Base directory containing file_* subfolders (default: {DEFAULT_DOCKING_BASE_DIR})"
    )
    parser.add_argument(
        "--file-range",
        default=None,
        help="Folder number range to process, e.g., '1-100' or a single number '7'. "
             "If omitted, all file_* subfolders under --docking-base-dir are processed."
    )
    args = parser.parse_args()

    try:
        sel = parse_file_range(args.file_range)
    except ValueError as ve:
        print(str(ve))
        exit(2)

    file_nums = discover_file_nums(args.docking_base_dir, sel)

    with pymol2.PyMOL() as pymol:
        cmd = pymol.cmd
        main(cmd, args.docking_base_dir, file_nums)
