#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Unified gating for catalytic triad hydrolases (Ser–His–Asp/Glu):
- Affinity gating + Ser–His geometry check
- Carbonyl-carbon gating around Ser-OG
- Save only passing poses into .pse
- TXT contains passing entries only
- Per-ligand CSV for downstream scoring
- Added: choose a unique Ser/His pair by your rule; record Ser–His distance
- Added: compute BD angle (∠OG–C–O) and record to CSV (not a gate)
- Added: write per-mode ligand atom PyMOL indices into the same CSV

NOTE (2025-08-18): Persist key selections into the .pse without changing any gating logic.
  Persisted selections include:
    - A1_Ser171
    - Ser_FIXED (picked Ser residue)
    - Ser_OG_FIXED (its OG atom)
    - His_FIXED (picked His residue)
    - His_NE2_FIXED (its NE2 atom)
    - TRIAD_FIXED (Ser_FIXED or His_FIXED)
    - Per-carbonyl passing selections like ligand_XXXX_Csel_N (already created by logic)
"""

import os
import re
import csv
import math
import argparse
import pymol
from pymol import cmd

# =========================
# Global config & thresholds
# =========================
PDBQT_DIR = '/public/agis/yanjianbin_group/liuyugeng/taxus_alphaflow_docking/dAC/data_output/1_PDBQT/protein/file_1'
A1_PATH   = '/public/agis/yanjianbin_group/liuyugeng/taxus_alphaflow_docking/dAC/data_output/dAC-T.pdb'
DOCKING_ROOT = '/public/agis/yanjianbin_group/liuyugeng/taxus_alphaflow_docking/dAC/data_output/4_docking_results/'
OUTPUT_DIR   = '/public/agis/yanjianbin_group/liuyugeng/taxus_alphaflow_docking/dAC/data_output/6_conformational_analysis_20250811'
CSV_DIR      = '/public/agis/yanjianbin_group/liuyugeng/taxus_alphaflow_docking/dAC/data_output/6_CSV_OUTPUT/gating_csv'
os.makedirs(CSV_DIR, exist_ok=True)

# ---- thresholds (edit here) ----
AFFINITY_THRESHOLD   = -3.0   # kcal/mol, affinity gating (<= passes)
SER_HIS_DIST_MAX     = 5.0    # Å, Ser OG — His NE2 must be within this
MAX_DIST_TO_SER_OG   = 5.0    # Å, carbonyl carbon to Ser OG
NEAR_O_RADIUS        = 2.0    # Å, to validate “C with two O neighbors”
# -------------------------------

# Append new columns at the end to keep compatibility with existing pipelines
CSV_HEADER = [
    "Ligand_id", "protein_id", "conformation", "mode",
    "affinity_kcal", "carbon_index", "carbon_id", "selection_name", "dist_to_SerOG_A",
    "dist_SerOG_HisNE2_A", "bd_angle_deg", "ligand_atom_ids"
]

# =========================
# Helpers
# =========================
def parse_out_affinity_modes(out_path, thr=AFFINITY_THRESHOLD):
    """Parse Vina .out table; return passing mode list (as str) and mode->affinity dict."""
    modes, mode2aff = [], {}
    if not os.path.exists(out_path):
        return modes, mode2aff
    data_started = False
    with open(out_path, 'r') as f:
        for line in f:
            s = line.rstrip('\n')
            if s.strip().startswith('-----+'):
                data_started = True
                continue
            if not data_started or not s.startswith('   '):
                continue
            m = re.match(r'\s{3}(\d+)\s+([-+]?\d*\.\d+|\d+)', s)
            if m:
                mode = int(m.group(1))
                aff  = float(m.group(2))
                mode2aff[mode] = aff
                if aff <= thr:
                    modes.append(str(mode))
    return modes, mode2aff

def keep_only_gated_ligand_states(list1):
    """Keep ligand_XXXX objects whose mode is in list1; delete others."""
    lig_objs = [o for o in cmd.get_names('objects') if o.startswith('ligand_')]
    for o in lig_objs:
        num = o.replace('ligand_', '').lstrip('0') or '1'
        if num not in list1:
            cmd.delete(o)

def _sel_nonempty(selname, expr):
    """Create selection and return True only if non-empty."""
    cmd.select(selname, expr)
    return cmd.count_atoms(selname) > 0

def _euclid(a, b):
    return math.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2)

def _angle_deg(a, b, c):
    # angle at b: ∠(a-b, c-b)
    bax = a[0]-b[0]; bay = a[1]-b[1]; baz = a[2]-b[2]
    bcx = c[0]-b[0]; bcy = c[1]-b[1]; bcz = c[2]-b[2]
    dot = bax*bcx + bay*bcy + baz*bcz
    na  = math.sqrt(bax*bax + bay*bay + baz*baz) + 1e-12
    nb  = math.sqrt(bcx*bcx + bcy*bcy + bcz*bcz) + 1e-12
    cos = max(-1.0, min(1.0, dot/(na*nb)))
    return math.degrees(math.acos(cos))

def pick_unique_ser_and_his():
    """
    Ser: from prot within 2 Å of A1:Ser171 -> pick the closest to A1:Ser171 CA.
    His: from byres (Ser around 5 Å) -> pick the His with minimal Ser-OG..His-NE2 distance (<= SER_HIS_DIST_MAX).
    Return (ser_chain, ser_resi, ser_og_xyz, his_chain, his_resi, dist_ser_his) or None.
    """
    if not _sel_nonempty('ser171', 'A1 and resn SER and resi 171'):
        return None
    if not _sel_nonempty('ser171_ca', 'ser171 and name CA'):
        return None
    ser171_ca_xyz = cmd.get_coords('ser171_ca', 1)[0]

    if not _sel_nonempty('near_res', 'prot and byres (ser171 around 2)'):
        return None
    if not _sel_nonempty('ser_near', 'near_res and resn SER'):
        return None

    ser_list = []
    cmd.iterate('ser_near and name CA', 'ser_list.append((chain, resi))', space={'ser_list': ser_list})
    if not ser_list:
        return None

    # pick closest Ser to A1:Ser171 CA
    best_ser = None  # (dist, chain, resi)
    for ser_chain, ser_resi in ser_list:
        ser_ca_expr = f"(prot and {'chain '+ser_chain+' and ' if ser_chain and ser_chain.strip() else ''}resi {ser_resi} and resn SER and name CA)"
        if not _sel_nonempty('ser_ca_tmp', ser_ca_expr):
            continue
        d = _euclid(cmd.get_coords('ser_ca_tmp', 1)[0], ser171_ca_xyz)
        if (best_ser is None) or (d < best_ser[0]):
            best_ser = (d, ser_chain, ser_resi)
    if best_ser is None:
        return None

    _, ser_chain, ser_resi = best_ser
    current_ser_expr = f"(prot and {'chain '+ser_chain+' and ' if ser_chain and ser_chain.strip() else ''}resi {ser_resi} and resn SER)"
    if not _sel_nonempty('current_ser', current_ser_expr):
        return None
    if not _sel_nonempty('ser_og', 'current_ser and name OG'):
        return None
    ser_og_xyz = cmd.get_coords('ser_og', 1)[0]

    # pick closest His (NE2) to Ser-OG within 5 Å byres and <= SER_HIS_DIST_MAX
    if not _sel_nonempty('around_ser', f'prot and byres ({current_ser_expr} around 5)'):
        return None
    if not _sel_nonempty('his_near', 'around_ser and resn HIS'):
        return None

    his_list = []
    cmd.iterate('his_near and name CA', 'his_list.append((chain, resi))', space={'his_list': his_list})
    if not his_list:
        return None

    best_his = None  # (dist_ser_his, his_chain, his_resi)
    for his_chain, his_resi in his_list:
        his_ne2_expr = f"(prot and {'chain '+his_chain+' and ' if his_chain and his_chain.strip() else ''}resi {his_resi} and name NE2)"
        if not _sel_nonempty('his_ne2', his_ne2_expr):
            continue
        try:
            d_sh = cmd.get_distance('ser_og', 'his_ne2')
        except:
            continue
        if d_sh <= SER_HIS_DIST_MAX:
            if (best_his is None) or (d_sh < best_his[0]):
                best_his = (d_sh, his_chain, his_resi)
    if best_his is None:
        return None

    dist_ser_his, his_chain, his_resi = best_his
    return (ser_chain, ser_resi, ser_og_xyz, his_chain, his_resi, dist_ser_his)

# =========================
# Main processing
# =========================
def process_one_protein(protein_name, docking_dir, output_dir, ligand_id, csv_writer):
    cmd.reinitialize()
    try:
        # Use 'prot' (not 'protein') to avoid conflict with PyMOL builtin selector
        cmd.load(os.path.join(PDBQT_DIR, protein_name + '.pdbqt'), 'prot')
        cmd.load(A1_PATH, 'A1')
    except Exception as e:
        print(f'[LoadError] {protein_name}: {e}')
        return

    # Choose unique Ser/His pair by your rule
    picked = pick_unique_ser_and_his()
    if picked is None:
        print(f'[NoSerHisPair] {protein_name}')
        return
    ser_chain, ser_resi, ser_og_xyz, his_chain, his_resi, dist_ser_his = picked

    # ---- Persist key selections for saving in PSE (no logic change) ----
    try:
        # Anchor selection for A1 Ser171
        if cmd.count_atoms('A1 and resn SER and resi 171') > 0:
            cmd.select('A1_Ser171', 'A1 and resn SER and resi 171')
        # Fixed Ser and its OG
        ser_fixed_expr = f"(prot and {'chain '+ser_chain+' and ' if ser_chain and ser_chain.strip() else ''}resi {ser_resi} and resn SER)"
        cmd.select('Ser_FIXED', ser_fixed_expr)
        cmd.select('Ser_OG_FIXED', 'Ser_FIXED and name OG')
        # Fixed His and its NE2
        his_fixed_expr = f"(prot and {'chain '+his_chain+' and ' if his_chain and his_chain.strip() else ''}resi {his_resi} and resn HIS)"
        cmd.select('His_FIXED', his_fixed_expr)
        cmd.select('His_NE2_FIXED', 'His_FIXED and name NE2')
        # Triad convenience selection
        cmd.select('TRIAD_FIXED', 'Ser_FIXED or His_FIXED')
        # (optional purely visual, does not affect logic)
        cmd.show('sticks', 'TRIAD_FIXED or A1_Ser171')
    except Exception as _e:
        # Selections are best-effort; do not interfere with gating
        pass
    # -------------------------------------------------------------------

    ligand_path = os.path.join(docking_dir, f"{ligand_id}@{protein_name}.pdbqt")
    out_path    = os.path.join(docking_dir, f"{ligand_id}_{protein_name}_cavity_1.out")
    list1, mode2aff = parse_out_affinity_modes(out_path, thr=AFFINITY_THRESHOLD)
    if not list1:
        print(f'[NoModesPassThr] {protein_name}')
        return

    try:
        cmd.load(ligand_path, 'ligand')
    except Exception as e:
        print(f'[LoadLigandError] {ligand_path}: {e}')
        return

    # split states -> ligand_0001, ...
    cmd.split_states('ligand')
    keep_only_gated_ligand_states(list1)

    passing_objs = set()
    txt_lines = []
    txt_path = os.path.join(output_dir, f"{protein_name}.txt")
    pse_path = os.path.join(output_dir, f"{protein_name}_filtered.pse")
    os.makedirs(output_dir, exist_ok=True)

    header_written = False

    # scan passing ligand objects only
    objs = [o for o in cmd.get_names('objects') if o.startswith('ligand_')]
    if not objs:
        return

    for obj in objs:
        raw_num = obj.replace('ligand_', '').lstrip('0') or '1'
        aff = mode2aff.get(int(raw_num), None)
        cmd.enable(obj)

        # capture per-mode ligand atom indices (PyMOL index)
        model = cmd.get_model(obj)
        ligand_atom_ids = ",".join(str(a.index) for a in model.atom)

        # Carbonyl carbon: C with O neighbor(s)
        if not _sel_nonempty('carbon_with_o', f"{obj} and elem C and neighbor (elem O)"):
            continue

        local_pass_lines = []
        for i, carbon in enumerate(cmd.get_model('carbon_with_o').atom, 1):
            carbon_index = carbon.index     # 保持原逻辑
            carbon_id    = carbon.id        # 新增：PDBQT 原子序号，可跨构象定位
            sele_name    = f"{obj}_Csel_{i}"
            if not _sel_nonempty(sele_name, f"{obj} and index {carbon_index}"):  # 不改动原逻辑
                continue

            if not _sel_nonempty('near_atoms', f"{obj} within {NEAR_O_RADIUS:.1f} of {sele_name}"):
                cmd.delete(sele_name)
                continue

            o_atoms = [a for a in cmd.get_model('near_atoms').atom if a.symbol == 'O']
            if len(o_atoms) != 2:
                cmd.delete(sele_name)
                continue

            c_xyz = cmd.get_coords(sele_name, 1)[0]
            dist = _euclid(c_xyz, ser_og_xyz)

            if dist <= MAX_DIST_TO_SER_OG:
                # pick nearest O as carbonyl O for BD angle measurement
                o_atom = min(
                    o_atoms,
                    key=lambda a: (a.coord[0]-c_xyz[0])**2 + (a.coord[1]-c_xyz[1])**2 + (a.coord[2]-c_xyz[2])**2
                )
                bd_angle = _angle_deg(ser_og_xyz, c_xyz, o_atom.coord)  # ∠(OG–C–O), degrees

                passing_objs.add(obj)
                local_pass_lines.append(f"{carbon_index}\t{sele_name}\t{dist:.2f}")

                protein_id = protein_name.split('_')[0]
                # write CSV row (append new fields at the end)
                csv_writer.writerow([
                    ligand_id, protein_id, protein_name,
                    int(raw_num), float(aff) if aff is not None else "",
                    int(carbon_index), int(carbon_id), sele_name, float(dist),
                    float(dist_ser_his) if dist_ser_his is not None else "",
                    float(bd_angle) if bd_angle is not None else "",
                    ligand_atom_ids
                ])
            else:
                cmd.delete(sele_name)

        if local_pass_lines:
            if not header_written:
                txt_lines.append(f"Carbonyl-carbon hits within {MAX_DIST_TO_SER_OG:.1f} Å of Ser-OG (fixed Ser/His chosen)\n")
                header_written = True
            txt_lines.append(f"\n--- Object: {obj} ---\n")
            txt_lines.append("index\tselection_name\tdist_to_SerOG(Å)\n")
            txt_lines.extend(line + "\n" for line in local_pass_lines)

    # keep only passing poses, save PSE/TXT for passing entries
    if passing_objs:
        for o in [x for x in cmd.get_names('objects') if x.startswith('ligand_')]:
            if o not in passing_objs:
                cmd.delete(o)

        # Ensure key selections are visible in the saved session (visual only; no logic change)
        try:
            cmd.show('sticks', 'TRIAD_FIXED or A1_Ser171 or (name *_Csel_*)')
        except Exception:
            pass

        try:
            cmd.save(pse_path)
        except Exception as e:
            print(f'[SavePSEError] {protein_name}: {e}')

        try:
            with open(txt_path, 'w') as ftxt:
                ftxt.writelines(txt_lines)
        except Exception as e:
            print(f'[WriteTxtError] {protein_name}: {e}')
    else:
        print(f"[NoPassingPose] {protein_name}")

# =========================
# Entry
# =========================
def _parse_range_arg(rstr):
    m = re.match(r'^\s*(\d+)\s*-\s*(\d+)\s*$', rstr)
    if not m:
        raise ValueError("Range must be in the form START-END (e.g., 1-100)")
    a, b = int(m.group(1)), int(m.group(2))
    if a > b:
        a, b = b, a
    return a, b

def main():
    parser = argparse.ArgumentParser(description="Unified gating runner")
    parser.add_argument(
        "--range",
        dest="range",
        type=str,
        default=None,
        help="Limit processed subfolders to an inclusive range, e.g. '1-100' for file_1..file_100"
    )
    args = parser.parse_args()

    pymol.finish_launching(['pymol', '-cq'])

    subdirs = [d for d in os.listdir(DOCKING_ROOT) if d.startswith('file_') and os.path.isdir(os.path.join(DOCKING_ROOT, d))]
    subdirs.sort(key=lambda x: int(x.split('_')[-1]))

    # Apply optional range filter
    if args.range:
        lo, hi = _parse_range_arg(args.range)
        subdirs = [d for d in subdirs if lo <= int(d.split('_')[-1]) <= hi]

    for sub in subdirs:
        ligand_id = sub.split('_')[-1]
        docking_dir = os.path.join(DOCKING_ROOT, sub)
        out_lig_dir = os.path.join(OUTPUT_DIR, sub)
        os.makedirs(out_lig_dir, exist_ok=True)

        csv_path = os.path.join(CSV_DIR, f"file_{ligand_id}.csv")
        with open(csv_path, 'w', newline='') as fh:
            writer = csv.writer(fh)
            writer.writerow(CSV_HEADER)

            pdbqt_files = [f for f in os.listdir(PDBQT_DIR) if f.endswith('.pdbqt')]
            for pdbqt_file in pdbqt_files:
                protein_name = os.path.splitext(pdbqt_file)[0]
                print(f'[Ligand {ligand_id}] Processing protein: {protein_name}')
                try:
                    process_one_protein(protein_name, docking_dir, out_lig_dir, ligand_id, writer)
                except Exception as e:
                    print(f'[Error] {protein_name}: {e}')
                finally:
                    cmd.delete('all')

        print(f'[CSV written] {csv_path}')

if __name__ == '__main__':
    main()
