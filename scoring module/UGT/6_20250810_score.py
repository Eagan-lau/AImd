import os
import glob
import csv
import argparse
import re
import pymol
from pymol import cmd, stored
import numpy as np

# =========================
# Config & thresholds
# =========================
AFFINITY_THRESHOLD = -3.0      # kcal/mol (coarse screen; tighten to -5.5 if needed)
NE2_O_DISTANCE_THRESHOLD = 5.0 # Å
C1_O_DISTANCE_THRESHOLD  = 5.0 # Å
CLASH_CUTOFF = 2.0             # Å (heavy-atom clash cutoff)
ANGLE_MIN = 100.0              # °
ANGLE_MAX = 180.1              # °

# Whether to save .pse (slow and large; default True to match original behavior)
SAVE_PSE = True

# Input/output directories
docking_base_dir = "/public/agis/yanjianbin_group/liuyugeng/taxus_alphaflow_docking/UDPGT/data_output/4_docking_results"
protein_dir      = "/public/agis/yanjianbin_group/liuyugeng/taxus_alphaflow_docking/UDPGT/data_output/1_PDBQT/protein/file_1"
output_base_dir  = "/public/agis/yanjianbin_group/liuyugeng/taxus_alphaflow_docking/UDPGT/data_output/6_conformational_analysis_20250810"

# Per-ligand CSV output directory
csv_dir = "/public/agis/yanjianbin_group/liuyugeng/taxus_alphaflow_docking/UDPGT/data_output/6_CSV_OUTPUT/gating_csv"
os.makedirs(csv_dir, exist_ok=True)

# CSV header (c4_atom_id -> o1_atom_id; protein_file -> conformation)
CSV_HEADER = [
    "Ligand_id",
    "protein_id",
    "conformation",
    "mode",
    "affinity_kcal",
    "dist_NE2_O",
    "dist_C1_O",
    "angle_deg",
    "ne2_atom_id",
    "o_atom_id",
    "c1_atom_id",
    "o1_atom_id",
    "clash_flag"
]

# =========================
# Helpers
# =========================
def distance(coord1, coord2):
    return np.linalg.norm(np.array(coord1) - np.array(coord2))

def calculate_angle(vec1, vec2):
    dot_product = np.dot(vec1, vec2)
    norm_vec1 = np.linalg.norm(vec1)
    norm_vec2 = np.linalg.norm(vec2)
    if norm_vec1 == 0 or norm_vec2 == 0:
        return None
    cos_angle = np.clip(dot_product / (norm_vec1 * norm_vec2), -1.0, 1.0)
    return np.degrees(np.arccos(cos_angle))

def detect_clashes(protein_sel, ligand_sel, cutoff):
    protein_sel_heavy = f"({protein_sel}) and not hydro"
    ligand_sel_heavy  = f"({ligand_sel}) and not hydro"
    clashes = cmd.find_pairs(protein_sel_heavy, ligand_sel_heavy, cutoff=cutoff, mode=0)
    return len(clashes) > 0

def parse_file_range_arg(arg_str):
    """
    Supports:
      - single integer: "7" -> [7]
      - inclusive range: "1-100" -> [1, 2, ..., 100]
    """
    if arg_str is None:
        return None
    s = arg_str.strip()
    if re.fullmatch(r"\d+", s):
        return [int(s)]
    m = re.fullmatch(r"(\d+)\-(\d+)", s)
    if m:
        start = int(m.group(1))
        end = int(m.group(2))
        if start > end:
            raise ValueError(f"--file-range start ({start}) is greater than end ({end}).")
        return list(range(start, end + 1))
    raise ValueError(f"Invalid --file-range format: '{arg_str}'. Expected 'N' or 'N-M'.")

def discover_all_file_nums(base_dir):
    """
    Finds all subfolders named 'file_<num>' directly under base_dir and returns the list of <num> (as strings) sorted.
    """
    subdirs = [d for d in glob.glob(os.path.join(base_dir, "file_*")) if os.path.isdir(d)]
    nums = []
    for d in subdirs:
        base = os.path.basename(d)
        m = re.fullmatch(r"file_(\d+)", base)
        if m:
            nums.append(m.group(1))
    return sorted(nums, key=lambda x: int(x))

# =========================
# Argparse
# =========================
parser = argparse.ArgumentParser(description="UGT gating and geometry screen with optional subfolder range control.")
parser.add_argument("--file-range", type=str, default=None,
                    help="Optional range to restrict processed subfolders. Examples: '7' or '1-100'. "
                         "This will process directories named file_<num> under docking_base_dir.")
args = parser.parse_args()

# =========================
# Prep file lists
# =========================
if args.file_range is None:
    # Default: process all subfolders under docking_base_dir
    file_nums = discover_all_file_nums(docking_base_dir)
else:
    # Use provided numeric range
    nums = parse_file_range_arg(args.file_range)
    # Keep as strings for later formatting
    file_nums = [str(n) for n in nums]

if not file_nums:
    print(f"No target subfolders resolved under {docking_base_dir}. Exiting.")
    raise SystemExit

print(f"Target file_nums: {file_nums}")

protein_files = glob.glob(os.path.join(protein_dir, "*.pdbqt"))
if not protein_files:
    print(f"No .pdbqt files found in {protein_dir}. Exiting.")
    raise SystemExit

print(f"Found {len(protein_files)} protein files to process.")

# =========================
# Init PyMOL
# =========================
pymol.finish_launching(['pymol', '-c'])

try:
    # =========================
    # Main loop
    # =========================
    for file_num in file_nums:
        print(f"\nProcessing file_num: {file_num}")
        docking_dir = os.path.join(docking_base_dir, f"file_{file_num}")
        if not os.path.isdir(docking_dir):
            print(f"Skipping file_{file_num}: directory does not exist under docking_base_dir.")
            continue

        output_dir  = os.path.join(output_base_dir,  f"file_{file_num}")
        os.makedirs(output_dir, exist_ok=True)

        # Per-ligand CSV (create/overwrite)
        csv_path_this = os.path.join(csv_dir, f"file_{file_num}.csv")
        csv_fh = open(csv_path_this, "w", newline="")
        csv_writer = csv.writer(csv_fh)
        csv_writer.writerow(CSV_HEADER)

        for protein_file in protein_files:
            protein1 = os.path.basename(protein_file).replace(".pdbqt", "")
            print(f"Processing {protein1} with file_num {file_num}...")
            try:
                # Docking outputs
                result_file = os.path.join(docking_dir, f"{file_num}_{protein1}_cavity_1.out")
                ligand_file = os.path.join(docking_dir, f"{file_num}@{protein1}.pdbqt")
                if not os.path.exists(result_file) or not os.path.exists(ligand_file):
                    print(f"Skipping {protein1}: missing result or ligand file.")
                    continue

                # Parse .out to obtain modes and energies; apply energy prefilter
                list1 = []
                affinity_dict = {}
                with open(result_file, "r") as f:
                    lines = f.readlines()
                for line in lines[2:]:  # skip header
                    parts = line.split()
                    if len(parts) >= 2:
                        try:
                            mode, affinity = int(parts[0]), float(parts[1])
                            affinity_dict[mode] = affinity
                            if affinity < AFFINITY_THRESHOLD:
                                list1.append(mode)
                        except ValueError:
                            continue

                if not list1:
                    print(f"Skipping {protein1}: no modes with affinity < {AFFINITY_THRESHOLD}.")
                    continue

                # Load receptor conformation
                cmd.reinitialize()
                cmd.load(protein_file, "protein1")

                # C1
                cmd.select("C1", "protein1 and resn U2F and name C1")
                if cmd.count_atoms("C1") == 0:
                    print(f"Skipping {protein1}: no C1 atom in U2F residue.")
                    continue
                zb1 = cmd.get_coords("C1", 1)[0]
                stored.c1 = {}
                cmd.iterate("C1", "stored.c1['name'] = name; stored.c1['id'] = ID")
                c1_info = f"Atom: {stored.c1['name']} (name: {stored.c1['name']}, id: {stored.c1['id']})"

                # Nearby His -> NE2
                cmd.select("near_His", "protein1 and resn HIS within 10 of C1")
                if cmd.count_atoms("near_His") == 0:
                    print(f"Skipping {protein1}: no His within 10 Å of C1.")
                    continue
                cmd.select("NE2_atoms", "near_His and name NE2")
                zb2_list = cmd.get_coords("NE2_atoms", 1)
                if zb2_list is None or len(zb2_list) == 0:
                    print(f"Skipping {protein1}: failed to get NE2 coordinates.")
                    continue
                stored.ne2_list = []
                cmd.iterate("NE2_atoms", "stored.ne2_list.append({'name': name, 'id': ID})")

                # Ligand poses (keep only energy-passing modes)
                cmd.load(ligand_file, "ligand")
                cmd.split_states("ligand")
                total_states = cmd.count_states("ligand")
                model_names = [f"ligand_{str(i).zfill(4)}" for i in range(1, total_states + 1)]
                for i in range(1, total_states + 1):
                    if i not in list1:
                        cmd.delete(model_names[i - 1])
                retained_models = [f"ligand_{str(mode).zfill(4)}" for mode in list1]

                results = {}
                for model_idx, model in enumerate(retained_models):
                    # Clash filter
                    if detect_clashes("protein1", model, CLASH_CUTOFF):
                        continue

                    # Select O bound to H (e.g., hydroxyl O)
                    cmd.select(f"oh_oxygens_{model}", f"({model} and name O) and (neighbor name H)")
                    zb3_list = cmd.get_coords(f"oh_oxygens_{model}", 1)
                    if zb3_list is None or len(zb3_list) == 0:
                        continue
                    stored.o_list = []
                    cmd.iterate(f"oh_oxygens_{model}", "stored.o_list.append({'name': name, 'id': ID})")

                    for o_idx, zb3 in enumerate(zb3_list):
                        for ne2_idx, zb2 in enumerate(zb2_list):
                            dist_NE2_O = distance(zb2, zb3)
                            if dist_NE2_O < NE2_O_DISTANCE_THRESHOLD:
                                dist_C1_O = distance(zb1, zb3)
                                if dist_C1_O < C1_O_DISTANCE_THRESHOLD:
                                    # Use O1 as the reference atom (instead of C4)
                                    cmd.select("O1", "protein1 and resn U2F and name O1")
                                    if cmd.count_atoms("O1") == 0:
                                        continue
                                    zbO1 = cmd.get_coords("O1", 1)[0]
                                    stored.o1 = {}
                                    cmd.iterate("O1", "stored.o1['name'] = name; stored.o1['id'] = ID")
                                    o1_info = f"Atom: {stored.o1['name']} (name: {stored.o1['name']}, id: {stored.o1['id']})"

                                    # Angle: reference axis C1→O1
                                    vec_C1_O1 = np.array(zbO1) - np.array(zb1)
                                    vec_C1_O  = np.array(zb3)  - np.array(zb1)
                                    angle = calculate_angle(vec_C1_O1, vec_C1_O)
                                    if angle is None or not (ANGLE_MIN < angle < ANGLE_MAX):
                                        continue

                                    if model not in results:
                                        results[model] = []
                                    results[model].append({
                                        "zb3": zb3,
                                        "dist_NE2_O": dist_NE2_O,
                                        "dist_C1_O": dist_C1_O,
                                        "angle": angle,
                                        "c1_info": c1_info,
                                        "o1_info": o1_info,
                                        "ne2_info": f"Atom: {stored.ne2_list[ne2_idx]['name']} (name: {stored.ne2_list[ne2_idx]['name']}, id: {stored.ne2_list[ne2_idx]['id']})",
                                        "o_info": f"Atom: {stored.o_list[o_idx]['name']} (name: {stored.o_list[o_idx]['name']}, id: {stored.o_list[o_idx]['id']})"
                                    })

                                    # Write CSV row
                                    mode_num = int(model.split('_')[1])
                                    affinity = affinity_dict.get(mode_num, np.nan)
                                    parsed_protein_id = protein1.split('_')[0]  # canonical protein ID
                                    conformation = protein1                    # already without .pdbqt

                                    csv_writer.writerow([
                                        file_num,
                                        parsed_protein_id,
                                        conformation,
                                        mode_num,
                                        float(affinity) if affinity is not None else np.nan,
                                        float(dist_NE2_O),
                                        float(dist_C1_O),
                                        float(angle),
                                        stored.ne2_list[ne2_idx]['id'],
                                        stored.o_list[o_idx]['id'],
                                        stored.c1['id'],
                                        stored.o1['id'],
                                        False
                                    ])

                # Optional: save .pse and .txt (for manual inspection)
                if results:
                    if SAVE_PSE:
                        passing_models = set(results.keys())
                        for m in retained_models:
                            if m not in passing_models:
                                cmd.delete(m)
                        pse_file = os.path.join(output_dir, f"{protein1}_file_{file_num}.pse")
                        cmd.save(pse_file)
                    txt_file = os.path.join(output_dir, f"{protein1}_file_{file_num}.txt")
                    with open(txt_file, "w") as f:
                        for model, infos in results.items():
                            mode_num = int(model.split('_')[1])
                            affinity = affinity_dict.get(mode_num, float('nan'))
                            f.write(f"Model: {model}, Affinity: {affinity:.3f} kcal/mol\n")
                            for info in infos:
                                f.write(
                                    f" dist_NE2_O: {info['dist_NE2_O']:.2f}, "
                                    f"dist_C1_O: {info['dist_C1_O']:.2f}, angle: {info['angle']:.2f}\n"
                                    f" C1: {info['c1_info']}\n"
                                    f" O1: {info['o1_info']}\n"
                                    f" NE2: {info['ne2_info']}\n"
                                    f" O(nuc): {info['o_info']}\n"
                                )
                    print(f"Saved gated results for {protein1} to {txt_file}{' and .pse' if SAVE_PSE else ''}")
                else:
                    print(f"Skipping {protein1}: no models satisfy distance/angle conditions.")

            except Exception as e:
                print(f"Error processing {protein1} with file_num {file_num}: {str(e)}")
                continue

        # Close CSV for this ligand
        csv_fh.close()
        print(f"Written CSV for file_num {file_num}: {csv_path_this}")

finally:
    # Close PyMOL
    try:
        cmd.quit()
    except Exception:
        pass

print("All protein files processed for all file_nums.")
print(f"Per-ligand gating CSVs written to: {csv_dir}")
