from pathlib import Path
import math

def get_pymol():
    try:
        import pymol2
        pm = pymol2.PyMOL()
        pm.start()
        return pm.cmd, pm
    except Exception:
        import pymol
        pymol.finish_launching(["pymol", "-cq"])
        return pymol.cmd, None

def parse_vina_out(out_file: str, affinity_threshold: float):
    modes = []
    mode2aff = {}
    p = Path(out_file)
    if not p.exists():
        return modes, mode2aff
    with open(p, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()
    table_started = False
    for line in lines:
        s = line.strip()
        if s.startswith("mode"):
            table_started = True
            continue
        if not table_started:
            continue
        if not s or s.startswith("---"):
            continue
        parts = s.split()
        if len(parts) < 2:
            continue
        try:
            mode = int(parts[0])
            aff = float(parts[1])
        except Exception:
            continue
        mode2aff[mode] = aff
        if aff <= affinity_threshold:
            modes.append(mode)
    return modes, mode2aff

def split_ligand_states(cmd, ligand_obj="ligand", prefix="ligand_"):
    cmd.split_states(ligand_obj)
    total = cmd.count_states(ligand_obj)
    out = []
    for i in range(1, total + 1):
        old = f"{ligand_obj}_{str(i).zfill(4)}"
        if old in cmd.get_names("objects"):
            new = f"{prefix}{str(i).zfill(4)}"
            cmd.set_name(old, new)
            out.append(new)
    try:
        cmd.delete(ligand_obj)
    except Exception:
        pass
    return out

def dist(a, b):
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)

def angle_three_points(a, b, c):
    bax = a[0] - b[0]
    bay = a[1] - b[1]
    baz = a[2] - b[2]
    bcx = c[0] - b[0]
    bcy = c[1] - b[1]
    bcz = c[2] - b[2]
    dot = bax * bcx + bay * bcy + baz * bcz
    na = math.sqrt(bax * bax + bay * bay + baz * baz) + 1e-12
    nb = math.sqrt(bcx * bcx + bcy * bcy + bcz * bcz) + 1e-12
    cosv = max(-1.0, min(1.0, dot / (na * nb)))
    return math.degrees(math.acos(cosv))

def calculate_angle(vec1, vec2):
    import numpy as np
    v1 = np.array(vec1)
    v2 = np.array(vec2)
    n1 = np.linalg.norm(v1)
    n2 = np.linalg.norm(v2)
    if n1 == 0 or n2 == 0:
        return None
    cosv = np.clip(np.dot(v1, v2) / (n1 * n2), -1.0, 1.0)
    return float(np.degrees(np.arccos(cosv)))

def detect_clashes_pairs(cmd, protein_sel, ligand_sel, cutoff):
    protein_sel_heavy = f"({protein_sel}) and not hydro"
    ligand_sel_heavy = f"({ligand_sel}) and not hydro"
    try:
        clashes = cmd.find_pairs(protein_sel_heavy, ligand_sel_heavy, cutoff=cutoff, mode=0)
        return len(clashes) > 0
    except Exception:
        return False

def has_clash_around(cmd, pose_obj, cutoff, protein_obj="protein1"):
    sel = f"{protein_obj} and not hydro and (({pose_obj} and not hydro) around {cutoff})"
    try:
        return cmd.count_atoms(sel) > 0
    except Exception:
        return False

def save_filtered_pse(cmd, passing_objects, out_path):
    current = set(cmd.get_names("objects"))
    for obj in list(current):
        if obj.startswith("ligand_") or obj.startswith("pose_") or obj.isdigit():
            if obj not in passing_objects:
                try:
                    cmd.delete(obj)
                except Exception:
                    pass
    if passing_objects:
        cmd.save(str(out_path))
