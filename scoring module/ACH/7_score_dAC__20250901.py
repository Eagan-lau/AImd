#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Batch scoring for Hydrolase (Ser–His–Asp/Glu) gating CSVs, UGT-style math.

For each input CSV under INPUT_DIR:
- Pose:    E(two-stage), D(weighted two-stage on two distances), A(Gaussian, BD angle)
           s_pose = 100 * (W_E*E + W_D*D + W_A*A)
- Conformation: s_r = max(s_pose) per (Ligand_id, protein_id, conformation)
- Protein: S_raw = best(s_r) * (1 + ALPHA * coverage)
           coverage = (# s_r >= COVER_T) / TOTAL_CONFS
- Normalize protein scores (min–max) to 0–100 within each CSV.
- Write 4 output CSVs per input, using the input file basename as prefix.

Expected columns in input CSVs produced by your hydrolase gating script:
  Ligand_id, protein_id, conformation, mode, affinity_kcal,
  carbon_index, selection_name,
  dist_to_SerOG_A, dist_SerOG_HisNE2_A, bd_angle_deg, ligand_atom_ids
"""

import os
import glob
import numpy as np
import pandas as pd

# ================= I/O =================
INPUT_DIR  = "/public/agis/yanjianbin_group/liuyugeng/taxus_alphaflow_docking/dAC/data_output/6_CSV_OUTPUT/gating_csv/"
OUTPUT_DIR = "/public/agis/yanjianbin_group/liuyugeng/taxus_alphaflow_docking/dAC/data_output/6_CSV_OUTPUT/"

# ================= Weights & global parameters (UPDATED) =================
W_E   = 0.30     # energy weight
W_D   = 0.50     # distance weight
W_A   = 0.20     # angle weight

# Energy scoring thresholds (same as UGT)
E_STRICT = -7.0  # <= -7 -> 1.0
E_LOOSE  = -3.0  # -7..-3 linear -> 0.0, >-3 -> 0.0

# Distance scoring (two-stage): default 3/5 Å window (UGT-style)
D_BEST   = 3.2   # <=3 -> 1.0
D_CUTOFF = 5.0   # 3..5 linear -> 0.0; >5 -> 0.0

# Optional tighter window for SerOG–HisNE2
D2_TIGHT_BEST   = 3.2
D2_TIGHT_CUTOFF = 5.0

# Angle scoring (Gaussian) for Bürgi–Dunitz angle at carbonyl addition
THETA0_BD = 107.0  # target BD angle (degrees)
SIGMA_BD  = 20.0   # <-- updated sigma

# Protein-layer coverage boost (same scheme as UGT)
ALPHA        = 0.30     # coverage boost factor
TOTAL_CONFS  = 6        # expected total conformations per protein
COVER_T      = 70.0     # per-conformation pass threshold (0–100 scale)

# Distance weights inside D (primary geometry vs. triad tightness)
W_D_CARBONYL = 0.90     # weight for d(C···Ser-OG)
W_D_TRIAD    = 0.10     # weight for d(Ser-OG···His-NE2)
W_D_SUM      = W_D_CARBONYL + W_D_TRIAD


# ================= Scoring helpers =================
def e_score_two_stage(dg: float) -> float:
    """Energy score in [0,1]: <=E_STRICT ->1; E_STRICT..E_LOOSE linear->0; >E_LOOSE ->0."""
    if pd.isna(dg): return 0.0
    x = float(dg)
    if x <= E_STRICT: return 1.0
    if x <= E_LOOSE:  return max(0.0, min(1.0, (E_LOOSE - x) / (E_LOOSE - E_STRICT)))
    return 0.0

def d_score_two_stage_default(d: float) -> float:
    """Distance score in [0,1]: <=3 ->1; 3..5 linear->0; >5 ->0."""
    if pd.isna(d): return 0.0
    x = float(d)
    if x <= D_BEST:   return 1.0
    if x <= D_CUTOFF: return max(0.0, min(1.0, (D_CUTOFF - x) / (D_CUTOFF - D_BEST)))
    return 0.0

def d_score_two_stage_tight(d: float) -> float:
    """Tighter two-stage (e.g., for SerOG–HisNE2): <=3.5 ->1; 3.5..5 ->0; >5 ->0."""
    if pd.isna(d): return 0.0
    x = float(d)
    if x <= D2_TIGHT_BEST:   return 1.0
    if x <= D2_TIGHT_CUTOFF: return max(0.0, min(1.0, (D2_TIGHT_CUTOFF - x) / (D2_TIGHT_CUTOFF - D2_TIGHT_BEST)))
    return 0.0

def D_from_two_dist_weighted(d_c_serog, d_serhis, use_tight_for_serhis=False) -> float:
    """Weighted average of two distance scores with optional tighter function for Ser–His."""
    s1 = d_score_two_stage_default(d_c_serog)
    s2 = (d_score_two_stage_tight(d_serhis) if use_tight_for_serhis
          else d_score_two_stage_default(d_serhis))
    wsum = W_D_SUM if W_D_SUM > 0 else 1.0
    return (W_D_CARBONYL * s1 + W_D_TRIAD * s2) / wsum

def a_score_gauss_bd(theta, theta0=THETA0_BD, sigma=SIGMA_BD, fold180=False) -> float:
    """
    Gaussian angle score in [0,1].
    For hydrolase BD angle, we DO NOT fold with min(d, 180-d) by default.
    """
    if pd.isna(theta): return 0.0
    t = float(theta)
    d = abs(t - theta0)
    if fold180:
        d = min(d, 180.0 - d)
    return float(np.exp(-0.5 * (d / sigma) ** 2))


# ================= Core scoring for one CSV =================
def score_one_csv(input_csv: str, out_dir: str, use_tight_serhis=False):
    os.makedirs(out_dir, exist_ok=True)

    base = os.path.splitext(os.path.basename(input_csv))[0]
    pose_out       = os.path.join(out_dir, f"{base}_pose_scores.csv")
    conf_out       = os.path.join(out_dir, f"{base}_conformation_scores.csv")
    enz_out        = os.path.join(out_dir, f"{base}_enzyme_scores.csv")
    best_pose_out  = os.path.join(out_dir, f"{base}_best_pose_per_protein.csv")

    # Read input
    df = pd.read_csv(input_csv)
    need_cols = [
        "Ligand_id","protein_id","conformation","mode","affinity_kcal",
        "dist_to_SerOG_A","dist_SerOG_HisNE2_A","bd_angle_deg"
    ]
    for c in need_cols:
        if c not in df.columns:
            raise ValueError(f"[{base}] Missing column in CSV: {c}")

    # --- Robust conformation normalization ---
    # 
    df["protein_id"]   = df["protein_id"].astype(str).str.strip()
    s = df["conformation"].astype(str).str.strip()

    #
    is_digit        = s.str.fullmatch(r"\d+")
    conf_from_self  = s.str.extract(r'_(\d+)$', expand=False)
    conf_from_pid   = df["protein_id"].str.extract(r'_(\d+)$', expand=False)

    df["conformation"] = (
        s.where(is_digit, np.nan)            # 
          .fillna(conf_from_self)            # 
          .fillna(conf_from_pid)             # 
          .fillna("0")                       # 
    ).astype(str)


    if df.empty:
        # Write empty shells for consistency
        pose_cols = need_cols + ["E_score","D_score","A_score","s_pose"]
        pd.DataFrame(columns=pose_cols).to_csv(pose_out, index=False)
        pd.DataFrame(columns=["Ligand_id","protein_id","conformation","s_r","n_pose"]).to_csv(conf_out, index=False)
        pd.DataFrame(columns=[
            "Ligand_id","protein_id","protein_score_raw","best_conformation",
            "coverage","n_pass_conformations","total_conformations",
            "max_s_r","mean_s_r","protein_score_norm"
        ]).to_csv(enz_out, index=False)
        pd.DataFrame(columns=[]).to_csv(best_pose_out, index=False)
        print(f"[{base}] Empty input; wrote empty outputs.")
        return

    # -------- Pose layer --------
    parts = []
    for conf, g in df.groupby("conformation", sort=False):
        E = g["affinity_kcal"].apply(e_score_two_stage).astype(float).values

        D = [
            D_from_two_dist_weighted(d1, d2, use_tight_for_serhis=use_tight_serhis)
            for d1, d2 in zip(g["dist_to_SerOG_A"], g["dist_SerOG_HisNE2_A"])
        ]

        # BD angle (no 180-folding)
        A = g["bd_angle_deg"].apply(lambda a: a_score_gauss_bd(a, THETA0_BD, SIGMA_BD, fold180=False)).astype(float).values

        s_pose = 100.0 * (W_E*np.asarray(E) + W_D*np.asarray(D) + W_A*np.asarray(A))

        gg = g.copy()
        gg["E_score"] = E
        gg["D_score"] = D
        gg["A_score"] = A
        gg["s_pose"]  = s_pose
        parts.append(gg)

    pose_df = (pd.concat(parts, ignore_index=True)
               if parts else pd.DataFrame(columns=df.columns.tolist()+["E_score","D_score","A_score","s_pose"]))
    pose_df.to_csv(pose_out, index=False)

    # -------- Conformation layer --------
    conf_df = (pose_df.groupby(["Ligand_id","protein_id","conformation"], as_index=False)
                      .agg(s_r=("s_pose","max"),
                           n_pose=("s_pose","size")))
    conf_df.to_csv(conf_out, index=False)

    # -------- Protein layer (raw + min–max) --------
    enz_rows = []
    best_pose_rows = []

    for (lig, pid), g in conf_df.groupby(["Ligand_id","protein_id"], sort=False):
        if g.empty:
            continue

        # Best conformation
        idx = g["s_r"].astype(float).idxmax()
        best_conf = g.loc[idx, "conformation"]
        best_sr   = float(g.loc[idx, "s_r"])

        # Coverage
        n_pass = int((g["s_r"].astype(float) >= COVER_T).sum())
        coverage = n_pass / float(TOTAL_CONFS) if TOTAL_CONFS > 0 else 0.0

        # Raw protein score
        protein_score_raw = best_sr * (1.0 + ALPHA * coverage)

        # Backtrack best pose within best conformation
        sub = pose_df[(pose_df["Ligand_id"]==lig) &
                      (pose_df["protein_id"]==pid) &
                      (pose_df["conformation"]==best_conf)]
        if not sub.empty:
            prow = sub.loc[sub["s_pose"].astype(float).idxmax()]
            best_pose_rows.append({
                "Ligand_id": lig,
                "protein_id": pid,
                "best_conformation": best_conf,
                "best_mode": int(prow["mode"]),
                "s_pose": float(prow["s_pose"]),
                "affinity_kcal": float(prow["affinity_kcal"]),
                "dist_to_SerOG_A": float(prow["dist_to_SerOG_A"]),
                "dist_SerOG_HisNE2_A": float(prow["dist_SerOG_HisNE2_A"]),
                "bd_angle_deg": float(prow["bd_angle_deg"]),
                "E_score": float(prow["E_score"]),
                "D_score": float(prow["D_score"]),
                "A_score": float(prow["A_score"])
            })

        enz_rows.append({
            "Ligand_id": lig,
            "protein_id": pid,
            "protein_score_raw": protein_score_raw,
            "best_conformation": best_conf,
            "coverage": coverage,
            "n_pass_conformations": n_pass,
            "total_conformations": int(TOTAL_CONFS),
            "max_s_r": best_sr,
            "mean_s_r": float(np.mean(g["s_r"]))
        })

    enzyme_df = pd.DataFrame(enz_rows)

    # Min–max normalization to 0–100 per input CSV
    if not enzyme_df.empty:
        v = enzyme_df["protein_score_raw"].astype(float)
        v_min, v_max = float(v.min()), float(v.max())
        if v_max > v_min:
            enzyme_df["protein_score_norm"] = 100.0 * (v - v_min) / (v_max - v_min)
        else:
            enzyme_df["protein_score_norm"] = 100.0
        enzyme_df = enzyme_df.sort_values(["protein_score_norm","max_s_r"], ascending=False)

    enzyme_df.to_csv(enz_out, index=False)
    pd.DataFrame(best_pose_rows).to_csv(best_pose_out, index=False)

    print(f"[OK] {base} ->")
    print("  -", pose_out)
    print("  -", conf_out)
    print("  -", enz_out)
    print("  -", best_pose_out)


# ================= Batch entry =================
def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    csv_list = sorted(glob.glob(os.path.join(INPUT_DIR, "*.csv")))
    if not csv_list:
        print(f"No CSV files found in: {INPUT_DIR}")
        return
    for csv_path in csv_list:
        try:
            # Set use_tight_serhis=True if you want a tighter Ser–His distance window
            score_one_csv(csv_path, OUTPUT_DIR, use_tight_serhis=False)
        except Exception as e:
            print(f"[Error] {os.path.basename(csv_path)}: {e}")

if __name__ == "__main__":
    main()
