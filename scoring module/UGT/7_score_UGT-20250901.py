#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Batch scoring for UGT gating CSVs.

For each input CSV under INPUT_DIR:
- Pose:    E(two-stage), D(two-stage on two distances), A(gaussian)
           s_pose = 100 * (W_E*E + W_D*D + W_A*A)
- Conformation: s_r = max(s_pose) per (Ligand_id, protein_id, conformation)
- Protein: S_raw = best(s_r) * (1 + ALPHA * coverage)
           coverage = (# s_r >= COVER_T) / TOTAL_CONFS
- Normalize protein scores (min–max) to 0–100 within each CSV.
- Write 4 output CSVs per input, using the input file basename as prefix.
"""

import os
import glob
import numpy as np
import pandas as pd

# ================= Fixed parameters (edit if needed) =================
# I/O
INPUT_DIR = "/public/agis/yanjianbin_group/liuyugeng/taxus_alphaflow_docking/UDPGT/data_output/6_CSV_OUTPUT/gating_csv"
OUTPUT_DIR = "/public/agis/yanjianbin_group/liuyugeng/taxus_alphaflow_docking/UDPGT/data_output/6_CSV_OUTPUT/"

# Weights and thresholds
W_E      = 0.30        # energy weight
W_D      = 0.40        # distance weight
W_A      = 0.30        # angle weight
THETA0   = 180         # target angle (SN2/backside-like)
SIGMA    = 30.0        # Gaussian sigma for angle
ALPHA    = 0.30        # coverage boost on protein layer
TOTAL_CONFS = 6        # expected total conformations per protein
COVER_T  = 70.0        # conformation pass threshold (0–100 scale)

# ================= Scoring helpers =================
def e_score_two_stage(dg: float) -> float:
    """Energy score in [0,1]: ≤-7 →1; -7..-3 linear→0; >-3 →0."""
    if pd.isna(dg): return 0.0
    x = float(dg)
    if x <= -7.0: return 1.0
    if x <= -3.0: return max(0.0, min(1.0, (-3.0 - x) / 4.0))
    return 0.0

def d_score_two_stage(d: float) -> float:
    """Distance score in [0,1]: ≤3 →1; 3..5 linear→0; >5 →0."""
    if pd.isna(d): return 0.0
    x = float(d)
    if x <= 3.2: return 1.0
    if x <= 5.0: return max(0.0, min(1.0, (5.0 - x) / 2.0))
    return 0.0

def D_from_two_dist(d1, d2) -> float:
    """Distance score: require BOTH distances valid, else 0."""
    if pd.isna(d1) or pd.isna(d2):
        return 0.0
    s1 = d_score_two_stage(d1)
    s2 = d_score_two_stage(d2)
    return (s1 + s2) / 2.0 if (s1 > 0 and s2 > 0) else 0.0

def a_score_gauss(theta, theta0=THETA0, sigma=SIGMA, flat=25.0) -> float:
    """
    Angle score with flat-top + Gaussian tail.
    - Within ±flat degrees of theta0: full score = 1.0
    - Outside ±flat: Gaussian decay
    """
    if pd.isna(theta):
        return 0.0
    t = float(theta)
    d = abs(t - theta0)
    d = min(d, 180.0 - d)   # 保持 0–180° 对称

    if d <= flat:
        return 1.0
    else:
        dd = d - flat
        return float(np.exp(-0.5 * (dd / sigma) ** 2))


def clamp01(x):
    return max(0.0, min(1.0, x))

# ================= Core scoring for one CSV =================
def score_one_csv(input_csv: str, out_dir: str):
    os.makedirs(out_dir, exist_ok=True)

    # Determine output file prefix based on input filename
    base = os.path.splitext(os.path.basename(input_csv))[0]
    pose_out = os.path.join(out_dir, f"{base}_pose_scores.csv")
    conf_out = os.path.join(out_dir, f"{base}_conformation_scores.csv")
    enz_out  = os.path.join(out_dir, f"{base}_enzyme_scores.csv")
    best_pose_out = os.path.join(out_dir, f"{base}_best_pose_per_protein.csv")

    # Read input
    df = pd.read_csv(input_csv)
    need_cols = ["Ligand_id","protein_id","conformation","mode","affinity_kcal",
                 "dist_NE2_O","dist_C1_O","angle_deg"]
    for c in need_cols:
        if c not in df.columns:
            raise ValueError(f"[{base}] Missing column in CSV: {c}")
    if df.empty:
        # Write empty shells for consistency
        pd.DataFrame(columns=need_cols + ["E_score","D_score","A_score","s_pose"]).to_csv(pose_out, index=False)
        pd.DataFrame(columns=["Ligand_id","protein_id","conformation","s_r","n_pose"]).to_csv(conf_out, index=False)
        pd.DataFrame(columns=["Ligand_id","protein_id","protein_score_raw","best_conformation",
                              "coverage","n_pass_conformations","total_conformations",
                              "max_s_r","mean_s_r","protein_score_norm"]).to_csv(enz_out, index=False)
        pd.DataFrame(columns=[]).to_csv(best_pose_out, index=False)
        print(f"[{base}] Empty input; wrote empty outputs.")
        return

    # -------- Pose layer --------
    parts = []
    for conf, g in df.groupby("conformation", sort=False):
        E = g["affinity_kcal"].apply(e_score_two_stage).astype(float).values
        D = [D_from_two_dist(d1, d2) for d1, d2 in zip(g["dist_NE2_O"], g["dist_C1_O"])]
        A = g["angle_deg"].apply(lambda a: a_score_gauss(a, THETA0, SIGMA)).astype(float).values
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
                "dist_NE2_O": float(prow["dist_NE2_O"]),
                "dist_C1_O": float(prow["dist_C1_O"]),
                "angle_deg": float(prow["angle_deg"]),
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
            score_one_csv(csv_path, OUTPUT_DIR)
        except Exception as e:
            print(f"[Error] {os.path.basename(csv_path)}: {e}")

if __name__ == "__main__":
    main()
