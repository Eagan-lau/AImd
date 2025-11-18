#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ACT scoring pipeline (CSV -> multi-level scores), aligned with UGT/CYP style.

Inputs (from your ACT gating script CSV):
  Columns expected:
    Ligand_id, protein_id, conformation, mode, affinity_kcal,
    nu_type, nu_atom_id, dist_NE2_Nu, dist_Nu_S1P, dist_Nu_C9,
    BD_angle_deg, sele_name, clash_flag

Outputs (per input CSV):
  - <base>__pose_scores.csv             (row-level poses with scores)
  - <base>__pose_scores_by_mode.csv     (aggregated per (Ligand,protein,conf,mode))
  - <base>__conformation_scores.csv     (max pose per conformation)
  - <base>__enzyme_scores.csv           (protein-level with min–max normalization)
  - <base>__best_pose_per_protein.csv   (trace-back of best protein pose)
"""

import os
import glob
import numpy as np
import pandas as pd

# ====================== Fixed paths (edit if needed) ======================
INPUT_DIR  = "/public/agis/yanjianbin_group/liuyugeng/taxus_alphaflow_docking/ACT/yixianji/data_output/6_CSV_OUTPUT/gating_csv"
OUTPUT_DIR = "/public/agis/yanjianbin_group/liuyugeng/taxus_alphaflow_docking/ACT/yixianji/data_output/6_CSV_OUTPUT/score_output"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ====================== Scoring hyper-parameters ==========================
# Energy (two-stage, same as UGT/CYP):
E_STAGE1 = -7.0      # ΔG <= -7  -> 1.0
E_STAGE2 = -3.0      # -7..-3    -> linear to 0

# Distances (two-stage for each, tuned to ACT gating windows)
# score_d = 1 when d <= D_OPT; linear to 0 at D_MAX; zero if d > D_MAX
NE2_NU_D_OPT = 3.2
NE2_NU_D_MAX = 5.0

NU_C9_D_OPT  = 3.2
NU_C9_D_MAX  = 5.0

# Angle target (kept for reference; piecewise scoring below does not use sigma)
THETA0 = 105.0

# Pose weighting (sum to 1.0)
W_E = 0.30           # energy
W_D = 0.40           # combined distance block
W_A = 0.30           # angle

# Internal mixing inside distance block (sum to 1.0) -- S1P removed
WD_NE2 = 0.50
WD_C9  = 0.50

# Protein-level mixing (same style as before)
ALPHA        = 0.30   # coverage bonus
TOTAL_CONFS  = 6      # expected conformations per protein (adjust if different)
COVER_T      = 70.0   # conformation score threshold for coverage (% of 100)

# ====================== Scoring functions ================================
def e_score_two_stage(dg: float) -> float:
    """ΔG two-stage: <= E_STAGE1 -> 1; E_STAGE1..E_STAGE2 -> linear; else 0."""
    if pd.isna(dg):
        return 0.0
    x = float(dg)
    if x <= E_STAGE1:
        return 1.0
    if x <= E_STAGE2:
        return max(0.0, min(1.0, (E_STAGE2 - x) / (E_STAGE2 - E_STAGE1)))
    return 0.0

def d_score_two_stage(d: float, d_opt: float, d_max: float) -> float:
    """Two-stage distance: <= d_opt -> 1; d_opt..d_max -> linear to 0; >d_max -> 0."""
    if pd.isna(d):
        return 0.0
    x = float(d)
    if x <= d_opt:
        return 1.0
    if x <= d_max and d_max > d_opt:
        return max(0.0, min(1.0, (d_max - x) / (d_max - d_opt)))
    return 0.0

def a_score_piecewise(theta: float) -> float:
    """
    Piecewise angle score:
      100 <= theta <= 110  -> 1.0
       85 <= theta < 100   -> linear 0..1  (85->0, 100->1)
      110 <  theta <= 125  -> linear 1..0  (110->1, 125->0)
      else                 -> 0.0
    No 180° wrapping, uses |theta - 105| implicitly via the intervals.
    """
    if pd.isna(theta):
        return 0.0
    t = float(theta)
    if 100.0 <= t <= 110.0:
        return 1.0
    if 85.0 <= t < 100.0:
        return (t - 85.0) / (100.0 - 85.0)
    if 110.0 < t <= 125.0:
        return (125.0 - t) / (125.0 - 110.0)
    return 0.0

# ====================== Core per-file scoring ============================
def score_one_csv(in_csv: str, out_dir: str):
    base = os.path.splitext(os.path.basename(in_csv))[0]  # e.g., file_157
    df = pd.read_csv(in_csv)

    needed = [
        "Ligand_id","protein_id","conformation","mode","affinity_kcal",
        "dist_NE2_Nu","dist_Nu_C9","BD_angle_deg"
    ]
    for c in needed:
        if c not in df.columns:
            raise ValueError(f"[{base}] Missing column: {c}")

    # ------------- Pose-level scores (row-wise) -------------
    E = df["affinity_kcal"].apply(e_score_two_stage).astype(float)
    D_ne2  = df["dist_NE2_Nu"].apply(lambda x: d_score_two_stage(x, NE2_NU_D_OPT, NE2_NU_D_MAX)).astype(float)
    D_c9   = df["dist_Nu_C9"].apply( lambda x: d_score_two_stage(x, NU_C9_D_OPT,  NU_C9_D_MAX )).astype(float)

    D_comb = WD_NE2 * D_ne2 + WD_C9 * D_c9
    A = df["BD_angle_deg"].apply(a_score_piecewise).astype(float)

    s_pose = 100.0 * (W_E * E + W_D * D_comb + W_A * A)

    pose_df = df.copy()
    pose_df["E_score"]      = E
    pose_df["D_NE2_score"]  = D_ne2
    pose_df["D_C9_score"]   = D_c9
    pose_df["D_comb"]       = D_comb
    pose_df["A_score"]      = A
    pose_df["s_pose"]       = s_pose

    pose_out = os.path.join(out_dir, f"{base}__pose_scores.csv")
    pose_df.to_csv(pose_out, index=False)

    # ------------- Aggregate per (Ligand, protein, conformation, mode) -------------
    bymode = (pose_df
              .groupby(["Ligand_id","protein_id","conformation","mode"], as_index=False)
              .agg(s_pose=("s_pose","max"),
                   affinity_kcal=("affinity_kcal","min")))
    bymode_out = os.path.join(out_dir, f"{base}__pose_scores_by_mode.csv")
    bymode.to_csv(bymode_out, index=False)

    # ------------- Conformation-level: s_r = max pose per (Lig,prot,conf) -------------
    conf_df = (bymode
               .groupby(["Ligand_id","protein_id","conformation"], as_index=False)
               .agg(s_r=("s_pose","max"),
                    n_pose=("s_pose","size")))
    conf_out = os.path.join(out_dir, f"{base}__conformation_scores.csv")
    conf_df.to_csv(conf_out, index=False)

    # ------------- Protein-level: S_raw + min–max normalization -------------
    enz_rows = []
    best_pose_rows = []

    for (lig, pid), g in conf_df.groupby(["Ligand_id","protein_id"], sort=False):
        if g.empty:
            continue

        idx = g["s_r"].astype(float).idxmax()
        best_conf = g.loc[idx, "conformation"]
        best_sr   = float(g.loc[idx, "s_r"])

        n_pass = int((g["s_r"].astype(float) >= COVER_T).sum())
        coverage = n_pass / float(TOTAL_CONFS) if TOTAL_CONFS > 0 else 0.0

        protein_score_raw = best_sr * (1.0 + ALPHA * coverage)

        sub = bymode[(bymode["Ligand_id"]==lig) &
                     (bymode["protein_id"]==pid) &
                     (bymode["conformation"]==best_conf)]
        if not sub.empty:
            prow = sub.loc[sub["s_pose"].astype(float).idxmax()]
            detail = (pose_df[(pose_df["Ligand_id"]==lig) &
                              (pose_df["protein_id"]==pid) &
                              (pose_df["conformation"]==best_conf) &
                              (pose_df["mode"]==prow["mode"])]
                      .sort_values("s_pose", ascending=False)
                      .head(1))
            if not detail.empty:
                drow = detail.iloc[0]
                best_pose_rows.append({
                    "Ligand_id": lig,
                    "protein_id": pid,
                    "best_conformation": best_conf,
                    "best_mode": int(prow["mode"]),
                    "s_pose": float(prow["s_pose"]),
                    "affinity_kcal": float(prow["affinity_kcal"]),
                    "dist_NE2_Nu": float(drow["dist_NE2_Nu"]),
                    "dist_Nu_C9":  float(drow["dist_Nu_C9"]),
                    "BD_angle_deg": float(drow["BD_angle_deg"]),
                    "E_score": float(drow["E_score"]),
                    "D_NE2_score": float(drow["D_NE2_score"]),
                    "D_C9_score": float(drow["D_C9_score"]),
                    "D_comb": float(drow["D_comb"]),
                    "A_score": float(drow["A_score"])
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

    if not enzyme_df.empty:
        v = enzyme_df["protein_score_raw"].astype(float)
        v_min, v_max = float(v.min()), float(v.max())
        if v_max > v_min:
            enzyme_df["protein_score_norm"] = 100.0 * (v - v_min) / (v_max - v_min)
        else:
            enzyme_df["protein_score_norm"] = 100.0
        enzyme_df = enzyme_df.sort_values(["protein_score_norm","max_s_r"], ascending=False)

    enz_out = os.path.join(out_dir, f"{base}__enzyme_scores.csv")
    enzyme_df.to_csv(enz_out, index=False)

    best_pose_df = pd.DataFrame(best_pose_rows)
    best_pose_out = os.path.join(out_dir, f"{base}__best_pose_per_protein.csv")
    best_pose_df.to_csv(best_pose_out, index=False)

    print(f"[OK] {base}:")
    print("  -", pose_out)
    print("  -", bymode_out)
    print("  -", conf_out)
    print("  -", enz_out)
    print("  -", best_pose_out)

# ====================== Batch over all gating CSVs =======================
def main():
    csv_files = sorted(glob.glob(os.path.join(INPUT_DIR, "*.csv")))
    if not csv_files:
        print(f"[Exit] No CSV found under: {INPUT_DIR}")
        return
    for csv_path in csv_files:
        try:
            score_one_csv(csv_path, OUTPUT_DIR)
        except Exception as e:
            print(f"[Error] {os.path.basename(csv_path)}: {e}")

if __name__ == "__main__":
    main()
