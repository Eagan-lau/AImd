#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Batch scoring for CYP450 C–H gating CSVs (UGT-style math)
- Processes all *.csv under INPUT_DIR, writes per-file outputs to OUTPUT_DIR.
- Math: E(two-segment linear), D(two-segment linear on Fe–H), A(gaussian, target 0°)
- Aggregation: atom -> pose(best) -> conformation(best) -> protein(coverage boost, TOTAL_CONFS=6)
- Normalization: per-file min–max to 0–100 (UGT style)

Required input columns (from gating step):
  Ligand_id, protein_id, conformation, mode, affinity_kcal,
  dist_FE_H, angle_deg, H_atom_id, C_atom_id, fe_atom_id,
  dist_fe_sg_A, prox_sg_found, clash_flag
"""

import os
import glob
import math
import numpy as np
import pandas as pd

# =========================
# >>>>>> EDIT PATHS HERE <<<<<<
# =========================
INPUT_DIR  = "/public/agis/yanjianbin_group/liuyugeng/taxus_alphaflow_docking/CYP450/data_output/6_CSV_OUTPUT/gating_csv"
OUTPUT_DIR = "/public/agis/yanjianbin_group/liuyugeng/taxus_alphaflow_docking/CYP450/data_output/6_CSV_OUTPUT/scoring_csv"

# =========================
# UGT-style Scoring Params
# =========================
# Energy E: two-segment linear (more negative ΔG is better)
E_MIN = -7.0   # ΔG <= E_MIN  → E=1
E_MAX = -3.0   # ΔG >= E_MAX  → E=0

# Distance D: two-segment linear on Fe–H
D_A = 4.0      # d <= 4.0  → D=1
D_B = 5.0      # d >= 5.0  → D=0

# Angle A: symmetric Gaussian (target 0°, sigma tunable)
ANGLE_TARGET    = 0.0
ANGLE_FLAT_DEG  = 10.0
ANGLE_SIGMA_DEG = 12.0
ANGLE_LINEAR_MAX = 30.0     # only used if switching to linear angle scoring
ANGLE_USE_GAUSSIAN = True

# Weights (atom-level)
W_E, W_A, W_D = 0.30, 0.40, 0.30

# Coverage / protein-level
SR_CONF_THRESH = 70.0   # S_conf ≥ 70 counts to coverage
ALPHA_COVER    = 0.30   # coverage boost coefficient
TOTAL_CONFS    = 6      # fixed denominator

EPS = 1e-8


# =========================
# Scoring kernels (UGT-style)
# =========================
def score_energy(dg: float) -> float:
    """Two-segment linear energy score."""
    if pd.isna(dg):
        return 0.0
    x = float(dg)
    if x <= E_MIN:
        return 1.0
    if x >= E_MAX:
        return 0.0
    return (E_MAX - x) / (E_MAX - E_MIN + EPS)

def score_distance_linear(d: float, a: float = D_A, b: float = D_B) -> float:
    """
    Two-segment linear distance score:
      d <= a  → 1.0
      a < d < b → linearly decreases to 0.0
      d >= b  → 0.0
    """
    if pd.isna(d):
        return 0.0
    x = float(d)
    if x <= a:
        return 1.0
    if x >= b:
        return 0.0
    return (b - x) / (b - a + EPS)

def score_angle(theta_deg: float) -> float:
    """Symmetric-Gaussian angle score with flat buffer near 0°."""
    if pd.isna(theta_deg):
        return 0.0
    theta = max(0.0, float(theta_deg))  # 已对称化
    if ANGLE_USE_GAUSSIAN:
        if theta <= ANGLE_FLAT_DEG:
            return 1.0
        else:
            d = theta - ANGLE_FLAT_DEG
            return math.exp(-0.5 * (d / ANGLE_SIGMA_DEG) ** 2)
    else:
        return max(0.0, min(1.0, (ANGLE_LINEAR_MAX - theta) / ANGLE_LINEAR_MAX))

def score_atom_row(row: pd.Series) -> dict:
    """Return E, A, D, and S_atom (0–100) for a single atom-row."""
    E = score_energy(row["affinity_kcal"])
    D = score_distance_linear(row["dist_FE_H"], D_A, D_B)
    A = score_angle(row["angle_deg"])
    S_atom = 100.0 * (W_E * E + W_A * A + W_D * D)
    return {"E": E, "A": A, "D": D, "S_atom": S_atom}


# =========================
# Per-file pipeline
# =========================
def score_one_csv(input_csv: str, output_dir: str):
    base = os.path.splitext(os.path.basename(input_csv))[0]
    pose_csv = os.path.join(output_dir, f"{base}_pose_scores.csv")
    conf_csv = os.path.join(output_dir, f"{base}_conf_scores.csv")
    prot_csv = os.path.join(output_dir, f"{base}_protein_scores.csv")

    try:
        df = pd.read_csv(input_csv)
    except Exception as e:
        print(f"[Error] Read fail: {input_csv} -> {e}")
        for p in (pose_csv, conf_csv, prot_csv):
            pd.DataFrame().to_csv(p, index=False)
        return

    if df.empty:
        print(f"[Warn] Empty CSV: {input_csv}")
        for p in (pose_csv, conf_csv, prot_csv):
            pd.DataFrame().to_csv(p, index=False)
        return

    need_cols = [
        "Ligand_id","protein_id","conformation","mode","affinity_kcal",
        "dist_FE_H","angle_deg","H_atom_id","C_atom_id","fe_atom_id",
        "dist_fe_sg_A","prox_sg_found","clash_flag"
    ]
    missing = [c for c in need_cols if c not in df.columns]
    if missing:
        print(f"[Error] Missing columns in {input_csv}: {missing}")
        for p in (pose_csv, conf_csv, prot_csv):
            pd.DataFrame().to_csv(p, index=False)
        return

    # Remove clashes if present (safety)
    if "clash_flag" in df.columns:
        df = df[df["clash_flag"] == False].copy()

    # Convert Fe–SG distance to numeric (may be 'NA')
    if "dist_fe_sg_A" in df.columns:
        df["dist_fe_sg_A"] = pd.to_numeric(df["dist_fe_sg_A"], errors="coerce")

    # ---- Atom-level scores ----
    atom_scores = df.apply(score_atom_row, axis=1, result_type="expand")
    df = pd.concat([df, atom_scores], axis=1)

    # ---- Pose level: best atom per pose; tie-break by shorter Fe–H, smaller angle ----
    pose_keys = ["Ligand_id", "protein_id", "conformation", "mode"]
    df_sorted = df.sort_values(
        by=pose_keys + ["S_atom", "dist_FE_H", "angle_deg"],
        ascending=[True, True, True, True, False, True, True]
    )
    pose_df = (df_sorted
               .drop_duplicates(subset=pose_keys, keep="first")
               .rename(columns={"S_atom": "S_pose_best"})
               .reset_index(drop=True))

    pose_out = pose_df[[
        "Ligand_id", "protein_id", "conformation", "mode",
        "H_atom_id", "C_atom_id", "fe_atom_id",
        "affinity_kcal", "dist_FE_H", "angle_deg",
        "E", "D", "A", "S_pose_best",
        "prox_sg_found", "dist_fe_sg_A"
    ]]
    pose_out.to_csv(pose_csv, index=False)

    # ---- Conformation level: best pose per (Ligand, protein, conformation) ----
    conf_keys = ["Ligand_id", "protein_id", "conformation"]
    conf_sorted = pose_out.sort_values(
        by=conf_keys + ["S_pose_best", "dist_FE_H", "angle_deg"],
        ascending=[True, True, True, False, True, True]
    )
    conf_df = (conf_sorted
               .drop_duplicates(subset=conf_keys, keep="first")
               .rename(columns={
                   "mode": "best_mode_in_conf",
                   "H_atom_id": "best_H_id_in_conf",
                   "C_atom_id": "best_C_id_in_conf",
                   "S_pose_best": "S_conf"
               })
               .reset_index(drop=True))
    conf_df.to_csv(conf_csv, index=False)

    # ---- Protein level: best conformation + coverage (TOTAL_CONFS=6) ----
    prot_keys = ["Ligand_id", "protein_id"]

    cover_tbl = (conf_df.assign(pass_conf=(conf_df["S_conf"] >= SR_CONF_THRESH).astype(int))
                 .groupby(prot_keys)["pass_conf"].sum().reset_index(name="n_cover"))

    prot_sorted = conf_df.sort_values(
        by=prot_keys + ["S_conf", "dist_FE_H", "angle_deg"],
        ascending=[True, True, False, True, True]
    )
    prot_best_conf = prot_sorted.drop_duplicates(subset=prot_keys, keep="first").copy()

    prot_df = prot_best_conf.merge(cover_tbl, on=prot_keys, how="left")
    prot_df["n_cover"] = prot_df["n_cover"].fillna(0).astype(int)
    prot_df["coverage"] = prot_df["n_cover"] / float(TOTAL_CONFS)

    prot_df["S_protein_raw"] = prot_df["S_conf"] * (1.0 + ALPHA_COVER * prot_df["coverage"])

    # ---- per-file normalization (UGT style) ----
    v = prot_df["S_protein_raw"].astype(float)
    vmin, vmax = float(v.min()), float(v.max())
    prot_df["S_protein_norm"] = 100.0 if vmax <= vmin else 100.0 * (v - vmin) / (vmax - vmin)

    protein_out = prot_df[[
        "Ligand_id", "protein_id",
        "conformation",            # best conformation
        "best_mode_in_conf",       # best pose
        "best_H_id_in_conf", "best_C_id_in_conf",
        "affinity_kcal", "dist_FE_H", "angle_deg",
        "S_conf", "n_cover", "coverage",
        "S_protein_raw", "S_protein_norm",
        "prox_sg_found", "dist_fe_sg_A"
    ]].rename(columns={
        "conformation": "best_conformation",
        "affinity_kcal": "best_affinity_in_protein",
        "dist_FE_H": "best_dist_FEH_in_protein",
        "angle_deg": "best_angle_in_protein"
    })
    protein_out.to_csv(prot_csv, index=False)

    print(f"[OK] {base} ->")
    print(f"   Pose scores : {pose_csv}")
    print(f"   Conf scores : {conf_csv}")
    print(f"   Protein     : {prot_csv}")


# =========================
# Batch main
# =========================
def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    csv_list = sorted(glob.glob(os.path.join(INPUT_DIR, "*.csv")))
    if not csv_list:
        print(f"[Warn] No CSV files found in: {INPUT_DIR}")
        return
    for csv_path in csv_list:
        try:
            score_one_csv(csv_path, OUTPUT_DIR)
        except Exception as e:
            print(f"[Error] {os.path.basename(csv_path)}: {e}")

if __name__ == "__main__":
    main()
