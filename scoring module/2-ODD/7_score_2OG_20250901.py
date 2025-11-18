#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Fe2OG scoring script (pose -> conformation -> protein).

Update:
- Process ALL CSV files under INPUT_DIR; for each input CSV, write three output
  tables with the CSV basename as a prefix, e.g.:
    <prefix>.fe2og_pose_scores.csv
    <prefix>.fe2og_conformation_scores.csv
    <prefix>.fe2og_protein_scores.csv

Scoring rules:
- Distance subscore: 1.0 for 2.5–3.5 Å; linearly decreases to 0.0 from 3.5–5.0 Å; >=5.0 Å -> 0.0.
- Angle subscore (axis deviation, deg): 1.0 for 0–30°; linearly decreases to 0.0 from 30–60°; >=60° -> 0.0.

Assumes the gating CSV has one row per passing H atom (may contain multiple rows per pose).
"""

import os
import glob
import math
import pandas as pd
import numpy as np

# =========================
# Fixed parameters (edit here)
# =========================
INPUT_DIR = "/public/agis/yanjianbin_group/liuyugeng/taxus_alphaflow_docking/Fe2OG/data_output/6_CSV_OUTPUT/gating_csv"
OUT_DIR   = "/public/agis/yanjianbin_group/liuyugeng/taxus_alphaflow_docking/Fe2OG/data_output/6_CSV_OUTPUT/"

# Total number of conformations per protein (for coverage calculation)
TOTAL_CONFS = 6

# Protein-level aggregation: best conformation * (1 + ALPHA * coverage)
ALPHA = 0.30

# Pose subscore weights
W_E = 0.30  # energy
W_D = 0.50  # distance
W_A = 0.20  # angle

# Energy subscore (kcal/mol), linear from E_FULL to E_ZERO
E_FULL = -7.0  # <= -> 1.0
E_ZERO = -3.0  # >= -> 0.0

# Distance subscore (Å), piecewise:
# 2.5–3.5 -> 1.0; 3.5–5.0 -> linear to 0.0; >=5.0 -> 0.0
D_FULL_LOW  = 2.5
D_FULL_HIGH = 3.5
D_ZERO      = 5.0

# Angle subscore (deg, axis deviation 0..90), piecewise:
# 0–30 -> 1.0; 30–60 -> linear to 0.0; >=60 -> 0.0
A_FULL = 30.0
A_ZERO = 60.0

# Only keep rows with clash_flag == False (safety)
FILTER_CLASH = True


# =========================
# Helpers (scoring functions)
# =========================
def score_energy(dg: float) -> float:
    """Energy subscore E in [0,1], linear between E_FULL and E_ZERO."""
    if pd.isna(dg):
        return 0.0
    if dg <= E_FULL:
        return 1.0
    if dg >= E_ZERO:
        return 0.0
    return (E_ZERO - dg) / (E_ZERO - E_FULL)


def score_distance(dist: float) -> float:
    """
    Distance subscore D in [0,1]:
      2.5–3.5 Å -> 1.0
      3.5–5.0 Å -> linear to 0.0
      >=5.0 Å   -> 0.0
    Values <2.5 Å are clamped to 1.0 (gating should exclude <2.5 anyway).
    """
    if pd.isna(dist):
        return 0.0
    d = float(dist)
    if d <= D_FULL_HIGH:
        return 1.0
    if d >= D_ZERO:
        return 0.0
    # linear from 3.5 -> 5.0
    return (D_ZERO - d) / (D_ZERO - D_FULL_HIGH)


def score_angle_axis_deviation(angle_deg: float) -> float:
    """
    Angle subscore A in [0,1], using axis deviation (0..90 deg):
      0–30  -> 1.0
      30–60 -> linear to 0.0
      >=60  -> 0.0
    """
    if pd.isna(angle_deg):
        return 0.0
    a = float(angle_deg)
    if a <= A_FULL:
        return 1.0
    if a >= A_ZERO:
        return 0.0
    # linear from 30 -> 60
    return (A_ZERO - a) / (A_ZERO - A_FULL)


def minmax_0_100(series: pd.Series) -> pd.Series:
    """Min-max normalize to 0-100. If all equal or empty, return 100s (or empty)."""
    if series is None or len(series) == 0:
        return series
    minv = series.min()
    maxv = series.max()
    if pd.isna(minv) or pd.isna(maxv):
        return series
    if maxv - minv < 1e-12:
        return pd.Series([100.0] * len(series), index=series.index)
    return 100.0 * (series - minv) / (maxv - minv)


# =========================
# Core scoring per CSV
# =========================
def score_single_csv(csv_path: str, out_dir: str):
    os.makedirs(out_dir, exist_ok=True)
    prefix = os.path.splitext(os.path.basename(csv_path))[0]

    df = pd.read_csv(csv_path)

    # Column compatibility: distance column name may vary by gating version
    distance_col = None
    for candidate in ("dist_FE_H", "dist_FE_atom", "dist_FE"):
        if candidate in df.columns:
            distance_col = candidate
            break
    if distance_col is None:
        raise ValueError(f"[{prefix}] No distance column found. Expected one of: dist_FE_H, dist_FE_atom, dist_FE.")

    # Clash column compatibility: prefer clash_flag; fallback to clash.flag
    clash_col = "clash_flag" if "clash_flag" in df.columns else ("clash.flag" if "clash.flag" in df.columns else None)

    required_cols = [
        "Ligand_id", "protein_id", "conformation", "mode",
        "affinity_kcal", distance_col, "angle_deg"
    ]
    if clash_col:
        required_cols.append(clash_col)

    for c in required_cols:
        if c not in df.columns:
            raise ValueError(f"[{prefix}] Missing required column in CSV: {c}")

    # Optional filter: keep only non-clashing rows
    if FILTER_CLASH and clash_col:
        df = df[df[clash_col] == False].copy()

    # Optional filter: H atoms only (if present)
    if "lig_atom_symbol" in df.columns:
        df = df[df["lig_atom_symbol"].astype(str).str.upper() == "H"].copy()

    # Handle empty after filtering
    out_pose_file = os.path.join(out_dir, f"{prefix}.fe2og_pose_scores.csv")
    out_conf_file = os.path.join(out_dir, f"{prefix}.fe2og_conformation_scores.csv")
    out_prot_file = os.path.join(out_dir, f"{prefix}.fe2og_protein_scores.csv")

    if df.empty:
        pd.DataFrame().to_csv(out_pose_file, index=False)
        pd.DataFrame().to_csv(out_conf_file, index=False)
        pd.DataFrame().to_csv(out_prot_file, index=False)
        print(f"[{prefix}] No rows to score after filtering. Wrote empty outputs:")
        print(" -", out_pose_file)
        print(" -", out_conf_file)
        print(" -", out_prot_file)
        return

    # Subscores
    df["E_sub"] = df["affinity_kcal"].apply(score_energy)
    df["D_sub"] = df[distance_col].apply(score_distance)
    df["A_sub"] = df["angle_deg"].apply(score_angle_axis_deviation)

    # Pose score per H-row
    df["S_pose_raw"] = (W_E * df["E_sub"] +
                        W_D * df["D_sub"] +
                        W_A * df["A_sub"])

    # Reduce to one row per pose: take the H-row with max S_pose_raw
    grp_keys = ["protein_id", "conformation", "mode"]
    idx_best_row = df.groupby(grp_keys)["S_pose_raw"].idxmax()
    df_pose = df.loc[idx_best_row].copy()

    # Conformation-level: best pose per (protein_id, conformation)
    idx_best_pose_per_conf = df_pose.groupby(["protein_id", "conformation"])["S_pose_raw"].idxmax()
    cols_keep = [
        "Ligand_id", "protein_id", "conformation", "mode",
        "affinity_kcal", distance_col, "angle_deg",
        "E_sub", "D_sub", "A_sub", "S_pose_raw"
    ]
    df_conf = df_pose.loc[idx_best_pose_per_conf, cols_keep].copy()
    df_conf.rename(columns={"S_pose_raw": "S_conf_raw", "mode": "best_mode"}, inplace=True)

    # Protein-level coverage and best conformation
    conf_counts = df_conf.groupby("protein_id")["conformation"].nunique().rename("n_confs_pass")
    df_protein = conf_counts.to_frame().reset_index()
    df_protein["coverage"] = (df_protein["n_confs_pass"] / float(TOTAL_CONFS)).clip(0.0, 1.0)

    best_conf_idx = df_conf.groupby("protein_id")["S_conf_raw"].idxmax()
    df_best_conf = df_conf.loc[best_conf_idx, ["protein_id", "conformation", "S_conf_raw"]].copy()
    df_best_conf.rename(columns={"conformation": "best_conformation"}, inplace=True)

    df_protein = df_protein.merge(df_best_conf, on="protein_id", how="left")
    df_protein["S_protein_raw"] = df_protein["S_conf_raw"] * (1.0 + ALPHA * df_protein["coverage"])

    # Normalize protein-level scores to 0–100 for ranking
    df_protein["S_protein_norm"] = minmax_0_100(df_protein["S_protein_raw"])

    # ========= Outputs per input CSV =========
    # 1) Pose-level (one row per pose)
    out_pose = df_pose.copy()
    out_pose["S_pose"] = (out_pose["S_pose_raw"] * 100.0).round(2)
    out_pose["E_sub"] = out_pose["E_sub"].round(3)
    out_pose["D_sub"] = out_pose["D_sub"].round(3)
    out_pose["A_sub"] = out_pose["A_sub"].round(3)
    out_pose = out_pose.sort_values(["protein_id", "conformation", "mode"])
    out_pose.to_csv(out_pose_file, index=False)

    # 2) Conformation-level (best pose per conformation)
    out_conf = df_conf.copy()
    out_conf["S_conf"] = (out_conf["S_conf_raw"] * 100.0).round(2)
    out_conf = out_conf.sort_values(["protein_id", "conformation"])
    out_conf.to_csv(out_conf_file, index=False)

    # 3) Protein-level
    out_prot = df_protein.copy()
    out_prot["S_conf"] = (out_prot["S_conf_raw"] * 100.0).round(2)
    out_prot["S_protein_raw_0_100"] = (out_prot["S_protein_raw"] * 100.0).round(2)
    out_prot["S_protein"] = out_prot["S_protein_norm"].round(2)
    out_prot = out_prot.sort_values("S_protein", ascending=False)
    cols = ["protein_id", "best_conformation", "n_confs_pass", "coverage",
            "S_conf", "S_protein_raw_0_100", "S_protein"]
    out_prot[cols].to_csv(out_prot_file, index=False)

    print(f"[{prefix}] Saved:")
    print(" -", out_pose_file)
    print(" -", out_conf_file)
    print(" -", out_prot_file)


# =========================
# Main: iterate all CSVs in INPUT_DIR
# =========================
def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    pattern = os.path.join(INPUT_DIR, "*.csv")
    csv_files = sorted(glob.glob(pattern))
    if not csv_files:
        print(f"No CSV files found under: {INPUT_DIR}")
        return

    for csv_path in csv_files:
        try:
            score_single_csv(csv_path, OUT_DIR)
        except Exception as e:
            prefix = os.path.splitext(os.path.basename(csv_path))[0]
            print(f"[{prefix}] ERROR: {e}")


if __name__ == "__main__":
    main()
