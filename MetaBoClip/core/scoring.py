import math
import pandas as pd
import numpy as np

def energy_linear(dg, e_min=-7.0, e_max=-3.0):
    if pd.isna(dg):
        return 0.0
    x = float(dg)
    if x <= e_min:
        return 1.0
    if x >= e_max:
        return 0.0
    return (e_max - x) / (e_max - e_min)

def distance_piecewise(d, best, cutoff):
    if pd.isna(d):
        return 0.0
    x = float(d)
    if x <= best:
        return 1.0
    if x >= cutoff:
        return 0.0
    return (cutoff - x) / (cutoff - best)

def gaussian_score(theta, target, sigma, flat=0.0, fold180=False):
    if pd.isna(theta):
        return 0.0
    x = float(theta)
    d = abs(x - target)
    if fold180:
        d = min(d, 180.0 - d)
    if d <= flat:
        return 1.0
    dd = d - flat
    return math.exp(-0.5 * (dd / sigma) ** 2)

def coupled_mean_positive(a, b):
    if pd.isna(a) or pd.isna(b):
        return 0.0
    a = float(a)
    b = float(b)
    return (a + b) / 2.0 if (a > 0 and b > 0) else 0.0

def aggregate_standard(df, score_col="s_pose", total_confs=6, cover_t=70.0, alpha=0.30):
    if df.empty:
        return pd.DataFrame(), pd.DataFrame()
    conf_df = df.groupby(["Ligand_id", "protein_id", "conformation"], as_index=False).agg(
        s_r=(score_col, "max"),
        n_rows=(score_col, "size"),
    )
    rows = []
    for (lig, pid), g in conf_df.groupby(["Ligand_id", "protein_id"], sort=False):
        idx = g["s_r"].astype(float).idxmax()
        best_conf = g.loc[idx, "conformation"]
        best_sr = float(g.loc[idx, "s_r"])
        n_pass = int((g["s_r"].astype(float) >= cover_t).sum())
        coverage = n_pass / float(total_confs) if total_confs > 0 else 0.0
        rows.append({
            "Ligand_id": lig,
            "protein_id": pid,
            "best_conformation": best_conf,
            "protein_score_raw": best_sr * (1.0 + alpha * coverage),
            "coverage": coverage,
            "n_pass_conformations": n_pass,
            "total_conformations": int(total_confs),
            "max_s_r": best_sr,
            "mean_s_r": float(np.mean(g["s_r"])),
        })
    prot_df = pd.DataFrame(rows)
    if not prot_df.empty:
        v = prot_df["protein_score_raw"].astype(float)
        vmin = float(v.min())
        vmax = float(v.max())
        prot_df["protein_score_norm"] = 100.0 if vmax <= vmin else 100.0 * (v - vmin) / (vmax - vmin)
        prot_df = prot_df.sort_values(["protein_score_norm", "max_s_r"], ascending=False)
    return conf_df, prot_df
