import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# =============================
# Basic settings
# =============================
INPUT_FILE = "protein_ligand_matrix_YewHook.xlsx"
SHEET_NAME = "Sheet1"

# Binding thresholds
strong_threshold = -7  # strong binding: e ≤ -7
weak_threshold   = -3  # weak binding: -7 < e ≤ -3

# Weights
alpha = 0.2  # weight for total weak-binding counts
beta  = 0.2  # weight for weak-binding ratio

# Output file names
ALL_CSV_OUT   = "cluster_binding_statistics_with_weights_alpha_0.2.csv"
ALL_PDF_OUT   = "all_clusters_with_weights_alpha_0.2.pdf"
TOP20_PDF_OUT = "top20_clusters_with_weights_alpha_0.2.pdf"
TOP20_CSV_OUT = "top20_clusters_with_weights_alpha_0.2.csv"

# =============================
# Read data
# =============================
df = pd.read_excel(INPUT_FILE, sheet_name=SHEET_NAME, header=0)

# Extract columns (assume col1=name, col2=cluster; columns from 3rd onward are substrates with binding energies)
protein_names = df['name']
clusters = df['cluster']
substrate_cols = df.columns[2:]  # columns starting from the 3rd column are substrate affinity energies

# =============================
# Statistics per cluster
# =============================
cluster_stats = []
grouped = df.groupby('cluster')

for clust, group in grouped:
    n_proteins = len(group)

    protein_stats = []
    for _, row in group.iterrows():
        data = row[substrate_cols].values  # binding energies of current protein against all substrates
        total_count  = data.size
        strong_count = (data <= strong_threshold).sum()
        weak_count   = ((data > strong_threshold) & (data <= weak_threshold)).sum()
        none_count   = (data > weak_threshold).sum()

        # Ratios
        strong_ratio = strong_count / total_count if total_count else 0.0
        weak_ratio   = weak_count   / total_count if total_count else 0.0
        none_ratio   = none_count   / total_count if total_count else 0.0

        protein_stats.append((strong_ratio, weak_ratio, none_ratio,
                              strong_count, weak_count, total_count))

    # Convert to array for aggregation
    protein_stats = np.array(protein_stats)
    avg_strong_ratio = protein_stats[:, 0].mean() if len(protein_stats) else 0.0
    avg_weak_ratio   = protein_stats[:, 1].mean() if len(protein_stats) else 0.0
    avg_none_ratio   = protein_stats[:, 2].mean() if len(protein_stats) else 0.0

    total_strong_count   = protein_stats[:, 3].sum() if len(protein_stats) else 0
    total_weak_count     = protein_stats[:, 4].sum() if len(protein_stats) else 0
    total_binding_count  = protein_stats[:, 5].sum() if len(protein_stats) else 0

    # Composite score (including weak-binding contributions)
    composite_score = (
        total_strong_count +
        alpha * total_weak_count +
        n_proteins * (avg_strong_ratio + beta * avg_weak_ratio)
    )

    cluster_stats.append({
        'cluster': clust,
        'avg_strong_ratio': avg_strong_ratio,
        'avg_weak_ratio':   avg_weak_ratio,
        'avg_none_ratio':   avg_none_ratio,
        'total_strong_count':  total_strong_count,
        'total_weak_count':    total_weak_count,
        'total_binding_count': total_binding_count,
        'composite_score':     composite_score,
        'total_proteins':      n_proteins
    })

# Create summary DataFrame and sort
cluster_df = pd.DataFrame(cluster_stats)
cluster_df = cluster_df.sort_values('composite_score', ascending=False)

# =============================
# Save full summary table
# =============================
cluster_df.to_csv(ALL_CSV_OUT, index=False, encoding='utf-8-sig')

# =============================
# Plot: all clusters
# =============================
plt.figure(figsize=(12, 8), dpi=300)
sns.barplot(data=cluster_df, x='cluster', y='composite_score', palette='Blues_r')
plt.xlabel("Cluster")
plt.ylabel("Composite Score")
plt.title("All Clusters by Composite Score (With Weights, α=0.2)")
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig(ALL_PDF_OUT, format="pdf")
plt.close()

# =============================
# Top 20: save + plot
# =============================
top_20_clusters = cluster_df.head(20).copy()
# Insert ranking column
top_20_clusters.insert(0, 'rank', np.arange(1, len(top_20_clusters) + 1))

# Keep only key columns (remove this list if you want all columns)
cols_to_save = [
    'rank', 'cluster', 'composite_score',
    'avg_strong_ratio', 'avg_weak_ratio', 'avg_none_ratio',
    'total_strong_count', 'total_weak_count', 'total_binding_count', 'total_proteins'
]
top_20_clusters.loc[:, cols_to_save].to_csv(TOP20_CSV_OUT, index=False, encoding='utf-8-sig')

# Plot Top 20 clusters
plt.figure(figsize=(10, 6), dpi=300)
sns.barplot(data=top_20_clusters, x='cluster', y='composite_score', palette='Blues_r')
plt.xlabel("Cluster")
plt.ylabel("Composite Score")
plt.title("Top 20 Clusters by Composite Score (With Weights, α=0.2)")
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig(TOP20_PDF_OUT, format="pdf")
plt.close()

print(f"Saved:\n- All clusters table: {ALL_CSV_OUT}\n- Top20 table: {TOP20_CSV_OUT}\n- All clusters plot: {ALL_PDF_OUT}\n- Top20 plot: {TOP20_PDF_OUT}")
