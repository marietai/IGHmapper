import pandas as pd
import matplotlib
matplotlib.use('Agg')  # For non-interactive backends
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

# === Load the results ===
path = "/mnt/c/Users/marietai/Desktop/GoalPost/phenotypic_analysis/assoc_results_with_fdr.tsv"
df = pd.read_csv(path, sep="\t")

# === Compute -log10(p) ===
df["pvalue"] = df["pvalue"].replace(0, 1e-300)
df["-log10_p"] = -np.log10(df["pvalue"])

# === Load gene BED annotation ===
gene_bed_path = "/mnt/c/Users/marietai/Desktop/GoalPost/LongRead/igh.bed"
gene_df = None
try:
    gene_df = pd.read_csv(gene_bed_path, sep="\t", header=None, names=["chr", "start", "end", "info"])
    gene_df["gene"] = gene_df["info"].str.extract(r'gene_name=([A-Za-z0-9\-]+)').fillna("")
    gene_df["mid"] = ((gene_df["start"] + gene_df["end"]) / 2) / 1000  # kb
except Exception as e:
    print(f"⚠️ BED track error: {e}")

# === SV Regions (not used in plot directly, but keep if needed later) ===
sv_regions = [
    (46.481, 56.100, 'lightcoral', 'Deletion'),
    (157.859, 167.416, 'lightcoral', 'Deletion'),
    (980.092, 1027.109, 'lightcoral', 'Deletion'),
    (218.472, 500.135, '#FFD580', 'Multi-allelic SV'),
    (485.425, 559.597, '#FFD580', 'Multi-allelic SV'),
    (681.645, 733.851, '#FFD580', 'Multi-allelic SV'),
    (214.579, 252.256, 'indigo', 'Complex'),
    (420.588, 424.103, "#006400", 'Duplication')
]

# === Create output folder ===
output_dir = r"/mnt/c/Users/marietai/Desktop/GoalPost/phenotypic_analysis/per_antibody_plots"
os.makedirs(output_dir, exist_ok=True)

# === Global Y-axis scale ===
global_ymax = df["-log10_p"].max()
antibodies = df["antibody"].unique()

# === Custom scatter legend handles for later ===
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], marker='o', color='w', label='False discovery rate (FDR) < 0.05 (baseline SNP array data)',
           markerfacecolor="#333333", markersize=6),
    Line2D([0], [0], marker='o', color='w', label='Did not reach FDR significance (baseline SNP array data)',
           markerfacecolor="#B0B0B0", markersize=6),
    Line2D([0], [0], marker='o', color='w', label='False discovery rate (FDR) < 0.05 (filtered SNP array data)',
           markerfacecolor="#1f77b4", markersize=6),
    Line2D([0], [0], marker='o', color='w', label='Did not reach FDR significance (filtered SNP array data)',
           markerfacecolor='lightblue', markersize=6),
    Line2D([0], [0], linestyle='--', color='black', label='$p$ = 0.05'),
    Patch(facecolor='lightcoral', edgecolor='none', alpha=0.2, label='Deletion'),
    Patch(facecolor='#FFD580', edgecolor='none', alpha=0.2, label='Multi-allelic SV'),
    Patch(facecolor='indigo', edgecolor='none', alpha=0.2, label='Complex'),
    Patch(facecolor="#006400", edgecolor='none', alpha=0.2, label='Duplication')
]

# === Plot per antibody ===
for ab in antibodies:
    sub = df[df["antibody"] == ab].copy()
    sub["snp_pos_kb"] = sub["snp_pos"] / 1000

    fig, ax = plt.subplots(figsize=(15, 5))
    # === Add SV region spans to the plot ===
    seen_sv_labels = set()
    for start, end, color, sv_type in sv_regions:
        ax.axvspan(start, end, color=color, alpha=0.2, label=sv_type if sv_type not in seen_sv_labels else None)
        seen_sv_labels.add(sv_type)

    # Plot significant
    sig = sub[sub["significant_fdr"]]
    ax.scatter(sig["snp_pos_kb"], sig["-log10_p"], s=14, color="#1f77b4", alpha=0.9)

    # Plot non-significant
    nonsig = sub[~sub["significant_fdr"]]
    ax.scatter(nonsig["snp_pos_kb"], nonsig["-log10_p"], s=12, color='lightblue', alpha=0.8)

    # Horizontal threshold
    ax.axhline(-np.log10(0.001), color="black", linestyle="--", linewidth=1)

    # Label top SNP
    top_snp = sub.loc[sub["pvalue"].idxmin()]
    top_x = top_snp["snp_pos"]
    ax.annotate(
        rf"SNP {int(top_x)} bp" + "\n" + rf"($p$ = {top_snp['pvalue']:.1e})",
        xy=(top_x / 1000, top_snp["-log10_p"]),
        xytext=(top_x / 1000, top_snp["-log10_p"] + 1),
        arrowprops=dict(facecolor='black', arrowstyle='->'),
        fontsize=8,
        ha='center'
    )

    # Gene annotations
    if gene_df is not None:
        for _, row in gene_df.iterrows():
            ax.axvline(row["mid"], color="gray", linestyle=":", linewidth=0.3, alpha=0.5)
            ax.text(
                row["mid"], -0.05, row["gene"],
                ha='center', va='top',
                rotation=90, fontsize=4,
                clip_on=False
            )

    ab_prefix = ab.split('_')[0]
    ax.set_title(f"{ab_prefix}")
    ax.set_xlabel("IGH locus position (kb)", labelpad=10)
    ax.tick_params(axis='x', pad=5, length=30, width=0.2)   # Moves 200, 400, 600 down
    ax.set_ylabel(r"$-\log_{10}(p)$")
    ax.set_ylim(bottom=0, top=global_ymax + 1.5)
    plt.tight_layout()

    # Save figure
    out_file = os.path.join(output_dir, f"{ab}_snp_assoc.png")
    plt.savefig(out_file, dpi=300)
    plt.close()

# === Generate separate legend ===
fig_legend = plt.figure(figsize=(4, 1.8))
fig_legend.legend(handles=legend_elements, loc='center', frameon=False, ncol=1, fontsize=8)
legend_out_path = os.path.join(output_dir, "legend_only.png")
fig_legend.savefig(legend_out_path, dpi=300, bbox_inches='tight')
plt.close()

print(f"✅ Finished. Plots saved to:\n{output_dir}")