import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import os

# === File paths ===
file_path = "/mnt/c/Users/marietai/Desktop/GoalPost/training_from_snp/3510493.txt"
gene_bed_path = "/mnt/c/Users/marietai/Desktop/GoalPost/LongRead/igh.bed"
outpath = "/mnt/c/Users/marietai/Desktop/GoalPost/training_from_snp/read_intensity_3510493.png"
label_col = 'label'

# === Load data ===
df = pd.read_csv(file_path, sep="\t", skipinitialspace=True)
df.columns = df.columns.str.strip()
df = df.dropna(subset=[label_col])  # use a list here
df[label_col] = df[label_col].astype(int)
positions = df["IGH_Pos"]
pos_kb = df["IGH_Pos"] / 1000

# === Load gene annotation ===
gene_df = pd.read_csv(gene_bed_path, sep="\t", header=None, names=["chr", "start", "end", "info"])
gene_df["gene"] = gene_df["info"].str.extract(r'gene_name=([A-Za-z0-9\-]+)').fillna("")
gene_df["mid"] = ((gene_df["start"] + gene_df["end"]) / 2) / 1000  # mid in kb

# === SV annotation (start, end, color, label), in kb ===
sv_regions = [
    (46.481, 56.100, 'lightcoral', 'Deletion'),
    (157.859, 167.416, 'lightcoral', 'Deletion'),
    (980.092, 1027.109, 'lightcoral', 'Deletion'),
    (218.472, 500.135, '#FFD580', 'mSV'),
    (485.425, 559.597, '#FFD580', 'mSV'),
    (681.645, 733.851, '#FFD580', 'mSV'),
    (214.579, 252.256, 'indigo', 'Complex'),
    (420.588, 424.103, "#006400", 'Duplication')
]

# === Create figure with custom heights ===
height_ratios = [2, 2, 2, 1]
fig = plt.figure(figsize=(36, 10))
gs = gridspec.GridSpec(nrows=4, ncols=1, height_ratios=height_ratios, hspace=0.05)
axs = [fig.add_subplot(gs[i]) for i in range(4)]
shown_labels = set()
sharex=True

# === SV shading helper ===
def add_sv_shading(ax):
    for start, end, color, label in sv_regions:
        if label not in shown_labels:
            ax.axvspan(start, end, color=color, alpha=0.3, label=label)
            shown_labels.add(label)
        else:
            ax.axvspan(start, end, color=color, alpha=0.3)

# === Plot 1: BAF ===
axs[0].scatter(pos_kb, df["baf"], s=5, alpha=0.7, color="purple")
axs[0].set_ylabel("B allele frequency\n(BAF)", fontsize=13)
axs[0].grid(True, linestyle='--', linewidth=0.5, alpha=0.6)
axs[0].tick_params(labelbottom=False)
add_sv_shading(axs[0])

# === Plot 2: Confidence (log scale) ===
axs[1].plot(pos_kb, df["confidence"], color="#707070")
axs[1].set_yscale("log")
axs[1].set_ylabel("Confidence\n(log scale)", fontsize=13)
axs[1].grid(True, linestyle='--', linewidth=0.2, alpha=0.6)
axs[1].tick_params(labelbottom=False)
add_sv_shading(axs[1])

# === Plot 3: Log2 Ratio ===
axs[2].plot(pos_kb, df["log2r"], color="red")
axs[2].set_ylabel("Log₂ ratio of\nprobe intensity", fontsize=13)
axs[2].grid(True, linestyle='--', linewidth=0.2, alpha=0.6)
axs[2].tick_params(labelbottom=False)
add_sv_shading(axs[2])

# === Plot 4: Labels (0 = Match, 1 = Mismatch) ===
axs[3].scatter(pos_kb[df["label"] == 1], [1] * sum(df["label"] == 1),
               color="#FF8000FF", s=12, label="Mismatch", zorder=3)
axs[3].scatter(pos_kb[df["label"] == 0], [0] * sum(df["label"] == 0),
               color="#1f76b4aa", s=7, label="Match", zorder=2)
axs[3].set_ylabel("Genotype\nconcordance", fontsize=13)
axs[3].set_xlabel("IGH locus position (kb)", labelpad=10, fontsize=18)
axs[3].tick_params(axis='x', pad=5, length=30, width=0.2)
axs[3].grid(True, linestyle='--', linewidth=0.5, alpha=0.6)
add_sv_shading(axs[3])

# === Add gene annotations to the bottom axis ===
for ax in axs:  # loop over axs[0], axs[1], axs[2], axs[3]
    for _, row in gene_df.iterrows():
        if 0 <= row["mid"] <= 1200:
            ax.axvline(row["mid"], color="gray", linestyle=":", linewidth=0.3, alpha=0.5)
            ax.text(
                row["mid"], -0.05, row["gene"],
                ha='center', va='top',
                rotation=90, fontsize=5,
                clip_on=False,
                transform=ax.get_xaxis_transform()
            )

# === Final styling & save ===
for ax in axs[:-1]:
    ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False)

plt.tight_layout()
plt.savefig(outpath, dpi=300)
plt.close()
print(f"✅ Saved: {outpath}")