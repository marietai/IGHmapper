import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# === File paths ===
input_path = "/mnt/c/Users/marietai/Desktop/GoalPost/ML/data/testing_features.48labeled.tsv"
gene_bed_path = "/mnt/c/Users/marietai/Desktop/GoalPost/LongRead/igh.bed"

# === Load data ===
df = pd.read_csv(input_path, sep="\t", low_memory=False)
label_col = 'label'

# === Drop rows where 'label' is missing ===
df = df.dropna(subset=[label_col])
df[label_col] = df[label_col].astype(int)

# === Extract IGH position from sample name ===
df["IGH_Pos"] = df["igh_sample"].astype(str).str.extract(r"^(\d+)").astype(float)

# === Calculate mismatch rate per position ===
grouped = df.groupby("IGH_Pos")["label"].value_counts().unstack().fillna(0)
grouped["mismatch_rate"] = (48 - grouped.get(0, 0)) / 48

# === Create plot ===
fig, ax = plt.subplots(figsize=(20, 6))
ax.scatter(grouped.index / 1000, grouped["mismatch_rate"], c="#FF8000", s=5, label="Mismatch rate")

# === SV shaded regions ===
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

shown_labels = set()
for start, end, color, label in sv_regions:
    show_label = label not in shown_labels
    ax.axvspan(start, end, color=color, alpha=0.3, label=label if show_label else None)
    shown_labels.add(label)

# === Gene BED annotation track ===
try:
    gene_df = pd.read_csv(gene_bed_path, sep="\t", header=None, names=["chr", "start", "end", "info"])
    gene_df["gene"] = gene_df["info"].str.extract(r'gene_name=([A-Za-z0-9\-]+)').fillna("")
    gene_df["mid"] = ((gene_df["start"] + gene_df["end"]) / 2)  # mid in raw bp
    gene_df["mid_kb"] = gene_df["mid"] / 1000

    for _, row in gene_df.iterrows():
        if 0 <= row["mid_kb"] <= 1200:
            # Full-height vertical line in data coordinates
            ax.axvline(
                x=row["mid_kb"],
                color="gray", linestyle=":", linewidth=0.3, alpha=0.6,
                zorder=1
            )

            # Extension below axis (into label area) in axis coordinates
            ax.plot(
                [row["mid_kb"], row["mid_kb"]],
                [0, -0.12],  # Extends downward
                transform=ax.get_xaxis_transform(),
                color="gray", linestyle=":", linewidth=0.3, alpha=0.6,
                clip_on=False,
                zorder=1
            )

            # Gene label below axis
            ax.text(
                row["mid_kb"], -0.05, row["gene"],
                ha='center', va='top',
                rotation=90, fontsize=5,
                transform=ax.get_xaxis_transform(),
                clip_on=False,
                zorder=2
            )

except Exception as e:
    print(f"⚠️ BED track error: {e}")

# === Final plot formatting ===
ax.set_xlim(-10, 1210)
ax.set_ylim(0, 1.0)
ax.set_xlabel("IGH locus position (kb)", labelpad=20, fontsize=15)
ax.xaxis.set_label_coords(0.5, -0.15)
ax.set_ylabel("Mismatch rate (1 - proportion of matches)")
handles, labels = ax.get_legend_handles_labels()
unique = dict(zip(labels, handles))
ax.legend(unique.values(), unique.keys(), loc="center left", bbox_to_anchor=(1, 0.5), fontsize=9)

plt.tight_layout()
plt.subplots_adjust(bottom=0.22)
plt.savefig("/mnt/c/Users/marietai/Desktop/GoalPost/training_from_snp/mismatch_rate_across_IGH_48.png", dpi=300)
plt.close()
