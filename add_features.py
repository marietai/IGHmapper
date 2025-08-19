import pandas as pd
import numpy as np
import re
from Bio import SeqIO

# === File paths
input_path = "/mnt/c/Users/marietai/Desktop/GoalPost/ML/data/training_features1.tsv"
igh_gene_file = "/mnt/c/Users/marietai/Desktop/GoalPost/training_sets/genotyping_out/IGH_Gene_Annotation_with_Uniqueness_Score.csv"
fasta_path = "/mnt/c/Users/marietai/Desktop/GoalPost/training_sets/igh_only.fasta"
output_file = "/mnt/c/Users/marietai/Desktop/GoalPost/ML/data/training_features.tsv"

# === Load base dataset and extract IGH position
df = pd.read_csv(input_path, sep="\t")
df["IGH_Pos"] = df["igh_sample"].astype(str).str.extract(r"^(\d+)").astype(float)

# === SV distance features
sv_regions = [
    (46481, 56100), (157859, 167416), (980092, 1027109),(420588, 424103),
    (218472, 500135), (485425, 559597), (681645, 733851), (214579, 252256)
]

# === Add SV type-specific binary features
def in_any(pos, ranges):
    return any(start <= pos <= end for start, end in ranges)

sv_del_regions = [(46481, 56100), (157859, 167416), (980092, 1027109)]
sv_dup_regions = [(420588, 424103)]
sv_msv_regions = [(218472, 500135), (485425, 559597), (681645, 733851)]
sv_complex_regions = [(214579, 252256)]

df["sv_del"] = df["IGH_Pos"].apply(lambda pos: 0 if in_any(pos, sv_del_regions) else 1)
df["sv_dup"] = df["IGH_Pos"].apply(lambda pos: 0 if in_any(pos, sv_dup_regions) else 1)
df["sv_msv"] = df["IGH_Pos"].apply(lambda pos: 0 if in_any(pos, sv_msv_regions) else 1)
df["sv_complex"] = df["IGH_Pos"].apply(lambda pos: 0 if in_any(pos, sv_complex_regions) else 1)

# === Compute sv_count: how many SV categories had a 0 (i.e. overlap)
df["sv_count"] = (
    (df["sv_del"] == 0).astype(int) +
    (df["sv_dup"] == 0).astype(int) +
    (df["sv_msv"] == 0).astype(int) +
    (df["sv_complex"] == 0).astype(int)
)

def compute_true_dist(pos, regions):
    distances = [0 if start <= pos <= end else min(abs(pos - start), abs(pos - end)) for start, end in regions]
    return min(distances) / 1000  # in kilobases

df["true_dist_to_sv"] = df["IGH_Pos"].apply(lambda pos: compute_true_dist(pos, sv_regions))
df["log_dist_to_sv"] = np.log1p(df["true_dist_to_sv"])

# === Absolute deviation features
df["abs_log2r"] = df["log2r"].abs()
df["abs_baf_diff"] = df["baf"].apply(lambda x: abs(x - 0.75) if x >= 0.5 else abs(x - 0.25))
df["abs_copy_diff"] = (df["copy"] - 2.0).abs()

# === IGH gene annotation
genes = pd.read_csv(igh_gene_file)
in_igh_gene, igh_family, igh_uniqueness_score = [], [], []

for pos in df["IGH_Pos"]:
    if pd.isna(pos):
        in_igh_gene.append(1)
        igh_family.append(17)
        igh_uniqueness_score.append(0)
        continue

    match = genes[(genes["start"] <= pos) & (genes["end"] >= pos)]
    if match.empty:
        in_igh_gene.append(1)
        igh_family.append(17)
        igh_uniqueness_score.append(0)
    else:
        in_igh_gene.append(0)
        gene = match.iloc[0]["gene"]
        score = match.iloc[0]["uniqueness_score"]
        igh_uniqueness_score.append(score)

        if isinstance(gene, str) and gene.startswith("IGHV"):
            m = re.findall(r"IGHV(\d+)", gene)
            igh_family.append(min(int(m[0]), 8) if m else 17)
        elif isinstance(gene, str) and gene.startswith("IGHD"):
            m = re.findall(r"IGHD(\d+)", gene)
            igh_family.append(min(int(m[0]) + 8, 15) if m else 17)
        elif isinstance(gene, str) and gene.startswith("IGHJ"):
            igh_family.append(16)
        else:
            igh_family.append(17)

df["in_igh_gene"] = in_igh_gene
df["igh_family"] = igh_family
df["igh_uniqueness_score"] = igh_uniqueness_score

# === Near IGH flag
def is_near_igh(pos, starts, ends, window=100):
    for start, end in zip(starts, ends):
        if abs(pos - start) <= window or abs(pos - end) <= window:
            return 0
    return 1

df["near_igh_gene"] = df["IGH_Pos"].apply(lambda pos: is_near_igh(pos, genes["start"], genes["end"]))

# === Compute GC content
record = next(SeqIO.parse(fasta_path, "fasta"))
window = 1000
gc_windows = []

for i in range(0, len(record.seq) - window + 1, window):
    win_seq = record.seq[i:i+window].upper()
    gc = (win_seq.count("G") + win_seq.count("C")) / window * 100
    gc_windows.append((i + 1, i + window, gc))

gc_df = pd.DataFrame(gc_windows, columns=["start", "end", "gc"])

def map_gc(pos):
    match = gc_df[(gc_df["start"] <= pos) & (gc_df["end"] >= pos)]
    return match["gc"].values[0] if not match.empty else np.nan

df["gc_content"] = df["IGH_Pos"].apply(map_gc)
df["abs_gc_diff"] = (df["gc_content"] - 50.0).abs()

# === Save final output
df.to_csv(output_file, sep="\t", index=False)
print(f"✅ Features saved to:\n{output_file}")
