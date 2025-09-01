#!/usr/bin/env python3
"""
add_features.py — Build engineered features for IGH SNPs.

Inputs
------
1) --input : TSV with at least columns
      - igh_sample (string; starts with numeric IGH position)
      - log2r, baf, copy  (float)
2) --genes : CSV with IGH gene annotation including columns
      - start, end, gene, uniqueness_score
3) --fasta : FASTA containing the IGH reference sequence (1-based)
4) --output: Path to write the enriched TSV

What it does
------------
- Parse IGH_Pos from 'igh_sample' (leading digits).
- Add IGH SV overlap flags (sv_del/sv_dup/sv_msv/sv_complex) and sv_count.
- Compute distance to the nearest SV (true_dist_to_sv in kb, plus log_dist_to_sv = ln(1+dist)).
- Absolute deviations:
    abs_log2r, abs_baf_diff (distance to {0.25, 0.75}), abs_copy_diff (from 2).
- IGH gene annotation:
    in_igh_gene (0=in gene, 1=outside), igh_family (V:1–8, D:9–15, J:16, other:17),
    igh_uniqueness_score (from annotation file).
- near_igh_gene (0 if within --near-window bp of gene boundary; 1 otherwise).
- GC features from FASTA in non-overlapping windows (--gc-window, default 1000 bp):
    gc_content (%), abs_gc_diff (from 50%).

SV regions
----------
Uses embedded IGH SV regions (bp). You can override with a JSON file via --sv-json,
providing a dict with keys: sv_del, sv_dup, sv_msv, sv_complex, each a list of [start, end].

Example
-------
python add_features.py \
  --input data/training_features1.tsv \
  --genes data/IGH_Gene_Annotation_with_Uniqueness_Score.csv \
  --fasta data/igh_only.fasta \
  --output data/training_features.tsv
"""

from __future__ import annotations
import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO


# ---- Default IGH SV regions (bp) ----
DEFAULT_SV = {
    "sv_del":    [(46481, 56100), (157859, 167416), (980092, 1027109)],
    "sv_dup":    [(420588, 424103)],
    "sv_msv":    [(218472, 500135), (485425, 559597), (681645, 733851)],
    "sv_complex": [(214579, 252256)],
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Add engineered features for IGH SNPs.")
    p.add_argument("--input", required=True, help="Input TSV with base features (includes igh_sample, log2r, baf, copy).")
    p.add_argument("--genes", required=True, help="IGH gene annotation CSV (start,end,gene,uniqueness_score).")
    p.add_argument("--fasta", required=True, help="FASTA of IGH reference sequence.")
    p.add_argument("--output", required=True, help="Output TSV path.")
    p.add_argument("--gc-window", type=int, default=1000, help="GC window size (bp). Default: 1000")
    p.add_argument("--near-window", type=int, default=100, help="Distance (bp) to count as 'near gene'. Default: 100")
    p.add_argument("--sv-json", default=None,
                   help="Optional JSON file with SV regions. Keys: sv_del, sv_dup, sv_msv, sv_complex; "
                        "each a list of [start, end].")
    return p.parse_args()


# ---------- Helpers ----------

def extract_igh_pos(series: pd.Series) -> pd.Series:
    """Extract leading digits from igh_sample -> IGH_Pos (float)."""
    pos = series.astype(str).str.extract(r"^(\d+)")[0]
    return pd.to_numeric(pos, errors="coerce").astype(float)


def in_any_ranges(pos: float, ranges: list[tuple[int, int]]) -> bool:
    return any(start <= pos <= end for start, end in ranges)


def min_dist_to_regions(pos: float, regions: list[tuple[int, int]]) -> float:
    """Min distance (bp) to any region, 0 if within."""
    if np.isnan(pos):
        return np.nan
    dists = []
    for (s, e) in regions:
        if s <= pos <= e:
            dists.append(0)
        else:
            dists.append(min(abs(pos - s), abs(pos - e)))
    return float(min(dists)) if dists else np.nan


def compute_sv_features(df: pd.DataFrame, sv_regions: dict) -> pd.DataFrame:
    # Binary flags: 0 if inside the SV category, 1 otherwise (matches your original encoding)
    for key, ranges in sv_regions.items():
        df[key] = df["IGH_Pos"].apply(lambda p: 0 if in_any_ranges(p, ranges) else 1)

    # sv_count = number of categories overlapped
    df["sv_count"] = sum((df[k] == 0).astype(int) for k in sv_regions.keys())

    # true_dist_to_sv (kb) and log_dist_to_sv
    all_regions = [r for v in sv_regions.values() for r in v]
    df["true_dist_to_sv"] = df["IGH_Pos"].apply(lambda p: min_dist_to_regions(p, all_regions) / 1000.0)
    df["log_dist_to_sv"] = np.log1p(df["true_dist_to_sv"])
    return df


def add_abs_deviation_features(df: pd.DataFrame) -> pd.DataFrame:
    df["abs_log2r"] = df["log2r"].abs()
    # BAF distance to {0.25, 0.75}
    df["abs_baf_diff"] = df["baf"].apply(lambda x: abs(x - 0.75) if x >= 0.5 else abs(x - 0.25))
    df["abs_copy_diff"] = (df["copy"] - 2.0).abs()
    return df


def map_gene_family(gene: str) -> int:
    """IGHV:1–8, IGHD:9–15, IGHJ:16, other/NA:17."""
    if isinstance(gene, str) and gene.startswith("IGHV"):
        m = re.findall(r"IGHV(\d+)", gene)
        return min(int(m[0]), 8) if m else 17
    if isinstance(gene, str) and gene.startswith("IGHD"):
        m = re.findall(r"IGHD(\d+)", gene)
        return min(int(m[0]) + 8, 15) if m else 17
    if isinstance(gene, str) and gene.startswith("IGHJ"):
        return 16
    return 17


def add_gene_features(df: pd.DataFrame, genes: pd.DataFrame, near_window: int) -> pd.DataFrame:
    """
    Add in_igh_gene, igh_family, igh_uniqueness_score, near_igh_gene.
    Vectorizes near_igh_gene using NumPy distances to gene boundaries.
    """
    # In-gene flags/attributes (loop is OK; genes table is modest)
    in_gene = []
    fam = []
    uniq = []

    for pos in df["IGH_Pos"].values:
        if np.isnan(pos):
            in_gene.append(1)
            fam.append(17)
            uniq.append(0)
            continue
        m = genes[(genes["start"] <= pos) & (genes["end"] >= pos)]
        if m.empty:
            in_gene.append(1)
            fam.append(17)
            uniq.append(0)
        else:
            in_gene.append(0)
            g = m.iloc[0]
            fam.append(map_gene_family(g.get("gene", "")))
            uniq.append(g.get("uniqueness_score", 0))

    df["in_igh_gene"] = in_gene
    df["igh_family"] = fam
    df["igh_uniqueness_score"] = uniq

    # near_igh_gene (0 if within window of any start/end, else 1)
    starts = genes["start"].to_numpy(dtype=float)
    ends = genes["end"].to_numpy(dtype=float)

    def near_flag(pos: float) -> int:
        if np.isnan(pos) or len(starts) == 0:
            return 1
        # distance to nearest boundary
        d = np.minimum(np.abs(pos - starts), np.abs(pos - ends)).min()
        return 0 if d <= near_window else 1

    df["near_igh_gene"] = df["IGH_Pos"].apply(near_flag)
    return df


def compute_gc_from_fasta(fasta_path: Path, window: int) -> pd.Series:
    """
    Return a pandas Series indexed by 0-based bin (pos-1)//window -> GC%.
    """
    record = next(SeqIO.parse(str(fasta_path), "fasta"))
    seq = str(record.seq).upper()
    n = len(seq)
    bins = (n + window - 1) // window

    gc_vals = np.zeros(bins, dtype=float)
    for i in range(bins):
        s = i * window
        e = min(s + window, n)
        win = seq[s:e]
        if not win:
            gc_vals[i] = np.nan
        else:
            gc_vals[i] = 100.0 * (win.count("G") + win.count("C")) / len(win)
    return pd.Series(gc_vals)


def add_gc_features(df: pd.DataFrame, fasta_path: Path, window: int) -> pd.DataFrame:
    gc_series = compute_gc_from_fasta(fasta_path, window)  # index = bin
    # For 1-based IGH_Pos: bin = floor((pos-1)/window)
    bins = np.floor((df["IGH_Pos"].values - 1) / window).astype("Int64")
    df["gc_content"] = pd.Series(bins).map(gc_series).values
    df["abs_gc_diff"] = (df["gc_content"] - 50.0).abs()
    return df


def main():
    args = parse_args()

    # Load base dataset
    df = pd.read_csv(args.input, sep="\t")
    if "igh_sample" not in df.columns:
        raise ValueError("Input must contain column 'igh_sample'.")
    for col in ("log2r", "baf", "copy"):
        if col not in df.columns:
            raise ValueError(f"Input must contain column '{col}'.")

    df["IGH_Pos"] = extract_igh_pos(df["igh_sample"])

    # SV regions: default or custom JSON
    if args.sv_json:
        sv_regions = json.loads(Path(args.sv_json).read_text())
        # Ensure tuples
        sv_regions = {k: [tuple(x) for x in v] for k, v in sv_regions.items()}
    else:
        sv_regions = DEFAULT_SV

    # Add features
    df = compute_sv_features(df, sv_regions)
    df = add_abs_deviation_features(df)

    # IGH gene annotation
    genes = pd.read_csv(args.genes)
    required_gene_cols = {"start", "end", "gene", "uniqueness_score"}
    if not required_gene_cols.issubset(genes.columns):
        missing = required_gene_cols - set(genes.columns)
        raise ValueError(f"Gene file missing required columns: {missing}")

    df = add_gene_features(df, genes, near_window=args.near_window)

    # GC features
    df = add_gc_features(df, Path(args.fasta), window=args.gc_window)

    # Save
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, sep="\t", index=False)
    print(f"✅ Features saved to: {out}")


if __name__ == "__main__":
    main()
