#!/usr/bin/env python3
"""
Plot read-intensity features across the IGH locus for a single sample.

Panels (top→bottom):
  1) BAF (B-allele frequency)
  2) Confidence (log scale)
  3) Log2 ratio of probe intensity
  4) Genotype concordance labels (0 = match, 1 = mismatch)

Gene midpoints (from BED) are shown as vertical dotted ticks with labels,
and structural-variant regions are shaded.

Usage:
  python plot_read_intensity.py \
    --input /path/to/3510493.txt \
    --bed /path/to/igh.bed \
    --output figures/read_intensity_3510493.png \
    --title "Sample 3510493 (IGH)"

Input requirements:
  - TSV with columns: IGH_Pos, baf, confidence, log2r, label
    * IGH_Pos is in bp (will be converted to kb on the x-axis)
    * label is integer (0 match, 1 mismatch)
BED requirements:
  - 4 columns: chr, start, end, info
  - 'info' contains 'gene_name=...'
"""

import argparse
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # headless
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# Default SV regions for IGH (in kb). Edit or replace with --no-sv to disable.
DEFAULT_SV_REGIONS = [
    (46.481, 56.100, "lightcoral", "Deletion"),
    (157.859, 167.416, "lightcoral", "Deletion"),
    (980.092, 1027.109, "lightcoral", "Deletion"),
    (218.472, 500.135, "#FFD580", "mSV"),
    (485.425, 559.597, "#FFD580", "mSV"),
    (681.645, 733.851, "#FFD580", "mSV"),
    (214.579, 252.256, "indigo", "Complex"),
    (420.588, 424.103, "#006400", "Duplication"),
]


def load_table(path: str, label_col: str = "label") -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", skipinitialspace=True)
    missing = [c for c in ["IGH_Pos", "baf", "confidence", "log2r", label_col] if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    df.columns = df.columns.str.strip()
    df = df.dropna(subset=[label_col]).copy()
    df[label_col] = df[label_col].astype(int)

    # Coordinates in kb for plotting
    df["pos_kb"] = df["IGH_Pos"] / 1000.0
    return df


def load_gene_bed(bed_path: str) -> pd.DataFrame:
    bed = pd.read_csv(bed_path, sep="\t", header=None, names=["chr", "start", "end", "info"])
    bed["gene"] = bed["info"].astype(str).str.extract(r"gene_name=([A-Za-z0-9\-]+)").fillna("")
    bed["mid_kb"] = ((bed["start"] + bed["end"]) / 2.0) / 1000.0
    return bed


def add_sv_shading(ax, sv_regions):
    shown = set()
    for start, end, color, label in sv_regions:
        ax.axvspan(start, end, color=color, alpha=0.3, label=(label if label not in shown el
