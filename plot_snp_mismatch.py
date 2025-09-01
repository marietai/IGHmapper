#!/usr/bin/env python3
"""
Plot SNP mismatch rate across the IGH locus.

This script:
1. Loads a labeled SNP dataset and gene BED file.
2. Computes per-position mismatch rates.
3. Generates a plot with shaded SV regions and gene annotations.

Usage:
    python plot_snp_mismatch.py --input data.tsv --bed genes.bed --output mismatch_plot.png
"""
import argparse
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def load_data(input_path: str, label_col: str = "label") -> pd.DataFrame:
    """
    Load and minimally tidy the input table.

    Requirements:
    - `label` (int): 0 = match, non-zero = mismatch
    - `igh_sample` (str): begins with numeric IGH position (e.g., '12345|...')
    """
    df = pd.read_csv(input_path, sep="\t", low_memory=False)
    if label_col not in df.columns:
        raise ValueError(f"Missing required column '{label_col}' in {input_path}")
    if "igh_sample" not in df.columns:
        raise ValueError("Missing required column 'igh_sample'")

    # Drop rows missing label and cast
    df = df.dropna(subset=[label_col]).copy()
    df[label_col] = df[label_col].astype(int)

    # Extract numeric IGH position (bp) from start of 'igh_sample'
    pos = df["igh_sample"].astype(str).str.extract(r"^(\d+)")[0]
    if pos.isna().any():
        raise ValueError("Could not extract IGH position from some 'igh_sample' values. "
                         "Expected a leading number, e.g., '12345|sampleX'.")
    df["IGH_Pos"] = pos.astype(float)
    return df

def plot_mismatch(df: pd.DataFrame, gene_bed_path: str, output_path: str) -> None:
    """
    Create a scatter plot of mismatch rate vs. IGH position (kb),
    shade SV regions, and annotate gene midpoints from BED.
    """
    # Compute per-position counts
    grouped = df.groupby("IGH_Pos")["label"].value_counts().unstack().fillna(0)
    total = grouped.sum(axis=1)  # total observations at each position
    matches = grouped.get(0, pd.Series(0, index=grouped.index))
    mismatch_rate = 1.0 - (matches / total.replace(0, float("nan")))
    mismatch_rate = mismatch_rate.fillna(0.0)

    # Plot
    fig, ax = plt.subplots(figsize=(20, 5))
    ax.scatter(grouped.index / 1000, mismatch_rate, s=12, alpha=0.4,
               label="Mismatch rate")

    # === Example SV shaded regions for IGH (in kb) ===
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
    shown = set()
    for start, end, color, label in sv_regions:
        ax.axvspan(start, end, color=color, alpha=0.3,
                   label=(label if label not in shown else None))
        shown.add(label)

    # === Gene BED annotation track ===
    if gene_bed_path:
        try:
            gene_df = pd.read_csv(gene_bed_path, sep="\t", header=None,
                                  names=["chr", "start", "end", "info"])
            # Extract gene name (edit regex if your BED differs)
            gene_df["gene"] = gene_df["info"].astype(str).str.extract(r'gene_name=([A-Za-z0-9\-]+)').fillna("")
            gene_df["mid_kb"] = ((gene_df["start"] + gene_df["end"]) / 2.0) / 1000.0

            # Thin labels for readability: show every Nth gene name
            show_every = max(int(len(gene_df) / 150), 1)
            for idx, row in gene_df.iterrows():
                if 0 <= row["mid_kb"] <= 1200:
                    ax.axvline(x=row["mid_kb"], color="gray", linestyle=":", linewidth=0.3, alpha=0.6, zorder=1)
                    if idx % show_every == 0 and row["gene"]:
                        ax.text(row["mid_kb"], -0.01, row["gene"], ha='center', va='top',
                                rotation=90, fontsize=7, transform=ax.get_xaxis_transform(),
                                clip_on=False, zorder=2)
        except Exception as e:
            print(f"⚠️ BED track error: {e}")

    # Final formatting
    ax.set_xlim(-10, 1210)
    ax.set_ylim(0, 1.0)
    ax.set_xlabel("IGH locus position (kb)", fontsize=15)
    ax.set_ylabel("Mismatch rate (1 - match proportion)")
    # Legend with unique labels
    handles, labels = ax.get_legend_handles_labels()
    unique = dict(zip(labels, handles))
    ax.legend(unique.values(), unique.keys(), loc="upper right", fontsize=9)

    plt.tight_layout()
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=500, bbox_inches="tight")
    plt.close()
    print(f"✅ Plot saved to: {output_path}")

def main():
    parser = argparse.ArgumentParser(description="Plot SNP mismatch rate across the IGH locus.")
    parser.add_argument("--input", required=True, help="Input TSV with columns: 'label' and 'igh_sample'")
    parser.add_argument("--bed", required=True, help="Gene BED file (expects 'gene_name=...' in 4th column)")
    parser.add_argument("--output", default="figures/mismatch_rate_across_IGH.png",
                        help="Output image path (PNG)")
    args = parser.parse_args()

    # Resolve paths
    in_path = Path(args.input)
    bed_path = Path(args.bed)
    out_path = Path(args.output)

    if not in_path.exists():
        raise FileNotFoundError(f"Input file not found: {in_path}")
    if not bed_path.exists():
        raise FileNotFoundError(f"BED file not found: {bed_path}")

    df = load_data(str(in_path))
    plot_mismatch(df, str(bed_path), str(out_path))

if __name__ == "__main__":
    main()

