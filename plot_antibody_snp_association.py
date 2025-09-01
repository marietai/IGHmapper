#!/usr/bin/env python3
"""
plot_antibody_snp_association.py — Per-antibody Manhattan-style plots over IGH

Inputs
------
--assoc    TSV from antibody_association_analysis.py (must include columns:
           antibody, snp_pos, pvalue, significant_fdr [bool or {0,1}])
--bed      (optional) IGH BED with 'gene_name=' in 4th col for gene ticks
--outdir   Output directory for figures

What it does
------------
- Computes -log10(p) (safe for p=0).
- Plots per-antibody points at snp_pos/1000 (kb), coloring FDR-significant vs others.
- Adds optional IGH gene vertical ticks/labels from BED.
- Adds optional SV shaded spans (IGH defaults provided, can disable).
- Writes one PNG per antibody plus a standalone legend PNG.

Examples
--------
python plot_antibody_snp_association.py \
  --assoc results/phenotypic_analysis/assoc_results_with_fdr.tsv \
  --bed   data/igh.bed \
  --outdir results/phenotypic_analysis/per_antibody_plots \
  --p-thresh 0.001 --no-sv
"""

from __future__ import annotations
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D


DEFAULT_SV_REGIONS = [
    (46.481, 56.100, 'lightcoral', 'Deletion'),
    (157.859, 167.416, 'lightcoral', 'Deletion'),
    (980.092, 1027.109, 'lightcoral', 'Deletion'),
    (218.472, 500.135, '#FFD580', 'Multi-allelic SV'),
    (485.425, 559.597, '#FFD580', 'Multi-allelic SV'),
    (681.645, 733.851, '#FFD580', 'Multi-allelic SV'),
    (214.579, 252.256, 'indigo', 'Complex'),
    (420.588, 424.103, '#006400', 'Duplication'),
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Per-antibody SNP association plots over IGH.")
    p.add_argument("--assoc", required=True, help="assoc_results_with_fdr.tsv")
    p.add_argument("--bed", default=None, help="IGH BED (optional; shows gene ticks/labels)")
    p.add_argument("--outdir", required=True, help="Output directory for per-antibody PNGs")
    p.add_argument("--p-thresh", type=float, default=0.001, help="Horizontal p-value threshold line (default 1e-3)")
    p.add_argument("--dpi", type=int, default=300, help="Figure DPI (default 300)")
    p.add_argument("--no-sv", action="store_true", help="Disable SV shaded regions")
    p.add_argument("--label-top", type=int, default=1, help="Annotate top N SNPs per antibody (default 1)")
    p.add_argument("--y-pad", type=float, default=1.5, help="Extra headroom added to global ymax")
    p.add_argument("--gene-fontsize", type=int, default=5, help="Font size for gene labels")
    return p.parse_args()


def load_assoc(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    if "antibody" not in df.columns or "snp_pos" not in df.columns or "pvalue" not in df.columns:
        raise ValueError("Association file must contain columns: antibody, snp_pos, pvalue")
    # ensure numeric
    df["snp_pos"] = pd.to_numeric(df["snp_pos"], errors="coerce").astype("Int64")
    # p->-log10(p) (guard p==0)
    df["pvalue"] = df["pvalue"].replace(0, 1e-300)
    df["minus_log10_p"] = -np.log10(df["pvalue"])
    # boolean FDR flag if present
    if "significant_fdr" in df.columns:
        df["significant_fdr"] = df["significant_fdr"].astype(int).astype(bool)
    else:
        df["significant_fdr"] = False
    return df.dropna(subset=["snp_pos", "pvalue"])


def load_bed(path: Path) -> pd.DataFrame | None:
    if path is None:
        return None
    try:
        gene_df = pd.read_csv(path, sep="\t", header=None, names=["chr", "start", "end", "info"])
        gene_df["gene"] = gene_df["info"].astype(str).str.extract(r"gene_name=([A-Za-z0-9\-]+)").fillna("")
        gene_df["mid_kb"] = ((gene_df["start"] + gene_df["end"]) / 2.0) / 1000.0
        return gene_df
    except Exception as e:
        print(f"⚠️ BED track error: {e}")
        return None


def make_legend(outdir: Path, dpi: int):
    legend_elements = [
        Line2D([0], [0], marker='o', color='w',
               label='FDR < 0.05', markerfacecolor="#1f77b4", markersize=6),
        Line2D([0], [0], marker='o', color='w',
               label='Not FDR-significant', markerfacecolor='lightblue', markersize=6),
        Line2D([0], [0], linestyle='--', color='black', label=r'$p$ threshold'),
        Patch(facecolor='lightcoral', edgecolor='none', alpha=0.2, label='Deletion'),
        Patch(facecolor='#FFD580', edgecolor='none', alpha=0.2, label='Multi-allelic SV'),
        Patch(facecolor='indigo', edgecolor='none', alpha=0.2, label='Complex'),
        Patch(facecolor="#006400", edgecolor='none', alpha=0.2, label='Duplication'),
    ]
    fig = plt.figure(figsize=(4, 1.8))
    fig.legend(handles=legend_elements, loc='center', frameon=False, ncol=1, fontsize=8)
    out = outdir / "legend_only.png"
    fig.savefig(out, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print(f"🗂  Saved legend: {out}")


def plot_per_antibody(df: pd.DataFrame, gene_df: pd.DataFrame | None, outdir: Path,
                      dpi: int, p_thresh: float, show_sv: bool, label_top: int,
                      y_pad: float, gene_fs: int):
    antibodies = df["antibody"].dropna().unique().tolist()
    global_ymax = float(df["minus_log10_p"].max()) if len(df) else 5.0

    for ab in antibodies:
        sub = df[df["antibody"] == ab].copy()
        if sub.empty:
            continue
        sub["pos_kb"] = sub["snp_pos"].astype(float) / 1000.0

        fig, ax = plt.subplots(figsize=(15, 5))

        # SV shading
        if show_sv:
            seen = set()
            for start, end, color, label in DEFAULT_SV_REGIONS:
                if label in seen:
                    lab = None
                else:
                    lab = label
                    seen.add(label)
                ax.axvspan(start, end, color=color, alpha=0.2, label=lab)

        # Points
        sig = sub[sub["significant_fdr"]]
        nonsig = sub[~sub["significant_fdr"]]
        ax.scatter(nonsig["pos_kb"], nonsig["minus_log10_p"], s=12, color="lightblue", alpha=0.8)
        ax.scatter(sig["pos_kb"],    sig["minus_log10_p"],    s=14, color="#1f77b4", alpha=0.9)

        # p-threshold
        ax.axhline(-np.log10(p_thresh), color="black", linestyle="--", linewidth=1)

        # Gene ticks/labels
        if gene_df is not None:
            for _, row in gene_df.iterrows():
                ax.axvline(row["mid_kb"], color="gray", linestyle=":", linewidth=0.3, alpha=0.5)
                if row["gene"]:
                    ax.text(row["mid_kb"], -0.05, row["gene"], ha='center', va='top',
                            rotation=90, fontsize=gene_fs,
                            clip_on=False, transform=ax.get_xaxis_transform())

        # Annotate top N SNPs by smallest p
        topN = sub.nsmallest(label_top, "pvalue")
        for _, r in topN.iterrows():
            ax.annotate(
                f"SNP {int(r['snp_pos'])} bp\n(p={r['pvalue']:.1e})",
                xy=(r["pos_kb"], r["minus_log10_p"]),
                xytext=(r["pos_kb"], r["minus_log10_p"] + 1),
                arrowprops=dict(facecolor='black', arrowstyle='->', lw=0.5),
                fontsize=8, ha='center'
            )

        # Axes & style
        ab_title = str(ab).split("_")[0]
        ax.set_title(ab_title)
        ax.set_xlabel("IGH locus position (kb)", labelpad=10)
        ax.set_ylabel(r"$-\log_{10}(p)$")
        ax.set_ylim(bottom=0, top=global_ymax + y_pad)
        ax.tick_params(axis='x', pad=5, length=30, width=0.2)
        plt.tight_layout()

        out_file = outdir / f"{ab}_snp_assoc.png"
        fig.savefig(out_file, dpi=dpi)
        plt.close(fig)
        print(f"🖼  Saved: {out_file}")


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    assoc = load_assoc(Path(args.assoc))
    gene_df = load_bed(Path(args.bed)) if args.bed else None

    plot_per_antibody(
        assoc, gene_df, outdir,
        dpi=args.dpi, p_thresh=args.p_thresh,
        show_sv=(not args.no_sv),
        label_top=args.label_top,
        y_pad=args.y_pad,
        gene_fs=args.gene_fontsize
    )

    make_legend(outdir, dpi=args.dpi)
    print(f"✅ Finished. Plots saved to: {outdir}")


if __name__ == "__main__":
    main()
