#!/usr/bin/env python3
"""
Analyze SNP array feature distributions vs. genotype labels.

Generates:
  1. Violin plots for continuous features vs. genotype label.
  2. Bar plots showing:
     - Proportions of each category normalized within genotype class.
     - Proportions normalized within feature category.

Usage:
  python analyze_features_grammar.py \
    --input data/all_features.tsv \
    --outdir results/feature_analysis
"""

import argparse
import os
from pathlib import Path
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# === Feature definitions ===
LINEAR_FEATURES = [
    "confidence", "log2r", "abs_log2r",
    "baf", "abs_baf_diff",
    "gc_content", "abs_gc_diff",
    "IGH_Pos", "true_dist_to_sv",
    "log_dist_to_sv", "igh_uniqueness_score"
]
CATEGORICAL_FEATURES = [
    "sv_del", "sv_dup", "sv_msv", "sv_complex",
    "sv_count", "in_igh_gene", "near_igh_gene", "igh_family"
]
LABEL_COL = "label"

FEATURE_LABEL_MAP = {
    "confidence": "Confidence score (log scale)",
    "log2r": "Log₂ ratio",
    "abs_log2r": "Log₂ ratio deviation",
    "baf": "BAF",
    "abs_baf_diff": "BAF deviation",
    "gc_content": "GC content (%)",
    "abs_gc_diff": "GC content deviation",
    "true_dist_to_sv": "Distance to nearest SV (kb)",
    "log_dist_to_sv": "Distance to nearest SV (ln-scaled)",
    "igh_uniqueness_score": "IGH gene positional redundancy score",
}

VALUE_LABEL_MAP = {
    "label=0": "Genotype match",
    "label=1": "Genotype mismatch",
    "sv_del=0": "SNP in deletion",
    "sv_dup=0": "SNP in duplication",
    "sv_msv=0": "SNP in multi-allelic SV",
    "sv_complex=0": "SNP in complex SV",
    "sv_count=0": "No SV overlap",
    "sv_count=1": "Overlap with 1 SV",
    "sv_count=2": "Overlap with 2 SVs",
    "in_igh_gene=0": "SNP in IGH gene",
    "near_igh_gene=0": "SNP within 100 bp of IGH gene",
    "near_igh_gene=1": "SNP beyond 100 bp of IGH gene",
    "igh_family=1": "SNP in IGHV1",
    "igh_family=2": "SNP in IGHV2",
    "igh_family=3": "SNP in IGHV3",
    "igh_family=4": "SNP in IGHV4",
    "igh_family=5": "SNP in IGHV5",
    "igh_family=6": "SNP in IGHV6",
    "igh_family=7": "SNP in IGHV7",
    "igh_family=8": "SNP in IGHV8",
    "igh_family=9": "SNP in IGHD1",
    "igh_family=10": "SNP in IGHD2",
    "igh_family=11": "SNP in IGHD3",
    "igh_family=12": "SNP in IGHD4",
    "igh_family=13": "SNP in IGHD5",
    "igh_family=14": "SNP in IGHD6",
    "igh_family=15": "SNP in IGHD7",
    "igh_family=16": "SNP in IGHJ",
}

ORDERED_FEATURE_VALUES = [
    "sv_del=0", "sv_dup=0", "sv_msv=0", "sv_complex=0",
    "sv_count=0", "sv_count=1", "sv_count=2",
    "in_igh_gene=0",
    "near_igh_gene=0", "near_igh_gene=1",
    "igh_family=1", "igh_family=2", "igh_family=3", "igh_family=4",
    "igh_family=5", "igh_family=6", "igh_family=7", "igh_family=8",
    "igh_family=9", "igh_family=10", "igh_family=11", "igh_family=12",
    "igh_family=13", "igh_family=14", "igh_family=15", "igh_family=16",
]


def load_data(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    if LABEL_COL not in df.columns:
        raise ValueError(f"Missing '{LABEL_COL}' in {path}")
    df = df.dropna(subset=[LABEL_COL]).copy()
    df[LABEL_COL] = df[LABEL_COL].astype(int)
    return df


def plot_linear_features(df: pd.DataFrame, outdir: Path):
    for feature in LINEAR_FEATURES:
        plt.figure(figsize=(6, 4))
        sns.violinplot(
            x=LABEL_COL, y=feature, data=df, inner="box",
            palette={"0": "lightblue", "1": "lightsalmon"}
        )
        plt.xlabel("Genotype (0 = match, 1 = mismatch)", fontsize=13)
        plt.ylabel(FEATURE_LABEL_MAP.get(feature, feature), fontsize=12)
        if feature == "confidence":
            plt.yscale("log")
        plt.tight_layout()
        plt.savefig(outdir / f"{feature}_vs_label_violin.png", dpi=300)
        plt.close()


def compute_proportions(df: pd.DataFrame):
    label_counts = {lbl: (df[LABEL_COL] == lbl).sum() for lbl in [0, 1]}
    prop_data = []

    for feature in CATEGORICAL_FEATURES:
        for label_val, total in label_counts.items():
            subset = df[df[LABEL_COL] == label_val]
            counts = subset[feature].value_counts().sort_index()
            for val, cnt in counts.items():
                key = f"{feature}={val}"
                if key in VALUE_LABEL_MAP:
                    prop_data.append({
                        "feature_value": VALUE_LABEL_MAP[key],
                        "label": f"Label {label_val}",
                        "proportion": cnt / total if total else 0,
                    })
    return prop_data, label_counts


def fill_missing(prop_data: list, by: str = "label") -> list:
    # Fill missing label/feature combinations with 0
    all_keys = {f"{VALUE_LABEL_MAP[k]}|Label {lbl}" for k in ORDERED_FEATURE_VALUES for lbl in [0, 1]}
    observed = {f"{d['feature_value']}|{d['label']}" for d in prop_data}
    for key in all_keys - observed:
        name, lbl = key.split("|")
        prop_data.append({"feature_value": name, "label": lbl, "proportion": 0.0})
    return prop_data


def plot_barplots(df: pd.DataFrame, outdir: Path):
    # 1) Normalized within genotype class
    prop_data, _ = compute_proportions(df)
    prop_data = fill_missing(prop_data)
    prop_df = pd.DataFrame(prop_data)
    plt.figure(figsize=(20, 5))
    sns.barplot(
        data=prop_df,
        x="feature_value", y="proportion", hue="label",
        order=[VALUE_LABEL_MAP[k] for k in ORDERED_FEATURE_VALUES if k in VALUE_LABEL_MAP],
        palette=["lightblue", "lightsalmon"],
    )
    plt.xticks(rotation=45, ha="right", fontsize=10)
    plt.ylabel("Proportion (normalized within genotype)", fontsize=12)
    plt.xlabel("Genomic context category", fontsize=12)
    plt.legend(title="Genotype", fontsize=9)
    plt.tight_layout()
    plt.savefig(outdir / "genotype_proportions_per_label.png", dpi=300)
    plt.close()

    # 2) Normalized within feature category
    cat_data = []
    for feature in CATEGORICAL_FEATURES:
        counts = df.groupby(feature)[LABEL_COL].value_counts(normalize=True).unstack().fillna(0)
        for val in counts.index:
            for label_val in [0, 1]:
                key = f"{feature}={val}"
                if key in VALUE_LABEL_MAP:
                    cat_data.append({
                        "feature_value": VALUE_LABEL_MAP[key],
                        "label": f"Label {label_val}",
                        "proportion": counts.loc[val, label_val],
                    })
    cat_data = fill_missing(cat_data)
    cat_df = pd.DataFrame(cat_data)
    plt.figure(figsize=(20, 6))
    sns.barplot(
        data=cat_df,
        x="feature_value", y="proportion", hue="label",
        order=[VALUE_LABEL_MAP[k] for k in ORDERED_FEATURE_VALUES if k in VALUE_LABEL_MAP],
        palette=["lightblue", "lightsalmon"],
    )
    plt.xticks(rotation=45, ha="right", fontsize=9)
    plt.ylabel("Proportion (normalized within feature)", fontsize=12)
    plt.xlabel("Genomic context category", fontsize=12)
    plt.legend().remove()
    plt.tight_layout()
    plt.savefig(outdir / "genotype_distribution_per_category.png", dpi=300)
    plt.close()


def main():
    parser = argparse.ArgumentParser(description="Analyze SNP features vs. genotype labels.")
    parser.add_argument("--input", required=True, help="Input TSV with features and labels")
    parser.add_argument("--outdir", default="feature_label_relationship", help="Output directory")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = load_data(args.input)
    plot_linear_features(df, outdir)
    plot_barplots(df, outdir)

    print(f"✅ Analysis complete. Plots saved in {outdir}")


if __name__ == "__main__":
    main()
