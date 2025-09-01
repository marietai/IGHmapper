#!/usr/bin/env python3
"""
antibody_association_analysis.py — Antibody ~ SNP association with QC + multiple testing

What it does
------------
1) Loads genotype predictions (wide-by-variant after pivot) and phenotype table.
2) Parses sample IDs (UID) and IGH SNP positions from `igh_sample`, builds a genotype matrix.
3) QC per SNP×antibody subset:
   - call rate >= MIN_CALLRATE
   - MAF >= MIN_MAF
   - N >= MIN_N
4) Fits OLS:  antibody_ranknorm ~ SNP (0/1/2) + covariates [status, plate, sex]
5) Multiple-testing correction (BH/FDR and Bonferroni) and exports:
   - assoc_results.tsv / assoc_results_with_fdr.tsv
   - top_hits_* tsvs plus per-sample genotype tables (wide & long) with REF/ALT mapping

Usage
-----
python antibody_association_analysis.py \
  --geno data/model_11_prediction_99.tsv \
  --pheno data/Multiplex_With_Sex.tsv \
  --outdir results/phenotypic_analysis \
  --min-callrate 0.95 --min-maf 0.05 --min-n 5
"""

from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

# ---------- CLI ----------
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Antibody ~ SNP association with QC and FDR.")
    p.add_argument("--geno", required=True, help="Genotype predictions TSV (must include igh_sample, gt).")
    p.add_argument("--pheno", required=True, help="Phenotype TSV (must include UID and *_ranknorm antibody columns).")
    p.add_argument("--outdir", required=True, help="Output directory for result files.")
    p.add_argument("--min-callrate", type=float, default=0.95, help="Min call rate per SNP×antibody subset. Default 0.95")
    p.add_argument("--min-maf", type=float, default=0.05, help="Min MAF per subset. Default 0.05")
    p.add_argument("--min-n", type=int, default=5, help="Min N per subset. Default 5")
    p.add_argument("--status-col", default="Final_disease_status", help="Phenotype status column name.")
    p.add_argument("--plate-col", default="Plate", help="Plate/batch column name.")
    p.add_argument("--sex-col", default="Sex", help="Sex column name.")
    p.add_argument("--filter-predicted", choices=["none","accurate","inaccurate"], default="none",
                   help="Optional filter: keep only predicted_label==1 (accurate) or ==0 (inaccurate).")
    return p.parse_args()

# ---------- Helpers ----------
def compute_callrate(geno_vec: np.ndarray) -> float:
    g = np.asarray(geno_vec, dtype=float)
    return float(np.isfinite(g).sum()) / g.size if g.size else 0.0

def compute_maf(geno_vec: np.ndarray) -> float:
    g = np.asarray(geno_vec, dtype=float)
    obs = np.isfinite(g)
    if obs.sum() == 0:
        return np.nan
    p = np.nanmean(g) / 2.0  # 0/1/2 coding
    if not np.isfinite(p):
        return np.nan
    return float(min(p, 1.0 - p))

def pick_col(df: pd.DataFrame, candidates: list[str]) -> str | None:
    for c in candidates:
        if c in df.columns:
            return c
    return None

def encode_gt(gt) -> int:
    if isinstance(gt, str):
        gt = gt.strip()
        if gt == "0/0": return 0
        if gt in {"0/1","1/0"}: return 1
        if gt == "1/1": return 2
    if isinstance(gt, (int, float)) and gt in (0,1,2):
        return int(gt)
    return -1  # missing/unrecognized

# ---------- Main ----------
def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Load
    geno_raw = pd.read_csv(args.geno, sep="\t")
    pheno_df = pd.read_csv(args.pheno, sep="\t")
    pheno_df["UID"] = pheno_df["UID"].astype(str)

    # Optional filter by predicted_label
    if args.filter-predicted != "none" and "predicted_label" in geno_raw.columns:
        keep_val = 1 if args.filter-predicted == "accurate" else 0
        geno_raw = geno_raw[geno_raw["predicted_label"] == keep_val].copy()

    # REF/ALT discovery
    ref_col = pick_col(geno_raw, ["ref","REF","Ref","A1","Allele1"])
    alt_col = pick_col(geno_raw, ["alt","ALT","Alt","A2","Allele2"])

    # Parse UID and SNP position from igh_sample
    if "igh_sample" not in geno_raw.columns or "gt" not in geno_raw.columns:
        raise ValueError("Genotype TSV must contain 'igh_sample' and 'gt' columns.")
    geno_raw["UID"] = geno_raw["igh_sample"].astype(str).str.extract(r"([0-9]+)$")  # sample id at end
    geno_raw["pos"] = geno_raw["igh_sample"].astype(str).str.extract(r"^([0-9]+)")  # position at start

    # SNP annotation (pos→REF/ALT), one row per pos
    if ref_col and alt_col:
        snp_info = (
            geno_raw[["pos", ref_col, alt_col]]
            .dropna(subset=["pos"])
            .drop_duplicates(subset=["pos"])
            .rename(columns={ref_col: "REF", alt_col: "ALT"})
            .assign(pos=lambda d: d["pos"].astype(str))
        )
    else:
        snp_info = pd.DataFrame(columns=["pos","REF","ALT"])

    # Minimal genotype table
    geno_df = geno_raw[["UID","pos","gt"]].copy()
    geno_df["gt"] = geno_df["gt"].apply(encode_gt)
    geno_df = geno_df[geno_df["gt"] >= 0]
    geno_df = geno_df.drop_duplicates(subset=["UID","pos"], keep="first")

    # Wide genotype matrix
    geno_wide = geno_df.pivot(index="UID", columns="pos", values="gt").reset_index()
    merged_df = pd.merge(pheno_df, geno_wide, on="UID", how="inner")
    merged_df.columns = merged_df.columns.map(str)

    # Identify antibodies and SNP columns
    antibodies = [c for c in merged_df.columns if c.endswith("_ranknorm")]
    snp_positions = [c for c in merged_df.columns if c.isdigit()]
    covariates = [f"C({args.status_col})", f"C({args.plate_col})", f"C({args.sex_col})"]
    raw_covariates = [args.status_col, args.plate_col, args.sex_col]

    print(f"📊 Antibodies: {len(antibodies)} | SNPs: {len(snp_positions)} | Samples merged: {merged_df.shape[0]}")

    results = []
    snps_with_nonmissing = set()
    snps_passing_qc = set()
    skipped_low_n = skipped_low_call = skipped_low_maf = 0

    # Association scanning
    for ab in antibodies:
        for snp in snp_positions:
            try:
                subset = merged_df[[ab, snp] + raw_covariates].copy()
                subset = subset[subset[snp] >= 0]                      # keep called genotypes
                subset = subset.dropna(subset=[ab] + raw_covariates)   # drop missing pheno/covs

                if subset.shape[0] < args.min_n:
                    skipped_low_n += 1
                    continue

                snps_with_nonmissing.add(snp)

                callrate = compute_callrate(subset[snp].values)
                if callrate < args.min_callrate:
                    skipped_low_call += 1
                    continue

                maf = compute_maf(subset[snp].values)
                if (not np.isfinite(maf)) or (maf < args.min_maf):
                    skipped_low_maf += 1
                    continue

                snps_passing_qc.add(snp)

                # OLS: antibody ~ SNP + covariates
                if covariates:
                    formula = f"{ab} ~ Q('{snp}') + " + " + ".join(covariates)
                else:
                    formula = f"{ab} ~ Q('{snp}')"
                model = smf.ols(formula=formula, data=subset).fit()

                beta = model.params.get(f"Q('{snp}')", np.nan)
                pval = model.pvalues.get(f"Q('{snp}')", np.nan)

                results.append({
                    "antibody": ab,
                    "snp_pos": snp,
                    "beta": float(beta) if np.isfinite(beta) else np.nan,
                    "pvalue": float(pval) if np.isfinite(pval) else np.nan,
                    "n": int(subset.shape[0]),
                    "maf": float(maf),
                    "callrate": float(callrate),
                })

            except Exception as e:
                # continue on model errors to keep the sweep going
                print(f"❌ {ab} ~ {snp}: {e}")
                continue

    print(f"SNPs with non-missing calls (any ab): {len(snps_with_nonmissing)}")
    print(f"SNPs passing QC (any ab):            {len(snps_passing_qc)}")
    print(f"Skipped low N: {skipped_low_n} | low call-rate: {skipped_low_call} | low MAF: {skipped_low_maf}")

    # Save raw results
    results_df = pd.DataFrame(results)
    results_df.to_csv(outdir / "assoc_results.tsv", sep="\t", index=False)

    if results_df.empty:
        print("❌ No valid associations after QC.")
        return

    # Multiple testing correction
    reject, p_fdr, _, _ = multipletests(results_df["pvalue"].values, alpha=0.05, method="fdr_bh")
    results_df["pvalue_fdr"] = p_fdr
    results_df["significant_fdr"] = reject
    results_df["pvalue_bonferroni"] = results_df["pvalue"] * len(results_df)
    results_df["significant_bonferroni"] = results_df["pvalue_bonferroni"] < 0.05

    results_df.to_csv(outdir / "assoc_results_with_fdr.tsv", sep="\t", index=False)
    results_df[results_df["significant_fdr"]].to_csv(outdir / "top_hits_fdr.tsv", sep="\t", index=False)
    results_df[results_df["significant_bonferroni"]].to_csv(outdir / "top_hits_bonferroni.tsv", sep="\t", index=False)
    results_df[results_df["pvalue"] < 1e-4].to_csv(outdir / "top_hits_p1e4.tsv", sep="\t", index=False)

    print(f"✅ FDR significant: {int(results_df['significant_fdr'].sum())} | "
          f"Bonferroni significant: {int(results_df['significant_bonferroni'].sum())}")

    # Add REF/ALT to top hits and export per-sample genotypes
