#!/usr/bin/env python3
"""
plot_top_SNP.py — For each antibody, plot the top-N SNP hits (titre vs genotype)

Inputs
------
--assoc       Assoc table from antibody_association_analysis.py
              (must include: antibody, snp_pos, beta, pvalue, pvalue_fdr)
--pheno       Phenotype table (must include: UID, raw titres and/or *_ranknorm)
--genotypes   Long genotype predictions with columns: igh_sample, gt
--outdir      Output directory for per-antibody top-N plots

Behavior
--------
- Builds per-sample genotype matrix (0/1/2) from the long file (infers UID and pos from igh_sample).
- For each antibody:
  - select up to --top-n SNPs: FDR<0.05 first (by p), then remaining by p (best-first).
  - plot raw titres (default) or rank-normalized values (--use-norm).
  - 1 row of up to N panels, each panel = scatter with jitter + mean±SEM per genotype (0/0,0/1,1/1).
  - Panel title includes beta, p, n0/n1/n2, and optional descriptive deltas on raw titres.
- Saves <Antibody>_topN.png in --outdir.

Examples
--------
# Raw titre mode (default)
python plot_top_SNP.py \
  --assoc results/per_antibody/assoc_results_with_fdr.tsv \
  --pheno data/Multiplex_With_Sex.tsv \
  --genotypes model_11/model_11_prediction_stratified_thr.tsv \
  --outdir results/per_antibody/top3 \
  --top-n 3

# Rank-normalized mode
python plot_top_SNP.py \
  --assoc results/per_antibody/assoc_results_with_fdr.tsv \
  --pheno data/Multiplex_Normalized.tsv \
  --genotypes model_11/model_11_prediction_stratified_thr.tsv \
  --outdir results/per_antibody/top3_norm \
  --use-norm --top-n 3
"""

from __future__ import annotations
import argparse
import math
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Plot per-antibody top-N SNP hits (titre vs genotype).")
    p.add_argument("--assoc", required=True, help="assoc_results_with_fdr.tsv")
    p.add_argument("--pheno", required=True, help="Phenotype table (with UID)")
    p.add_argument("--genotypes", required=True, help="Long genotype predictions (with igh_sample, gt)")
    p.add_argument("--outdir", required=True, help="Output directory for plots")
    p.add_argument("--use-norm", action="store_true", help="Use *_ranknorm columns for y-axis instead of raw titres")
    p.add_argument("--raw-cols", nargs="*", default=None,
                   help="Explicit raw titre column names to consider (default: auto-detect common antibody names if present)")
    p.add_argument("--top-n", type=int, default=3, help="Max SNPs per antibody to plot (default: 3)")
    p.add_argument("--dpi", type=int, default=500, help="Figure DPI (default: 500)")
    p.add_argument("--seed", type=int, default=1234, help="Jitter RNG seed")
    return p.parse_args()


# ----------------------- helpers -----------------------
def encode_gt(gt):
    """Map '0/0','0/1','1/1' to 0/1/2; other/NA -> -1"""
    if isinstance(gt, str):
        s = gt.strip()
        if s == "0/0": return 0
        if s in {"0/1", "1/0"}: return 1
        if s == "1/1": return 2
    if isinstance(gt, (int, float)) and gt in (0, 1, 2):
        return int(gt)
    return -1

def jitter(values, amount=0.08, rng=None):
    if rng is None:
        rng = np.random.default_rng(1234)
    return values + rng.normal(0, amount, size=len(values))

def group_summary(y):
    """Return (mean, SEM) for error bars."""
    y = np.asarray(y, dtype=float)
    y = y[np.isfinite(y)]
    if y.size == 0:
        return np.nan, np.nan
    mean = y.mean()
    sem = y.std(ddof=1) / math.sqrt(y.size) if y.size > 1 else np.nan
    return mean, sem

def top_hits_for_antibody(df_ab: pd.DataFrame, top_n: int) -> pd.DataFrame:
    """FDR-significant first (by p), then remaining by p; cap at top_n; order by position for readability."""
    df_ab = df_ab.copy()
    df_ab["snp_pos"] = df_ab["snp_pos"].astype(str)
    sig = df_ab[df_ab.get("pvalue_fdr", 1.0) < 0.05].sort_values("pvalue")
    rest = df_ab[df_ab.get("pvalue_fdr", 1.0) >= 0.05].sort_values("pvalue")
    top = pd.concat([sig, rest]).head(top_n)
    if top.empty:
        return top
    return top.sort_values("snp_pos", key=lambda x: x.astype(int))

def pick_raw_antibody_columns(pheno: pd.DataFrame) -> set[str]:
    # If user not given, try a common panel used earlier; else any non-ranknorm columns that look like titres
    candidates = {"Spy0843","MC2C3","SpeA","DNaseB","SCPA","GAC","SpyCEP","SpyAD","SLO","SpnA"}
    present = candidates & set(pheno.columns)
    if present:
        return present
    # fallback: any numeric, non-UID, not *_ranknorm
    numeric = {c for c in pheno.columns if pd.api.types.is_numeric_dtype(pheno[c])}
    return {c for c in numeric if not str(c).endswith("_ranknorm") and c != "UID"}


# ----------------------- main -----------------------
def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Load assoc + phenotype
    assoc = pd.read_csv(args.assoc, sep="\t")
    pheno = pd.read_csv(args.pheno, sep="\t")
    if "UID" not in pheno.columns:
        raise ValueError("Phenotype file must contain 'UID' column.")
    pheno["UID"] = pheno["UID"].astype(str)

    # Decide which antibody columns we will plot
    if args.use_norm:
        antibody_cols = [c for c in pheno.columns if str(c).endswith("_ranknorm")]
        y_label = "Antibody (rank-normalized)"
    else:
        antibody_cols = list(args.raw_cols) if args.raw_cols else list(pick_raw_antibody_columns(pheno))
        y_label = "Antibody (raw titre)"
    if not antibody_cols:
        raise ValueError("No antibody columns found to plot (check --use-norm or --raw-cols).")

    # Build genotype matrix from long file
    geno_long = pd.read_csv(args.genotypes, sep="\t")
    if "igh_sample" not in geno_long.columns or "gt" not in geno_long.columns:
        raise ValueError("Genotype long file must include 'igh_sample' and 'gt'.")
    geno_long["UID"] = geno_long["igh_sample"].astype(str).str.extract(r"([0-9]+)$")
    geno_long["pos"] = geno_long["igh_sample"].astype(str).str.extract(r"^([0-9]+)")
    geno_long = geno_long.dropna(subset=["UID", "pos"]).copy()
    geno_long["UID"] = geno_long["UID"].astype(str)
    geno_long["pos"] = geno_long["pos"].astype(str)
    geno_long["gt_code"] = geno_long["gt"].apply(encode_gt)
    geno_long = geno_long[geno_long["gt_code"] >= 0].drop_duplicates(subset=["UID","pos"], keep="first")

    # Optional REF/ALT lookup (if present)
    ref_alt = None
    for ref_name in ("ref","REF","A1","Allele1"):
        if ref_name in geno_long.columns:
            for alt_name in ("alt","ALT","A2","Allele2"):
                if alt_name in geno_long.columns:
                    ref_alt = (geno_long.drop_duplicates("pos")[["pos", ref_name, alt_name]]
                               .rename(columns={ref_name:"REF", alt_name:"ALT"}))
                    ref_alt["pos"] = ref_alt["pos"].astype(str)
                    break
            break

    geno_wide = geno_long.pivot(index="UID", columns="pos", values="gt_code").reset_index()
    geno_wide.columns.name = None

    merged = pd.merge(pheno, geno_wide, on="UID", how="inner")
    merged.columns = merged.columns.map(str)

    # Normalize assoc table types
    assoc["antibody"] = assoc["antibody"].astype(str)
    assoc["snp_pos"] = assoc["snp_pos"].astype(str)

    # RNG for jitter
    rng = np.random.default_rng(args.seed)

    # Iterate antibodies by assoc groups
    for ab_name, df_ab in assoc.groupby("antibody", sort=False):
        # choose phenotype column name matching mode
        pheno_col = ab_name if args.use_norm else (ab_name[:-9] if ab_name.endswith("_ranknorm") else ab_name)
        if pheno_col not in antibody_cols or pheno_col not in merged.columns:
            continue

        hits = top_hits_for_antibody(df_ab, args.top_n)
        if hits.empty:
            continue

        snps = hits["snp_pos"].astype(str).tolist()
        meta = hits.set_index("snp_pos")[["beta","pvalue"]].to_dict(orient="index")

        nplots = min(args.top_n, len(snps))
        fig, axes = plt.subplots(1, nplots, figsize=(nplots * 4.8, 4), squeeze=False)
        axes = axes[0]
        fig.suptitle(f"{pheno_col}: top {nplots} SNP hits (titre vs genotype)", fontsize=12)

        for i, snp in enumerate(snps[:nplots]):
            ax = axes[i]

            dfi = merged[["UID", pheno_col, snp]].copy()
            dfi = dfi.dropna(subset=[pheno_col])
            dfi = dfi[dfi[snp].isin([0,1,2])]

            if dfi.empty:
                ax.set_visible(False)
                continue

            # Descriptive stats on raw titres (if in raw mode)
            pct_11_vs_00 = np.nan
            per_allele_pct = np.nan
            if not args.use_norm:
                g0 = dfi.loc[dfi[snp] == 0, pheno_col].astype(float)
                g2 = dfi.loc[dfi[snp] == 2, pheno_col].astype(float)
                if g0.size and g2.size:
                    med0, med2 = float(np.median(g0)), float(np.median(g2))
                    if np.isfinite(med0) and med0 > 0 and np.isfinite(med2):
                        pct_11_vs_00 = 100.0 * (med2 - med0) / med0
                # crude log-linear slope ~ percent per allele
                try:
                    sub = dfi.loc[dfi[pheno_col] > 0, [pheno_col, snp]].dropna()
                    if sub.shape[0] >= 10:
                        X = sub[[snp]].values.astype(float)
                        y = np.log(sub[pheno_col].values.astype(float))
                        X_ = np.c_[np.ones(len(X)), X]
                        beta_hat = np.linalg.lstsq(X_, y, rcond=None)[0][1]
                        per_allele_pct = 100.0 * (np.exp(beta_hat) - 1.0)
                except Exception:
                    pass

            # Optional REF/ALT
            ref, alt = None, None
            if ref_alt is not None:
                m = ref_alt.loc[ref_alt["pos"] == snp]
                if not m.empty:
                    ref = str(m["REF"].iloc[0])
                    alt = str(m["ALT"].iloc[0])

            # Panel title
            bp = meta.get(snp, {"beta": np.nan, "pvalue": np.nan})
            n0 = int((dfi[snp] == 0).sum())
            n1 = int((dfi[snp] == 1).sum())
            n2 = int((dfi[snp] == 2).sum())
            title = f"SNP {snp}"
            if ref and alt:
                title += f" ({ref}>{alt})"
            title += f"\nβ={bp['beta']:.3g}, p={bp['pvalue']:.2e} | n0={n0}, n1={n1}, n2={n2}"

            extra = []
            if np.isfinite(pct_11_vs_00):
                extra.append(f"Δ% med(1/1 vs 0/0)={pct_11_vs_00:+.1f}%")
            if np.isfinite(per_allele_pct):
                extra.append(f"%/allele≈{per_allele_pct:+.1f}%")
            if extra:
                title += "\n" + " | ".join(extra)
            ax.set_title(title, fontsize=9)

            # Scatter with jitter + mean±SEM
            order = [0,1,2]
            labels = {0: "0/0", 1: "0/1", 2: "1/1"}
            xs, ys = [], []
            for g in order:
                vals = dfi.loc[dfi[snp] == g, pheno_col].values.astype(float)
                if vals.size:
                    xs.extend(jitter(np.full(vals.size, g), rng=rng))
                    ys.extend(vals)
            ax.scatter(xs, ys, s=16, alpha=0.3, color="#3F3F3F")

            for g in order:
                vals = dfi.loc[dfi[snp] == g, pheno_col].values.astype(float)
                mean, sem = group_summary(vals)
                if np.isfinite(mean):
                    ax.errorbar(g, mean, yerr=sem, fmt="o", capsize=2, color="black")

            ax.set_xticks(order)
            ax.set_xticklabels([labels[g] for g in order])
            ax.set_xlabel("Genotype")
            ax.set_ylabel(y_label)

        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        out_path = outdir / f"{pheno_col}_top{nplots}.png"
        plt.savefig(out_path, dpi=args.dpi)
        plt.close(fig)
        print(f"🖼  Saved: {out_path}")

    print(f"✅ Done. Plots saved in: {outdir}")


if __name__ == "__main__":
    main()
