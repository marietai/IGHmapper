import os
import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
import statsmodels.formula.api as smf

def compute_callrate(geno_vec):
    g = np.asarray(geno_vec, dtype=float)
    return np.isfinite(g).sum() / g.size if g.size else 0.0

def compute_maf(geno_vec):
    g = np.asarray(geno_vec, dtype=float)
    obs = np.isfinite(g)
    if obs.sum() == 0:
        return np.nan
    p = np.nanmean(g) / 2.0  # allele freq with 0/1/2 coding
    if not np.isfinite(p):
        return np.nan
    return min(p, 1.0 - p)

# === File paths ===
geno_path = "/mnt/c/Users/marietai/Desktop/GoalPost/ML/model_11_prediction_99.tsv"
pheno_path = "/mnt/c/Users/marietai/Desktop/GoalPost/phenotypic_analysis/Multiplex_With_Sex.tsv"
output_dir = "/mnt/c/Users/marietai/Desktop/GoalPost/phenotypic_analysis"

# === Load genotype and phenotype ===
geno_raw = pd.read_csv(geno_path, sep="\t")  # keep a copy with REF/ALT
pheno_df = pd.read_csv(pheno_path, sep="\t")
pheno_df["UID"] = pheno_df["UID"].astype(str)

# QC thresholds (new)
MIN_CALLRATE = 0.95     # >=95% observed genotypes in the merged analysis subset
MIN_MAF = 0.05          # >=5% minor allele frequency in the analysis subset
MIN_N = 1        

# === Filter out predicted_label == 1 BEFORE selecting columns ===
#geno_raw = geno_raw[geno_raw["predicted_label"] == 0]

# --- Detect REF/ALT column names robustly ---
def pick_col(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    return None

ref_col = pick_col(geno_raw, ["ref","REF","Ref","A1","Allele1"])
alt_col = pick_col(geno_raw, ["alt","ALT","Alt","A2","Allele2"])

# === Parse position and sample ID from 'igh_sample' ===
geno_raw["UID"] = geno_raw["igh_sample"].str.extract(r"([0-9]+)$")  # sample ID at end
geno_raw["pos"] = geno_raw["igh_sample"].str.extract(r"^([0-9]+)")  # position at start

# Keep only needed columns for matrix
geno_df = geno_raw[["UID", "pos", "gt"]].copy()

# === Build SNP annotation (REF/ALT) once per position ===
if ref_col is not None and alt_col is not None:
    snp_info = (
        geno_raw[["pos", ref_col, alt_col]]
        .dropna(subset=["pos"])
        .drop_duplicates(subset=["pos"])
        .rename(columns={ref_col: "REF", alt_col: "ALT"})
        .assign(pos=lambda d: d["pos"].astype(str))
    )
else:
    # If REF/ALT absent, still create an empty table so later merges don't fail
    snp_info = pd.DataFrame(columns=["pos","REF","ALT"])

# === Encode GT as integer ===
def encode_gt(gt):
    if isinstance(gt, str):
        gt = gt.strip()
        if gt == "0/0":
            return 0
        elif gt in {"0/1", "1/0"}:
            return 1
        elif gt == "1/1":
            return 2
    return -1  # unrecognized or missing

geno_df["gt"] = geno_df["gt"].apply(encode_gt)

# Remove unrecognized genotypes
geno_df = geno_df[geno_df["gt"] >= 0]

# === Create genotype matrix ===
geno_df = geno_df.drop_duplicates(subset=["UID", "pos"], keep="first")
geno_wide = (
    geno_df.pivot(index="UID", columns="pos", values="gt")
    .fillna(-2).astype(int).reset_index()
)

# === Merge with phenotype ===
merged_df = pd.merge(pheno_df, geno_wide, on="UID", how="inner")
merged_df.columns = merged_df.columns.map(str)

# === Identify antibodies and SNPs ===
antibodies = [c for c in merged_df.columns if c.endswith("_ranknorm")]
snp_positions = [c for c in merged_df.columns if c.isdigit()]
non_missing_calls = (merged_df[snp_positions] >= 0).sum().sum()

covariates = ["C(Final_disease_status)", "C(Plate)", "C(Sex)"]
raw_covariates = ["Final_disease_status", "Plate", "Sex"]

print(f"📊 Starting association: {len(antibodies)} antibodies, {len(snp_positions)} SNPs")
print(f"🔎 Total merged samples: {merged_df.shape[0]}")
print(f"Total non-missing genotype calls: {non_missing_calls:,}")

# === Run association tests ===
results = []
snps_with_nonmissing = set()   # has at least one analyzable row
snps_passing_qc = set()        # passes callrate + MAF + MIN_N at least once

skipped_low_n = 0
skipped_low_callrate = 0
skipped_low_maf = 0


for ab in antibodies:
    for snp in snp_positions: 
        try:
            subset = merged_df[[ab, snp] + raw_covariates].copy()

            # keep rows with called genotype
            subset = subset[subset[snp] >= 0]

            # drop rows missing ab/covariates only
            subset = subset.dropna(subset=[ab] + raw_covariates)

            # require a minimal N (e.g., 5)
            if subset.shape[0] < 5:
                continue
            # require a minimal N (keep your current rule or use MIN_N)
            if subset.shape[0] < MIN_N:
                skipped_low_n += 1
                continue

            snps_with_nonmissing.add(snp)

            # compute call rate and MAF on the analysis subset genotypes
            callrate = compute_callrate(subset[snp].values)
            if callrate < MIN_CALLRATE:
                skipped_low_callrate += 1
                continue

            maf = compute_maf(subset[snp].values)
            if (not np.isfinite(maf)) or (maf < MIN_MAF):
                skipped_low_maf += 1
                continue

            # if we get here, this SNP passed QC (for this antibody)
            snps_passing_qc.add(snp)

            # fit model
            if covariates:
                formula = f"{ab} ~ Q('{snp}') + " + " + ".join(covariates)
            else:
                formula = f"{ab} ~ Q('{snp}')"

            model = smf.ols(formula=formula, data=subset).fit()

            beta = model.params.get(f"Q('{snp}')", float("nan"))
            pval = model.pvalues.get(f"Q('{snp}')", float("nan"))

            results.append({
                "antibody": ab,
                "snp_pos": snp,
                "beta": beta,
                "pvalue": pval,
                "n": subset.shape[0],
            })

        except Exception as e:
            print(f"❌ Error on {ab} vs {snp}: {e}")
            continue

print(f"Number of SNPs with at least one non-missing genotype in analysis: {len(snps_with_nonmissing)}")
print(f"SNPs passing call-rate (≥{MIN_CALLRATE}) and MAF (≥{MIN_MAF}) at least once: {len(snps_passing_qc)}")
print(f"  Skipped due to low N (<{MIN_N}): {skipped_low_n}")
print(f"  Skipped due to low call-rate:   {skipped_low_callrate}")
print(f"  Skipped due to low MAF:         {skipped_low_maf}")

# write lists (optional)
pd.Series(sorted(snps_passing_qc)).to_csv(
    f"{output_dir}/snps_passing_callrate_maf.tsv", index=False, header=False
)
pd.Series(sorted(snps_with_nonmissing)).to_csv(
    f"{output_dir}/snps_with_nonmissing.tsv", index=False, header=False
)

# === Save raw results ===
results_df = pd.DataFrame(results)
results_df.to_csv(f"{output_dir}/assoc_results.tsv", sep="\t", index=False)

# === Multiple testing correction ===
if results_df.empty:
    print("❌ No valid associations found.")
else:
    pvals = results_df["pvalue"].values
    reject, pvals_fdr, _, _ = multipletests(pvals, alpha=0.05, method="fdr_bh")

    results_df["pvalue_fdr"] = pvals_fdr
    results_df["significant_fdr"] = reject
    results_df["pvalue_bonferroni"] = results_df["pvalue"] * len(results_df)
    results_df["significant_bonferroni"] = results_df["pvalue_bonferroni"] < 0.05

    # === Save corrected results ===
    results_df.to_csv(f"{output_dir}/assoc_results_with_fdr.tsv", sep="\t", index=False)
    results_df[results_df["significant_fdr"]].to_csv(f"{output_dir}/top_hits_fdr.tsv", sep="\t", index=False)
    results_df[results_df["significant_bonferroni"]].to_csv(f"{output_dir}/top_hits_bonferroni.tsv", sep="\t", index=False)
    results_df[results_df["pvalue"] < 1e-4].to_csv(f"{output_dir}/top_hits_p1e4.tsv", sep="\t", index=False)

    print(f"✅ Completed. Significant (FDR): {results_df['significant_fdr'].sum()} | Bonferroni: {results_df['significant_bonferroni'].sum()}")

# === View genotypes for top variants ===
if not results_df.empty:
    # pick your "top" set
    top_hits = results_df[results_df["significant_fdr"]].copy()
    if top_hits.empty:
        top_hits = results_df.nsmallest(50, "pvalue").copy()  # fallback: top 50 by p

    # add REF/ALT to the summary of top hits
    top_hits = (
        top_hits.merge(snp_info.rename(columns={"pos":"snp_pos"}), on="snp_pos", how="left")
        .sort_values("pvalue")
    )
    top_hits.to_csv(f"{output_dir}/top_hits_with_ref_alt.tsv", sep="\t", index=False)

    # --- Per-sample genotypes for the top SNPs (WIDE) ---
    top_snp_cols = top_hits["snp_pos"].astype(str).unique().tolist()
    decode_map = {0: "0/0", 1: "0/1", 2: "1/1", -2: "./."}
    top_geno_wide = merged_df[["UID"] + raw_covariates + top_snp_cols].copy()
    for snp in top_snp_cols:
        top_geno_wide[snp] = top_geno_wide[snp].map(decode_map)

    # Save a header block with REF/ALT mapping for readability
    refalt_map = snp_info[snp_info["pos"].isin(top_snp_cols)].copy()
    refalt_map = refalt_map.rename(columns={"pos":"snp_pos"})
    refalt_map.to_csv(f"{output_dir}/top_hits_ref_alt_map.tsv", sep="\t", index=False)
    top_geno_wide.to_csv(f"{output_dir}/top_hits_genotypes_wide.tsv", sep="\t", index=False)

    # --- Per-sample genotypes for the top SNPs (LONG; includes REF/ALT inline) ---
    long = (
        top_geno_wide
        .melt(id_vars=["UID"] + raw_covariates, value_vars=top_snp_cols,
              var_name="snp_pos", value_name="genotype")
        .merge(refalt_map, on="snp_pos", how="left")
        .sort_values(["snp_pos","UID"])
    )
    long.to_csv(f"{output_dir}/top_hits_genotypes_long.tsv", sep="\t", index=False)

    print(f"🧬 Exported:")
    print(f"  • top_hits_with_ref_alt.tsv (summary with REF/ALT)")
    print(f"  • top_hits_ref_alt_map.tsv (SNP → REF/ALT lookup)")
    print(f"  • top_hits_genotypes_wide.tsv (per-sample wide table)")
    print(f"  • top_hits_genotypes_long.tsv (per-sample long table with REF/ALT)")