Antibody Analysis Workflow

This folder contains a lightweight workflow for antibody (multiplex serology) analysis over the IGH locus:

Preprocessing – rank-normalize antibody titres

Association – run antibody ~ SNP tests with QC and multiple testing correction

Genome-wide style plots – per-antibody Manhattan-style plots across IGH

Top-hit inspection – per-antibody panels of titre vs genotype for top SNPs

0. Setup

Install dependencies:

pandas
numpy
scipy
statsmodels
matplotlib


(If your project uses a global requirements.txt, ensure these are included.)

1) Rank-normalize antibody titres

rank_normalization.py applies a rank-based inverse normal transform (INT) to a set of antibody columns.

Usage

python rank_normalization.py \
  --input data/Multiplex_Serology_GoalPost_Data.xlsx \
  --output results/Multiplex_Normalized.tsv
# Optional: specify explicit columns
# --cols Antibody1 Antibody2 Antibody3


Output

Adds *_ranknorm columns to the table and writes a TSV.

2) Antibody ~ SNP association with QC

antibody_association_analysis.py builds a genotype matrix from a long file,
merges with phenotypes, runs OLS per antibody×SNP with covariates, applies QC, and corrects for multiple testing.

Inputs

--geno: long genotype predictions (must include igh_sample, gt, and ideally ref/alt for labeling)

--pheno: phenotype TSV (must include UID and *_ranknorm columns from step 1)

Covariates default: Final_disease_status, Plate, Sex (customizable)

QC defaults

call rate ≥ 0.95

MAF ≥ 0.05

N ≥ 5

Usage

python antibody_association_analysis.py \
  --geno results/model_11_prediction_99.tsv \
  --pheno results/Multiplex_Normalized.tsv \
  --outdir results/phenotypic_analysis \
  --min-callrate 0.95 --min-maf 0.05 --min-n 5


Outputs (under --outdir)

assoc_results.tsv

assoc_results_with_fdr.tsv (adds FDR/BH & Bonferroni columns)

top_hits_fdr.tsv, top_hits_bonferroni.tsv, top_hits_p1e4.tsv

Per-sample exports for top hits:

top_hits_with_ref_alt.tsv (summary)

top_hits_ref_alt_map.tsv (SNP→REF/ALT)

top_hits_genotypes_wide.tsv

top_hits_genotypes_long.tsv

Tip: You can optionally filter to predicted “accurate” SNPs with --filter-predicted accurate if your genotype file contains a predicted_label column produced by the machine-learning pipeline.

3) Per-antibody IGH Manhattan-style plots

plot_antibody_snp_association.py makes a figure per antibody across IGH with optional SV shading and gene ticks.

Inputs

--assoc: assoc_results_with_fdr.tsv (must include antibody, snp_pos, pvalue; significant_fdr optional)

--bed: optional BED with gene_name= info for IGH gene ticks/labels

Usage

python plot_antibody_snp_association.py \
  --assoc results/phenotypic_analysis/assoc_results_with_fdr.tsv \
  --bed   data/igh.bed \
  --outdir results/phenotypic_analysis/per_antibody_plots \
  --p-thresh 0.001
# add --no-sv to disable SV shading


Outputs

One PNG per antibody (consistent y-axis across antibodies)

legend_only.png with the plotted symbol/region legend.

4) Top-N SNP per-antibody titre–genotype panels

plot_top_SNP.py draws side-by-side panels showing titre vs genotype (0/0, 0/1, 1/1) for the top-N SNPs per antibody.

Default uses raw titres; pass --use-norm to plot *_ranknorm instead.

Titles include β, p, genotype counts (n0/n1/n2), and (in raw mode) descriptive % changes.

Inputs

--assoc: assoc_results_with_fdr.tsv

--pheno: phenotype table (UID, plus raw titles and/or *_ranknorm)

--genotypes: long predictions (igh_sample, gt, optionally ref/alt)

--outdir: output directory

Usage (raw titres)

python plot_top_SNP.py \
  --assoc results/phenotypic_analysis/assoc_results_with_fdr.tsv \
  --pheno data/Multiplex_With_Sex.tsv \
  --genotypes model_11/model_11_prediction_stratified_thr.tsv \
  --outdir results/phenotypic_analysis/top3 \
  --top-n 3


Usage (rank-normalized titres)

python plot_top_SNP.py \
  --assoc results/phenotypic_analysis/assoc_results_with_fdr.tsv \
  --pheno results/Multiplex_Normalized.tsv \
  --genotypes model_11/model_11_prediction_stratified_thr.tsv \
  --outdir results/phenotypic_analysis/top3_norm \
  --use-norm --top-n 3


Output

<Antibody>_topN.png per antibody in --outdir.

Notes & tips

The workflow expects sample IDs to be embedded at the end of igh_sample and positions at the start (e.g., "105860505|...|UID123"). If your format differs, adjust the regex in the scripts.

If your phenotype file uses different column names for status/plate/sex, change them via --status-col, --plate-col, --sex-col.

To use only “high-confidence” SNPs, filter by predicted_label upstream (or use the --filter-predicted flag in the association script).
