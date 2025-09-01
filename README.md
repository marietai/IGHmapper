IGH Toolkit

A collection of scripts for working with the Immunoglobulin Heavy Chain (IGH) locus:

Data conversion: Convert PLINK datasets from GRCh38 to a custom IGH coordinate system.

Plotting: Visualize mismatch rates, read intensities, and feature distributions.

Analysis: Explore genotype–feature relationships across the IGH locus.

Installation

Clone this repository and install Python dependencies:

git clone https://github.com/<your-username>/igh-toolkit.git
cd igh-toolkit
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt


Dependencies:

pandas, matplotlib, seaborn, statsmodels

plink2 (for dataset extraction and regeneration)

Scripts
1. Coordinate Conversion (igh_convert_and_check.py)

Manage conversion between GRCh38 and custom IGH reference coordinates.

Identify IGH locus coverage
python igh_convert_and_check.py identify-locus --map data/igh_to_hg38.txt

Extract IGH region from PLINK dataset
python igh_convert_and_check.py extract-plink \
  --bfile-prefix data/GoalPost_final_QCpassed \
  --from-bp 105860505 --to-bp 106883718 \
  --out-prefix data/GoalPost_IGH_loci

Convert .bim to IGH coordinates
python igh_convert_and_check.py convert-bim \
  --bim data/GoalPost_IGH_loci.bim \
  --map data/igh_to_hg38.txt \
  --out-bim data/GoalPost_IGH_reference.bim --overwrite-alleles

Regenerate .bed/.fam with new .bim
python igh_convert_and_check.py regen-bedfam \
  --bfile-prefix data/GoalPost_IGH_loci \
  --bim data/GoalPost_IGH_reference.bim \
  --out-prefix data/GoalPost_IGH_reference --run

Confirm mapping
python igh_convert_and_check.py confirm-mapping \
  --custom-bim data/GoalPost_IGH_reference.bim \
  --map data/igh_to_hg38.txt \
  --out-success results/successfully_mapped.txt \
  --out-fail results/unsuccessfully_mapped.txt

2. Plot SNP Mismatch (scripts/plot_snp_mismatch.py)

Plot mismatch rates across IGH with SV region shading and gene annotations.

python scripts/plot_snp_mismatch.py \
  --input data/all_features.tsv \
  --bed data/igh.bed \
  --output figures/mismatch_rate.png

3. Plot Read Intensity (scripts/plot_read_intensity.py)

Visualize BAF, confidence, log2 ratio, and genotype concordance across IGH.

python scripts/plot_read_intensity.py \
  --input data/sample_3510493.txt \
  --bed data/igh.bed \
  --output figures/read_intensity.png

4. Analyze Features vs Labels (scripts/analyze_features_grammar.py)

Explore continuous and categorical features vs genotype mismatch labels.

python scripts/analyze_features_grammar.py \
  --input data/all_features.tsv \
  --outdir results/feature_analysis


This produces:

Violin plots for continuous features

Bar plots of category proportions
