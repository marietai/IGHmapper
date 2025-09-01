📦 Installation

Clone the repository and install dependencies:

git clone https://github.com/<your-username>/igh-toolkit.git
cd igh-toolkit
python -m venv .venv && source .venv/bin/activate  # Windows: .venv\Scripts\activate
pip install -r requirements.txt


Dependencies include:

pandas, numpy, scikit-learn, seaborn, matplotlib

biopython (for FASTA GC content calculation)

joblib (model persistence)

xgboost (optional; if training with --models XGBoost)

mlflow (optional; if using --mlflow logging)

⚙️ Workflow Overview

Feature engineering (add_features.py)
Create SNP-level features from base genotype calls, IGH gene annotations, SV regions, and reference FASTA.

Model training (train_model.py)
Train classifiers (RandomForest, XGBoost, DecisionTree) on labeled SNPs.
Select thresholds that prioritize precision (PPV) for “accurate” calls.

Model application (run_model.py)
Apply the trained model to new SNP data.
Choose between a unique global threshold or stratified thresholds (higher cutoff in structurally complex IGH regions).

🚀 Scripts
1. add_features.py

Generate engineered features for SNPs in the IGH locus.

Inputs:

--input TSV with base columns: igh_sample, log2r, baf, copy.

--genes IGH gene annotation CSV with start,end,gene,uniqueness_score.

--fasta FASTA sequence covering IGH locus.

Optional: --sv-json custom JSON with SV regions.

Outputs:

TSV with additional features:

SV overlap flags (sv_del, sv_dup, sv_msv, sv_complex, sv_count)

Distance features (true_dist_to_sv, log_dist_to_sv)

Absolute deviations (abs_log2r, abs_baf_diff, abs_copy_diff)

Gene context (in_igh_gene, igh_family, igh_uniqueness_score, near_igh_gene)

GC content features (gc_content, abs_gc_diff)

Example:

python add_features.py \
  --input data/training_features1.tsv \
  --genes data/IGH_Gene_Annotation_with_Uniqueness_Score.csv \
  --fasta data/igh_only.fasta \
  --output data/training_features.tsv

2. train_model.py

Train classifiers to predict SNP accuracy labels.

Inputs:

--train Labeled training TSV (must include label column: 0=accurate, 1=error).

--test Labeled test TSV.

--drop-cols List of non-predictive columns to exclude.

Outputs (in --outdir):

Trained model(s) (*_model.joblib)

Threshold metadata (*_threshold.json)

ROC and PR curves

Feature importance plots (for tree models)

(Optional) MLflow logs if --mlflow enabled.

Example:

python train_model.py \
  --train data/training_features.127labeled.tsv \
  --test  data/testing_features.48labeled.tsv \
  --outdir model_11 \
  --experiment IGH_11 \
  --models RandomForest XGBoost \
  --target-ppv 0.99 \
  --mlflow

3. run_model.py

Apply trained models to new SNPs and generate predictions.

Inputs:

--model Trained .joblib pipeline.

--threshold-json Threshold file from training.

--input TSV with feature columns.

--drop-cols Columns to exclude (same as training).

Outputs:

TSV with:

p_good: probability SNP is accurate (label 0)

predicted_label: binary prediction (1=accurate, 0=inaccurate)

Thresholding modes:

Unique: apply global threshold from JSON.

Stratified (--use-stratified):

HIGH stratum: (sv_del==0 OR sv_complex==0 OR igh_family in {1,3,4})
→ apply --thr-high (default 0.99)

ELSE stratum: apply JSON threshold or --thr-else.

Example (unique):

python run_model.py \
  --model model_11/RandomForest_model.joblib \
  --threshold-json model_11/RandomForest_threshold.json \
  --input testing_sets/testing_features_allsamples.tsv \
  --output model_11/predictions_unique_thr.tsv


Example (stratified):

python run_model.py \
  --model model_11/RandomForest_model.joblib \
  --threshold-json model_11/RandomForest_threshold.json \
  --input testing_sets/testing_features_allsamples.tsv \
  --output model_11/predictions_stratified_thr.tsv \
  --use-stratified --thr-high 0.99
