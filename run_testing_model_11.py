import pandas as pd
import numpy as np
import joblib
import os
import json

# === File paths ===
model_path = r"/mnt/c/Users/marietai/Desktop/GoalPost/ML/model_11/RandomForest_model.joblib"
test_path  = r"/mnt/c/Users/marietai/Desktop/GoalPost/testing_sets/testing_features_allsamples.tsv"
output_path = r"/mnt/c/Users/marietai/Desktop/GoalPost/ML/model_11_prediction_99.tsv"  # <- 99%
thr_path = r"/mnt/c/Users/marietai/Desktop/GoalPost/ML/model_11/RandomForest_threshold.json"

with open(thr_path, "r") as f:
    thr_info = json.load(f)

# === Load model ===
model = joblib.load(model_path)

# === Load test data ===
df_test = pd.read_csv(test_path, sep="\t").dropna()

# === Columns to drop (mirror training as closely as possible) ===
drop_cols = [
    "igh_sample","ref","alt","gt","position",
    "abs_baf_diff","abs_gc_diff","abs_log2r","log_dist_to_sv",
    # add any other non-features you dropped in training (e.g., 'label', 'lr_*') if present:
    "label","lr_REF","lr_ALT","lr_GT"
]
drop_cols = [c for c in drop_cols if c in df_test.columns]
X = df_test.drop(columns=drop_cols)

# IMPORTANT: Do NOT fit new encoders here.
# Your model pipeline (SimpleImputer + estimator) expects the same numeric features you trained on.
# If you had categorical encodings during training, load/reuse those encoders instead of fitting new ones.

# === Align columns to the model's expected feature order if available ===
# Many sklearn models store feature_names_in_ when fitted with a DataFrame.
feature_names = None
if hasattr(model, "feature_names_in_"):
    feature_names = list(model.feature_names_in_)
else:
    # Try last step in a pipeline
    try:
        last = model.named_steps[list(model.named_steps.keys())[-1]]
        if hasattr(last, "feature_names_in_"):
            feature_names = list(last.feature_names_in_)
    except Exception:
        pass

if feature_names is not None:
    # Keep only expected columns, in order; add missing as NaN
    missing = [c for c in feature_names if c not in X.columns]
    if missing:
        # create missing columns with NaN (imputer in the pipeline will handle them)
        for c in missing:
            X[c] = np.nan
    X = X[feature_names]
else:
    # Fallback: ensure stable column order
    X = X.reindex(sorted(X.columns), axis=1)

# === Probability of POSITIVE class = label 0 ===
# Be robust to class ordering
def get_pos_index(clf):
    if hasattr(clf, "classes_"):
        return int(np.where(clf.classes_ == 0)[0][0])
    # Pipeline: look at final estimator
    if hasattr(clf, "named_steps"):
        last = clf.named_steps[list(clf.named_steps.keys())[-1]]
        return int(np.where(last.classes_ == 0)[0][0])
    raise RuntimeError("Could not locate classes_ on model.")

pos_idx = get_pos_index(model)
proba_pos = model.predict_proba(X)[:, pos_idx]  # P(label==0)

# Sanity: ensure we use the same proba column as in training
def get_pos_index_from_model_or_file(clf, fallback_idx):
    # Prefer the saved index; fall back to model introspection if needed
    if fallback_idx is not None:
        return int(fallback_idx)
    if hasattr(clf, "classes_"):
        return int(np.where(clf.classes_ == 0)[0][0])
    if hasattr(clf, "named_steps"):
        last = clf.named_steps[list(clf.named_steps.keys())[-1]]
        return int(np.where(last.classes_ == 0)[0][0])
    raise RuntimeError("Could not locate classes_ on model.")

pos_idx = get_pos_index_from_model_or_file(model, thr_info.get("class_index_used"))

proba_pos = model.predict_proba(X)[:, pos_idx]  # P(label==0)

THRESH = float(thr_info["chosen_threshold"])
preds_pos = (proba_pos >= THRESH).astype(int)   # 1 = predicted accurate (label 0)

print(f"Using threshold from training: THRESH={THRESH:.4f} (target PPV={thr_info.get('target_ppv')})")

# === Save output ===
out = df_test.copy()
out["probability"] = proba_pos
out["predicted_label"] = preds_pos  # 1 = predicted accurate (label 0), 0 = predicted inaccurate (label 1)
os.makedirs(os.path.dirname(output_path), exist_ok=True)
out.to_csv(output_path, sep="\t", index=False)
print(f"✅ Predictions saved to: {output_path}")

# Optional: quick summary
n_pos = int(preds_pos.sum())
print(f"Predicted 'accurate' (label 0) at THRESH={THRESH}: {n_pos} / {len(out)} ({n_pos/len(out):.2%})")
