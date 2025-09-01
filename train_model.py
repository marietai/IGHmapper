#!/usr/bin/env python3
"""
train_model.py — Train genotype-accuracy classifiers for IGH features.

What it does
------------
- Loads train/test TSVs (expects a `label` column: 0=accurate, 1=error).
- Encodes categorical features consistently across train/test.
- Trains one or more classifiers with a SimpleImputer (median).
- Selects an operating threshold on TEST to prioritize precision (PPV) for the positive class (= label 0).
- Saves: fitted model (.joblib), threshold (.json), ROC/PR curves (.png), feature importance (.png, when available).
- Optionally logs metrics/artifacts to MLflow.

Examples
--------
python train_model.py \
  --train data/training_features.127labeled.tsv \
  --test  data/testing_features.48labeled.tsv \
  --outdir model_11 \
  --experiment IGH_11 \
  --models RandomForest XGBoost \
  --target-ppv 0.99 \
  --mlflow

Notes
-----
- "Positive" is defined as label==0 (accurate) to match your downstream usage.
- Dropping of non-predictive columns is configurable via --drop-cols.
"""

from __future__ import annotations
import argparse
import json
import os
from datetime import datetime
from pathlib import Path
import warnings

import joblib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.impute import SimpleImputer
from sklearn.metrics import (
    average_precision_score,
    confusion_matrix,
    precision_recall_curve,
    precision_score,
    recall_score,
    roc_auc_score,
    roc_curve,
)
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import LabelEncoder
from sklearn.ensemble import RandomForestClassifier
from sklearn.tree import DecisionTreeClassifier

try:
    from xgboost import XGBClassifier
    _HAS_XGB = True
except Exception:
    _HAS_XGB = False

# ----------------------------
# Pretty names for feature importance plot
# ----------------------------
FEATURE_NAME_MAP = {
    "IGH_Pos": "IGH custom reference coordinate",
    "gc_content": "GC content (%)",
    "baf": "BAF",
    "log2r": "Log2 ratio",
    "sv_del": "SNP overlaps deletion SV",
    "sv_dup": "SNP overlaps duplication SV",
    "sv_msv": "SNP overlaps multi-allelic SV",
    "sv_complex": "SNP overlaps complex SV",
    "sv_count": "SV overlap count",
    "in_igh_gene": "SNP within IGH gene",
    "igh_family": "IGH gene family",
    "igh_uniqueness_score": "IGH gene positional redundancy score",
    "near_igh_gene": "SNP near IGH gene (±100 bp)",
    "confidence": "Confidence score",
    "true_dist_to_sv": "Distance to nearest SV (kb)",
}

# ----------------------------
# CLI
# ----------------------------
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Train genotype-accuracy models for IGH features.")
    p.add_argument("--train", required=True, help="Training TSV (must include 'label').")
    p.add_argument("--test",  required=True, help="Test TSV (must include 'label').")
    p.add_argument("--outdir", default="model_11", help="Directory for saved artifacts.")
    p.add_argument("--experiment", default="IGH_11", help="(Optional) MLflow experiment name.")
    p.add_argument("--mlflow", action="store_true", help="If set, log metrics/artifacts to MLflow.")
    p.add_argument("--target-ppv", type=float, default=0.99, help="Target PPV for threshold selection on TEST.")
    p.add_argument(
        "--drop-cols",
        nargs="*",
        default=["igh_sample","ref","alt","gt","position","log_dist_to_sv","abs_log2r",
                 "abs_baf_diff","abs_gc_diff","lr_REF","lr_ALT","lr_GT","label"],
        help="Columns to drop from features."
    )
    p.add_argument("--random-state", type=int, default=42, help="Random seed.")
    p.add_argument("--models", nargs="+", default=["RandomForest","XGBoost"],
                   choices=["RandomForest","XGBoost","DecisionTree"],
                   help="Which models to train.")
    p.add_argument("--rf-estimators", type=int, default=500, help="n_estimators for RandomForest.")
    p.add_argument("--xgb-estimators", type=int, default=500, help="n_estimators for XGBoost.")
    return p.parse_args()

# ----------------------------
# Data utils
# ----------------------------
def load_split(train_path: Path, test_path: Path, drop_cols: list[str]) -> tuple[pd.DataFrame, pd.Series, pd.DataFrame, pd.Series]:
    df_tr = pd.read_csv(train_path, sep="\t")
    df_te = pd.read_csv(test_path,  sep="\t")

    if "label" not in df_tr.columns or "label" not in df_te.columns:
        raise ValueError("Both train and test must contain a 'label' column (0=accurate, 1=error).")

    y_tr = df_tr["label"].astype(int)
    y_te = df_te["label"].astype(int)

    # Do NOT drop rows globally; imputer will handle missing values.
    X_tr = df_tr.drop(columns=[c for c in drop_cols if c in df_tr.columns], errors="ignore")
    X_te = df_te.drop(columns=[c for c in drop_cols if c in df_te.columns], errors="ignore")

    # Ensure same columns/order across train/test
    common_cols = sorted(set(X_tr.columns) & set(X_te.columns))
    X_tr = X_tr[common_cols].copy()
    X_te = X_te[common_cols].copy()

    # Consistent encoding for any object/categorical columns
    cat_cols = list(X_tr.select_dtypes(include=["object","category"]).columns)
    for col in cat_cols:
        le = LabelEncoder()
        all_vals = pd.concat([X_tr[col].astype(str), X_te[col].astype(str)], axis=0)
        le.fit(all_vals)
        X_tr[col] = le.transform(X_tr[col].astype(str))
        X_te[col] = le.transform(X_te[col].astype(str))

    return X_tr, y_tr, X_te, y_te

# ----------------------------
# Models
# ----------------------------
def build_model(name: str, args: argparse.Namespace):
    if name == "RandomForest":
        clf = RandomForestClassifier(
            n_estimators=args.rf_estimators,
            random_state=args.random_state,
            n_jobs=-1,
        )
    elif name == "DecisionTree":
        clf = DecisionTreeClassifier(random_state=args.random_state)
    elif name == "XGBoost":
        if not _HAS_XGB:
            raise RuntimeError("XGBoost is not installed; pip install xgboost or remove from --models.")
        clf = XGBClassifier(
            n_estimators=args.xgb_estimators,
            random_state=args.random_state,
            eval_metric="logloss",
            tree_method="hist",
            use_label_encoder=False,
        )
    else:
        raise ValueError(f"Unknown model: {name}")

    # Impute missing numeric/categorical with median (simple & robust for trees)
    return make_pipeline(SimpleImputer(strategy="median"), clf)

# ----------------------------
# Thresholding where positive = label 0
# ----------------------------
def pick_threshold_for_target_ppv(y_pos: np.ndarray, proba_pos: np.ndarray, target_ppv: float) -> tuple[float,int]:
    """Return (chosen_threshold, op_idx) maximizing recall among points with precision >= target_ppv.
       If none meet target, fall back to best F0.5."""
    prec, rec, thr = precision_recall_curve(y_pos, proba_pos)  # prec/rec: n+1; thr: n
    prec_trim, rec_trim = prec[1:], rec[1:]                    # align with thr
    mask = prec_trim >= target_ppv
    if mask.any():
        idx_local = np.argmax(rec_trim[mask])
        thr_choice = float(thr[mask][idx_local])
        op_idx = int(np.where(mask)[0][idx_local])
        return thr_choice, op_idx
    # fallback: F0.5
    beta = 0.5
    fbeta = (1 + beta**2) * (prec_trim * rec_trim) / (beta**2 * prec_trim + rec_trim + 1e-12)
    op_idx = int(np.nanargmax(fbeta))
    return float(thr[op_idx]), op_idx

# ----------------------------
# Plots
# ----------------------------
def plot_roc_pr(
    outdir: Path, name: str,
    y_pos_tr: np.ndarray, proba_tr: np.ndarray,
    y_pos_te: np.ndarray, proba_te: np.ndarray,
    chosen_thr: float, op_idx: int
):
    # ROC
    fpr_tr, tpr_tr, _ = roc_curve(y_pos_tr, proba_tr)
    fpr_te, tpr_te, _ = roc_curve(y_pos_te, proba_te)
    auc_tr = roc_auc_score(y_pos_tr, proba_tr)
    auc_te = roc_auc_score(y_pos_te, proba_te)

    plt.figure()
    plt.plot(fpr_tr, tpr_tr, label=f"Train AUC = {auc_tr:.2f}", color="black")
    plt.plot(fpr_te, tpr_te, label=f"Test  AUC = {auc_te:.2f}", color="gray")
    plt.plot([0,1],[0,1],"--", color="lightgray")
    plt.xlabel("False positive rate"); plt.ylabel("True positive rate")
    plt.legend(loc="lower right", frameon=False)
    plt.tight_layout()
    plt.savefig(outdir / f"{name}_ROC_curve_train_test.png", dpi=300)
    plt.close()

    # PR
    prec_tr, rec_tr, _ = precision_recall_curve(y_pos_tr, proba_tr)
    prec_te, rec_te, thr_te = precision_recall_curve(y_pos_te, proba_te)

    ap_tr = average_precision_score(y_pos_tr, proba_tr)
    ap_te = average_precision_score(y_pos_te, proba_te)

    plt.figure()
    plt.plot(rec_tr, prec_tr, label=f"Train AP = {ap_tr:.2f}", color="black")
    plt.plot(rec_te, prec_te, label=f"Test  AP = {ap_te:.2f}", color="gray")
    # Operating point (remember op_idx corresponds to prec/rec index 1+ for PR curve)
    plt.scatter(rec_te[op_idx + 1], prec_te[op_idx + 1], s=35)
    plt.text(rec_te[op_idx + 1], prec_te[op_idx + 1],
             f"thr={chosen_thr:.2f}", fontsize=8, ha="right", va="top")
    plt.xlabel("Recall"); plt.ylabel("Precision")
    plt.legend(loc="lower left", frameon=False)
    plt.tight_layout()
    plt.savefig(outdir / f"{name}_PR_curve_train_test.png", dpi=300)
    plt.close()

def plot_feature_importance(outdir: Path, name: str, pipeline, feature_names: list[str]):
    try:
        final_est = pipeline.named_steps[list(pipeline.named_steps.keys())[-1]]
        if hasattr(final_est, "feature_importances_"):
            importances = np.asarray(final_est.feature_importances_)
            labels = np.array([FEATURE_NAME_MAP.get(c, c) for c in feature_names])
            order = np.argsort(importances)[::-1][:20]

            plt.figure(figsize=(6, 5))
            sns.barplot(x=importances[order], y=labels[order], palette=["black"] + ["gray"] * 19)
            plt.xlabel("Feature importance")
            plt.tight_layout()
            plt.savefig(outdir / f"{name}_feature_importance.png", dpi=300)
            plt.close()
    except Exception as e:
        print(f"[WARN] Could not plot feature importance for {name}: {e}")

# ----------------------------
# Main
# ----------------------------
def main():
    args = parse_args()
    warnings.filterwarnings("ignore")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Optional MLflow
    if args.mlflow:
        import mlflow
        import mlflow.sklearn
        mlflow.set_experiment(args.experiment)

    # Load data
    X_tr, y_tr, X_te, y_te = load_split(Path(args.train), Path(args.test), args.drop_cols)
    feature_cols = list(X_tr.columns)

    # Positive class is label==0
    y_pos_tr = (y_tr == 0).astype(int).to_numpy()
    y_pos_te = (y_te == 0).astype(int).to_numpy()

    # Train requested models
    for name in args.models:
        pipeline = build_model(name, args)

        if args.mlflow:
            mlflow.end_run() if "mlflow" in globals() else None  # be safe if nested
            run_ctx = mlflow.start_run(run_name=f"{name}_Split")
        else:
            run_ctx = None

        try:
            # Fit
            pipeline.fit(X_tr, y_tr)

            # Save model
            model_path = outdir / f"{name}_model.joblib"
            joblib.dump(pipeline, model_path)
            print(f"✅ Saved model to: {model_path}")

            # Probabilities for positive class (label==0)
            # Pipeline supports predict_proba; classes_ on final estimator
            classes_ = pipeline.named_steps[list(pipeline.named_steps.keys())[-1]].classes_
            pos_idx = int(np.where(classes_ == 0)[0][0])

            proba_tr = pipeline.predict_proba(X_tr)[:, pos_idx]
            proba_te = pipeline.predict_proba(X_te)[:, pos_idx]

            # Threshold selection on TEST
            chosen_thr, op_idx = pick_threshold_for_target_ppv(y_pos_te, proba_te, args.target_ppv)

            # Thresholded predictions (TEST)
            preds_te = (proba_te >= chosen_thr).astype(int)
            ppv_te = precision_score(y_pos_te, preds_te, zero_division=0)
            rec_te = recall_score(y_pos_te, preds_te, zero_division=0)

            # Persist threshold info
            thr_record = {
                "model_name": name,
                "target_ppv": float(args.target_ppv),
                "chosen_threshold": float(chosen_thr),
                "pos_label": 0,                     # positive = label 0 (accurate)
                "class_index_used": int(pos_idx),   # column in predict_proba
                "timestamp": datetime.now().isoformat(),
            }
            thr_path = outdir / f"{name}_threshold.json"
            thr_path.write_text(json.dumps(thr_record, indent=2))
            print(f"💾 Saved threshold to: {thr_path}")

            # Metrics
            ap_tr   = average_precision_score(y_pos_tr, proba_tr)
            roc_tr  = roc_auc_score(y_pos_tr, proba_tr)
            ap_te   = average_precision_score(y_pos_te, proba_te)
            roc_te  = roc_auc_score(y_pos_te, proba_te)

            tn, fp, fn, tp = confusion_matrix(y_pos_te, preds_te, labels=[0, 1]).ravel()

            # Plots
            plot_roc_pr(outdir, name, y_pos_tr, proba_tr, y_pos_te, proba_te, chosen_thr, op_idx)
            plot_feature_importance(outdir, name, pipeline, feature_cols)

            # MLflow logging (optional)
            if args.mlflow:
                mlflow.log_artifact(str(model_path))
                mlflow.log_artifact(str(thr_path))
                mlflow.log_metrics({
                    "pr_auc_pos0_train":  ap_tr,
                    "roc_auc_pos0_train": roc_tr,
                    "pr_auc_pos0_test":   ap_te,
                    "roc_auc_pos0_test":  roc_te,
                    "precision_at_thr_pos0_test": ppv_te,
                    "recall_at_thr_pos0_test":    rec_te,
                    "chosen_threshold_pos0_test": chosen_thr,
                    "TP_pos0_test": int(tp), "FP_pos0_test": int(fp),
                    "TN_pos0_test": int(tn), "FN_pos0_test": int(fn),
                })
        finally:
            if run_ctx is not None:
                import mlflow
                mlflow.end_run()

if __name__ == "__main__":
    main()

